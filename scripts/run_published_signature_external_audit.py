"""Source-aware external audit of published HNSCC immunotherapy signatures.

This pilot deliberately does not fit a predictive model. It applies exact published
gene lists to pretreatment T-cell/CD8 pseudobulk (or whole-tumour bulk where noted),
using a common within-cohort mean-z scoring rule so direction and transportability
can be compared without information leaking from outcomes.
"""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import rdata
from scipy import sparse
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "published_signature_external_audit"
STRICT_IMCISION_EXCLUSIONS = {"Pat21", "Pat27", "Pat34"}
RNG = np.random.default_rng(20260910)


GZMK40 = [
    "GZMK", "CD27", "LAG3", "TNFRSF9", "ITM2C", "CD81", "GZMH", "AC034238.1",
    "CXCL13", "ZBED2", "TNFSF4", "NKG7", "NAB1", "SEMA4A", "TMEM155", "MTRNR2L12",
    "TSC22D1", "DUSP4", "TOX", "SH2D1A", "APOBEC3G", "VCAM1", "DGKZ", "CRTAM",
    "ITGB2", "HAVCR2", "LRRN3", "AL391056.1", "TRBV21-1", "PTMS", "ENPEP",
    "HLA-DPA1", "CD44", "SAMSN1", "AL606807.1", "BHLHE40", "IKZF3", "LINC00158",
    "CD74", "CTLA4",
]


def load_signatures() -> tuple[dict[str, list[str]], dict[str, dict[str, str]]]:
    korea = pd.read_excel(
        ROOT / "data/external/GSE286827/supplements_extracted/mmc7.xlsx", header=None
    ).iloc[1:, 0].dropna().astype(str).tolist()
    imc_dir = ROOT / "data/external/biomarker_transport/PMC10551666_supplements"
    dysf = pd.read_excel(imc_dir / "cd-22-0851_supp9.xlsx", header=3)
    dysf = dysf.loc[
        dysf["Dysfunctional trajectory signature"].astype(str).str.upper().eq("YES"),
        "gene_short_name",
    ].dropna().astype(str).tolist()
    shared_raw = pd.read_excel(imc_dir / "cd-22-0851_supp10.xlsx", header=None)
    shared = shared_raw.iloc[4:, 0].dropna().astype(str).tolist()
    signatures = {
        "GZMK40": GZMK40,
        "Cytotoxic2": ["GZMA", "PRF1"],
        "Exhausted5": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
        "IFNG6": ["IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG"],
        "TLS9": ["BCL6", "CD86", "CXCR4", "LAMP3", "SELL", "CCR7", "CXCL13", "CCL21", "CCL19"],
        "TcellInflamedGEP18": [
            "CD276", "HLA-DQA1", "CD274", "IDO1", "HLA-DRB1", "HLA-E", "CMKLR1",
            "PDCD1LG2", "PSMB10", "LAG3", "CXCL9", "STAT1", "CD8A", "CCL5", "NKG7",
            "TIGIT", "CD27", "CXCR6",
        ],
        "KoreaResponderNetwork30": korea,
        "IMCISION_DysfunctionalTrajectory238": dysf,
        "IMCISION_SharedTrajectory84": shared,
    }
    provenance = {
        "GZMK40": {"source_cohort": "GSE296954", "context": "pretreatment T-cell pseudobulk", "direction": "higher_favorable"},
        "Cytotoxic2": {"source_cohort": "external_literature", "context": "published comparator", "direction": "higher_favorable"},
        "Exhausted5": {"source_cohort": "external_literature", "context": "published comparator", "direction": "unspecified"},
        "IFNG6": {"source_cohort": "external_literature", "context": "published comparator", "direction": "higher_favorable"},
        "TLS9": {"source_cohort": "external_literature", "context": "published comparator", "direction": "higher_favorable"},
        "TcellInflamedGEP18": {"source_cohort": "external_literature", "context": "published comparator", "direction": "higher_favorable"},
        "KoreaResponderNetwork30": {"source_cohort": "GSE286827", "context": "pretreatment CXCL13+ Tex network", "direction": "higher_favorable"},
        "IMCISION_DysfunctionalTrajectory238": {"source_cohort": "GSE232240_IMCISION", "context": "dysfunctional CD8 trajectory", "direction": "higher_activation_dysfunction"},
        "IMCISION_SharedTrajectory84": {"source_cohort": "GSE232240_IMCISION", "context": "shared Treg/CD8 activation trajectory", "direction": "higher_activation"},
    }
    for name, genes in signatures.items():
        if len(genes) != len(set(genes)):
            raise ValueError(f"Duplicate genes in {name}")
    if len(korea) != 30 or len(dysf) != 238 or len(shared) != 84 or len(GZMK40) != 40:
        raise ValueError("Published signature length mismatch")
    return signatures, provenance


def bh_adjust(p: pd.Series) -> np.ndarray:
    values = p.to_numpy(float)
    order = np.argsort(values)
    ranked = values[order]
    adj = np.minimum.accumulate((ranked * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.minimum(adj, 1.0)
    return out


def matrix_market_pseudobulk(matrix_path: Path, row_to_gene: dict[int, str]) -> tuple[dict[str, float], float, int]:
    selected = {g: 0.0 for g in row_to_gene.values()}
    total = 0.0
    n_cells = 0
    with gzip.open(matrix_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            n_rows, n_cells, _ = map(int, line.split())
            break
        for line in handle:
            row, _, value = line.split()
            val = float(value)
            total += val
            gene = row_to_gene.get(int(row))
            if gene is not None:
                selected[gene] += val
    return selected, total, n_cells


def read_feature_symbols(path: Path) -> list[str]:
    frame = pd.read_csv(path, sep="\t", header=None, compression="gzip", dtype=str)
    column = 1 if frame.shape[1] > 1 else 0
    return frame.iloc[:, column].str.upper().tolist()


def finalize_expression(rows: list[dict], union: list[str]) -> pd.DataFrame:
    frame = pd.DataFrame(rows).set_index("sample_key")
    counts = frame.reindex(columns=union, fill_value=0.0).astype(float)
    scale = 1e6 / frame["total_umi"].replace(0, np.nan)
    expr = np.log1p(counts.mul(scale, axis=0))
    meta_cols = [c for c in frame.columns if c not in union]
    return frame[meta_cols].join(expr)


def read_main(union: list[str]) -> pd.DataFrame:
    raw = ROOT / "data/raw/GSE288199_decoded"
    clinical = pd.read_csv(ROOT / "results/clinical_incremental_pilot/cohort_and_features.tsv", sep="\t")
    clinical = clinical[clinical["clinical_eligible"].astype(bool)].set_index("patient_id")
    cache_path = OUT / "main_pretreatment_pseudobulk_cache.tsv"
    if cache_path.exists():
        cached = pd.read_csv(cache_path, sep="\t")
        if set(clinical.index).issubset(set(cached["patient_id"])):
            return finalize_expression(cached.to_dict("records"), union)
    rows = []
    for matrix_path in sorted(raw.glob("GSE288199_HN*_pre_T_matrix.mtx.gz")):
        match = re.search(r"_(HN\d+)_pre_", matrix_path.name)
        if not match or match.group(1) not in clinical.index:
            continue
        patient = match.group(1)
        feature_path = Path(str(matrix_path).replace("_matrix.mtx.gz", "_features.tsv.gz"))
        symbols = read_feature_symbols(feature_path)
        row_to_gene = {i + 1: g for i, g in enumerate(symbols) if g in union}
        selected, total, n_cells = matrix_market_pseudobulk(matrix_path, row_to_gene)
        row = {"sample_key": patient, "patient_id": patient, "cohort": "GSE288199",
               "regimen": clinical.loc[patient, "regimen"], "timepoint": "pre",
               "continuous_outcome": float(clinical.loc[patient, "pathologic_response_percent"]),
               "binary_outcome": int(float(clinical.loc[patient, "pathologic_response_percent"]) >= 50),
               "endpoint": "continuous_pTR", "compartment": "sorted_T_cells",
               "n_cells": n_cells, "total_umi": total, **selected}
        rows.append(row)
        pd.DataFrame(rows).to_csv(cache_path, sep="\t", index=False)
    return finalize_expression(rows, union)


def aggregate_sparse_by_group(counts: sparse.spmatrix, groups: np.ndarray, n_groups: int) -> np.ndarray:
    indicator = sparse.csr_matrix((np.ones(len(groups)), (np.arange(len(groups)), groups)), shape=(len(groups), n_groups))
    aggregated = counts @ indicator
    return aggregated.toarray().astype(float) if sparse.issparse(aggregated) else np.asarray(aggregated, dtype=float)


def read_korea(union: list[str]) -> pd.DataFrame:
    obj = rdata.read_rds(ROOT / "data/external/biomarker_transport/GSE286827_countdata_CD8T.rds.gz")
    counts = sparse.csc_matrix((np.asarray(obj.x), np.asarray(obj.i), np.asarray(obj.p)), shape=tuple(np.asarray(obj.Dim, int)))
    genes = np.asarray(obj.Dimnames[0]).astype(str)
    meta = rdata.read_rds(ROOT / "data/external/biomarker_transport/GSE286827_metadata_CD8T_TCR.rds.gz")
    pre = meta["Time"].astype(str).str.lower().eq("pre").to_numpy()
    meta = meta.iloc[np.flatnonzero(pre)].copy()
    counts = counts[:, pre]
    patients = pd.Index(meta["Patient"].astype(str).unique())
    group = patients.get_indexer(meta["Patient"].astype(str))
    totals = aggregate_sparse_by_group(counts, group, len(patients)).sum(axis=0)
    chosen = np.flatnonzero(np.isin(genes, union))
    selected = aggregate_sparse_by_group(counts[chosen, :], group, len(patients))
    lookup = {str(genes[i]): selected[j] for j, i in enumerate(chosen)}
    rows = []
    for j, patient in enumerate(patients):
        one = meta.loc[meta["Patient"].astype(str).eq(patient)].iloc[0]
        row = {"sample_key": patient, "patient_id": patient, "cohort": "GSE286827",
               "regimen": str(one["Neoadj_type"]), "timepoint": "pre",
               "continuous_outcome": float(one["Tumor_regression"]),
               "binary_outcome": int(float(one["Tumor_regression"]) >= 50),
               "endpoint": "continuous_tumor_regression", "compartment": "CD8_T_cells",
               "n_cells": int(np.sum(group == j)), "total_umi": totals[j]}
        row.update({g: float(lookup[g][j]) if g in lookup else 0.0 for g in union})
        rows.append(row)
    out = finalize_expression(rows, union)
    return out.loc[out["n_cells"].ge(50)]


def read_imcision(union: list[str]) -> pd.DataFrame:
    meta = pd.read_csv(ROOT / "data/external/state_metadata/GSM7324295_Meta_data_IMCISION.txt.gz", sep="\t")
    eligible = meta["timepoint"].astype(str).str.lower().eq("pre") & meta["cell_type"].isin(["CD4", "CD8"])
    sub = meta.loc[eligible].copy()
    patients = pd.Index(sub["patient"].astype(str).unique())
    group = patients.get_indexer(sub["patient"].astype(str))
    selected = {g: np.zeros(len(patients), dtype=float) for g in union}
    totals = np.zeros(len(patients), dtype=float)
    count_path = ROOT / "data/external/biomarker_transport/GSM7324294_Count_data_IMCISION.txt.gz"
    with gzip.open(count_path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if header != meta["cell_id"].astype(str).tolist():
            raise ValueError("IMCISION count/metadata order mismatch")
        idx = np.flatnonzero(eligible.to_numpy())
        for line in handle:
            gene, values = line.rstrip("\r\n").split("\t", 1)
            arr = np.fromstring(values, sep="\t", dtype=float)[idx]
            agg = np.bincount(group, weights=arr, minlength=len(patients))
            totals += agg
            if gene in selected:
                selected[gene] += agg
    rows = []
    for j, patient in enumerate(patients):
        one = sub.loc[sub["patient"].astype(str).eq(patient)].iloc[0]
        row = {"sample_key": patient, "patient_id": patient, "cohort": "GSE232240_IMCISION",
               "regimen": "nivolumab_ipilimumab", "timepoint": "pre",
               "continuous_outcome": np.nan, "binary_outcome": int(str(one["response"]) == "RE"),
               "endpoint": "strict_MPR_vs_pathologic_0_to_20", "compartment": "T_cells",
               "n_cells": int(np.sum(group == j)), "total_umi": totals[j]}
        row.update({g: float(selected[g][j]) for g in union})
        rows.append(row)
    out = finalize_expression(rows, union)
    return out.loc[out["n_cells"].ge(50) & ~out["patient_id"].isin(STRICT_IMCISION_EXCLUSIONS)]


def read_gse296954(union: list[str]) -> pd.DataFrame:
    source = ROOT / "data/external/biomarker_transport/GSE296954_GEX"
    pathology = {"1": 27, "2": 41, "5": 29, "8": 23, "12": 43, "14": 46,
                 "17": 68, "18": 89, "19": 60, "20": 31, "21": 27}
    rows = []
    for path in sorted(source.glob("*_filtered_feature_bc_matrix.h5")):
        token = path.name.split("_")[2]
        if not token.endswith("A"):
            continue
        patient = token[:-1]
        with h5py.File(path, "r") as handle:
            m = handle["matrix"]
            counts = sparse.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:]))
            genes = np.asarray([x.decode() for x in m["features"]["name"][:]])
        totals_cell = np.asarray(counts.sum(axis=0)).ravel()
        detected = np.diff(counts.indptr)
        mito = np.asarray(counts[np.char.startswith(genes, "MT-"), :].sum(axis=0)).ravel()
        mito_pct = np.divide(100 * mito, totals_cell, out=np.zeros_like(totals_cell, dtype=float), where=totals_cell > 0)
        qc = (totals_cell >= 500) & (detected >= 200) & (mito_pct < 25)
        tcell = np.asarray((counts[np.isin(genes, ["CD3D", "CD3E", "TRAC"]), :] > 0).sum(axis=0)).ravel() > 0
        cd8 = np.asarray((counts[genes == "CD8B", :] > 0).sum(axis=0)).ravel() > 0
        keep = qc & tcell & cd8
        chosen = np.flatnonzero(np.isin(genes, union))
        sums = np.asarray(counts[chosen, :][:, keep].sum(axis=1)).ravel()
        lookup = {str(genes[i]): float(sums[j]) for j, i in enumerate(chosen)}
        row = {"sample_key": patient, "patient_id": patient, "cohort": "GSE296954",
               "regimen": "bintrafusp_alfa", "timepoint": "pre",
               "continuous_outcome": float(pathology[patient]), "binary_outcome": int(pathology[patient] >= 50),
               "endpoint": "continuous_pathologic_regression", "compartment": "CD8_proxy",
               "n_cells": int(keep.sum()), "total_umi": float(totals_cell[keep].sum())}
        row.update({g: lookup.get(g, 0.0) for g in union})
        rows.append(row)
    out = finalize_expression(rows, union)
    return out.loc[out["n_cells"].ge(50)]


def read_gse281729(union: list[str]) -> pd.DataFrame:
    raw = pd.read_csv(ROOT / "data/external/GSE281729/GSE281729_expression_logTPM.txt.gz", sep="\t", header=None, dtype=str)
    sample_ids = raw.iloc[13, 2:].astype(str).tolist()
    expression = raw.iloc[14:, 2:].apply(pd.to_numeric, errors="coerce")
    expression.index = raw.iloc[14:, 1].astype(str).str.upper()
    expression.columns = sample_ids
    expression = expression.groupby(level=0).mean().T
    endpoint = pd.read_csv(ROOT / "results/transportability_screen/gse281729_scores.tsv", sep="\t").set_index("sample_key")
    common = expression.index.intersection(endpoint.index)
    rows = expression.loc[common].reindex(columns=union, fill_value=0.0)
    rows["patient_id"] = rows.index
    rows["cohort"] = "GSE281729_bulk"
    rows["regimen"] = "neoadjuvant_immunotherapy"
    rows["timepoint"] = "pre"
    rows["continuous_outcome"] = np.nan
    rows["binary_outcome"] = endpoint.loc[common, "OVERALL_PATH_GE20"].astype(int).to_numpy()
    rows["endpoint"] = "overall_pathologic_response_ge20"
    rows["compartment"] = "whole_tumour_bulk"
    rows["n_cells"] = np.nan
    rows["total_umi"] = np.nan
    return rows


def score_cohort(frame: pd.DataFrame, signatures: dict[str, list[str]], union: list[str]) -> tuple[pd.DataFrame, list[dict]]:
    expr = frame[union].astype(float)
    mean = expr.mean(axis=0)
    sd = expr.std(axis=0, ddof=0).replace(0, np.nan)
    z = expr.sub(mean, axis=1).div(sd, axis=1)
    output = frame.drop(columns=union).copy()
    coverage = []
    measured = set(expr.columns[(expr.max(axis=0) > 0) | (expr.min(axis=0) < 0)])
    for name, genes in signatures.items():
        present = [g for g in genes if g in measured]
        output[name] = z[present].mean(axis=1) if present else np.nan
        coverage.append({"cohort": str(frame["cohort"].iloc[0]), "signature": name,
                         "n_genes": len(genes), "n_present": len(present),
                         "coverage": len(present) / len(genes),
                         "missing_genes": ";".join(sorted(set(genes) - set(present)))})
    return output, coverage


def bootstrap_effect(x: np.ndarray, y: np.ndarray, continuous: bool, n_boot: int = 2000) -> tuple[float, float]:
    values = []
    if continuous:
        for _ in range(n_boot):
            idx = RNG.integers(0, len(x), len(x))
            if np.unique(y[idx]).size > 1:
                values.append(spearmanr(x[idx], y[idx]).statistic)
    else:
        pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
        for _ in range(n_boot):
            idx = np.r_[RNG.choice(pos, len(pos), True), RNG.choice(neg, len(neg), True)]
            values.append(roc_auc_score(y[idx], x[idx]))
    return tuple(np.nanquantile(values, [0.025, 0.975]))


def metrics(scores: pd.DataFrame, signatures: dict[str, list[str]], provenance: dict[str, dict[str, str]]) -> pd.DataFrame:
    rows = []
    for cohort, frame in scores.groupby("cohort", sort=False):
        continuous = frame["continuous_outcome"].notna().all()
        y = frame["continuous_outcome" if continuous else "binary_outcome"].to_numpy(float)
        for name in signatures:
            x = frame[name].to_numpy(float)
            keep = np.isfinite(x) & np.isfinite(y)
            xx, yy = x[keep], y[keep]
            if len(xx) < 8 or np.unique(yy).size < 2:
                continue
            if continuous:
                test = spearmanr(xx, yy)
                effect, p, effect_name = test.statistic, test.pvalue, "spearman_rho"
            else:
                effect = roc_auc_score(yy.astype(int), xx)
                p = mannwhitneyu(xx[yy == 1], xx[yy == 0], alternative="two-sided").pvalue
                effect_name = "auroc"
            low, high = bootstrap_effect(xx, yy, continuous)
            source = provenance[name]["source_cohort"]
            rows.append({"cohort": cohort, "signature": name, "n": len(xx),
                         "events": np.nan if continuous else int(yy.sum()), "effect": effect,
                         "effect_name": effect_name, "ci_low": low, "ci_high": high, "p": p,
                         "source_cohort": source, "source_excluded": source != cohort,
                         "standardized_direction": effect if continuous else 2 * (effect - 0.5)})
    out = pd.DataFrame(rows)
    out["q_within_cohort"] = out.groupby("cohort")["p"].transform(bh_adjust)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    signatures, provenance = load_signatures()
    union = sorted(set().union(*map(set, signatures.values())))
    cohort_frames = [read_main(union), read_korea(union), read_imcision(union), read_gse296954(union), read_gse281729(union)]
    scored, coverage = [], []
    for frame in cohort_frames:
        one, cov = score_cohort(frame, signatures, union)
        scored.append(one)
        coverage.extend(cov)
    scores = pd.concat(scored, ignore_index=False).reset_index(names="sample_key")
    coverage_frame = pd.DataFrame(coverage)
    result = metrics(scores, signatures, provenance)
    evaluable = coverage_frame.loc[coverage_frame["coverage"].ge(0.8)].merge(
        result[["cohort", "signature", "source_excluded"]], on=["cohort", "signature"], how="inner"
    )
    external_counts = evaluable.loc[evaluable["source_excluded"]].groupby("signature")["cohort"].nunique()
    gate_signatures = external_counts[external_counts.ge(2)].index.tolist()
    summary = {
        "scoring": "mean of within-cohort gene-level z scores on pretreatment pseudobulk",
        "n_signatures": len(signatures),
        "signature_sizes": {k: len(v) for k, v in signatures.items()},
        "cohort_n": scores.groupby("cohort").size().astype(int).to_dict(),
        "hard_gate": {"rule": ">=4 exact signatures with >=80% coverage in >=2 non-source cohorts",
                      "n_passing_signatures": len(gate_signatures), "passing_signatures": gate_signatures,
                      "pass": len(gate_signatures) >= 4},
        "known_literal_source_issue": "Korea Table S6 contains CSTW; it was retained literally and CTSW was not substituted.",
    }
    pd.DataFrame([{"signature": k, "genes": ";".join(v), **provenance[k]} for k, v in signatures.items()]).to_csv(
        OUT / "signature_catalog.tsv", sep="\t", index=False
    )
    coverage_frame.to_csv(OUT / "signature_coverage.tsv", sep="\t", index=False)
    scores.to_csv(OUT / "patient_signature_scores.tsv", sep="\t", index=False)
    result.to_csv(OUT / "external_replication_metrics.tsv", sep="\t", index=False)
    (OUT / "audit.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    print(result.sort_values(["cohort", "p"])[["cohort", "signature", "n", "effect_name", "effect", "ci_low", "ci_high", "p", "q_within_cohort", "source_excluded"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
