"""Unified ssGSEA sensitivity analysis across five HNSCC cohorts."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/published_signature_external_audit"
PACKAGES = ROOT / "tmp/python_packages"
if PACKAGES.exists():
    sys.path.insert(0, str(PACKAGES))
import gseapy as gp  # noqa: E402


def load_catalog() -> tuple[dict[str, list[str]], dict[str, str]]:
    catalog = pd.read_csv(OUT / "signature_catalog.tsv", sep="\t")
    sets = {row.signature: str(row.genes).split(";") for row in catalog.itertuples()}
    source = {row.signature: str(row.source_cohort) for row in catalog.itertuples()}
    return sets, source


def read_bulk() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(ROOT / "data/external/GSE281729/GSE281729_expression_logTPM.txt.gz", sep="\t", header=None, dtype=str)
    samples = raw.iloc[13, 2:].astype(str).tolist()
    expr = raw.iloc[14:, 2:].apply(pd.to_numeric, errors="coerce")
    expr.index = raw.iloc[14:, 1].astype(str).str.upper()
    expr.columns = samples
    expr = expr.groupby(level=0).mean()
    endpoint = pd.read_csv(ROOT / "results/transportability_screen/gse281729_scores.tsv", sep="\t").set_index("sample_key")
    common = expr.columns.intersection(endpoint.index)
    meta = pd.DataFrame(index=common)
    meta["cohort"] = "GSE281729_bulk"
    meta["regimen"] = "neoadjuvant_immunotherapy"
    meta["continuous_outcome"] = np.nan
    meta["binary_outcome"] = endpoint.loc[common, "OVERALL_PATH_GE20"].astype(int)
    return expr.loc[:, common], meta


def load_cohorts() -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    clinical = pd.read_csv(ROOT / "results/clinical_incremental_pilot/cohort_and_features.tsv", sep="\t").set_index("patient_id")
    main_expr = pd.read_csv(OUT / "main_full_pretreatment_T_pseudobulk_logcpm.tsv.gz", sep="\t", index_col=0).fillna(0.0)
    main_meta = clinical.loc[main_expr.columns, ["regimen", "pathologic_response_percent"]].rename(
        columns={"pathologic_response_percent": "continuous_outcome"})
    main_meta["binary_outcome"] = main_meta["continuous_outcome"].ge(50).astype(int)
    main_meta["cohort"] = "GSE288199"
    cohorts = {"GSE288199": (main_expr, main_meta)}
    cache = OUT / "full_pseudobulk"
    for cohort in ["GSE286827", "GSE232240_IMCISION", "GSE296954"]:
        expr = pd.read_csv(cache / f"{cohort}_logcpm.tsv.gz", sep="\t", index_col=0).fillna(0.0)
        meta = pd.read_csv(cache / f"{cohort}_metadata.tsv", sep="\t", index_col=0)
        meta.index = meta.index.astype(str)
        expr.columns = expr.columns.astype(str)
        cohorts[cohort] = (expr.loc[:, meta.index], meta)
    cohorts["GSE281729_bulk"] = read_bulk()
    return cohorts


def main() -> None:
    signatures, sources = load_catalog()
    cohorts = load_cohorts()
    score_frames, coverage_rows = [], []
    for cohort, (expr, meta) in cohorts.items():
        input_frame = expr.reset_index(names="Gene")
        run = gp.ssgsea(data=input_frame, gene_sets=signatures, outdir=None, no_plot=True,
                        sample_norm_method="rank", min_size=2, max_size=500,
                        permutation_num=0, threads=1, seed=20260910, verbose=False)
        score = run.res2d.pivot(index="Name", columns="Term", values="NES").astype(float)
        score = score.join(meta, how="inner")
        score["sample_key"] = score.index
        score_frames.append(score.reset_index(drop=True))
        measured = set(expr.index[expr.max(axis=1).gt(0)])
        for name, genes in signatures.items():
            coverage_rows.append({"cohort": cohort, "signature": name, "n_genes": len(genes),
                                  "n_present": len(set(genes) & measured),
                                  "coverage": len(set(genes) & measured) / len(genes)})
    scores = pd.concat(score_frames, ignore_index=True)
    rows = []
    for cohort, frame in scores.groupby("cohort"):
        continuous = frame["continuous_outcome"].notna().all()
        y = frame["continuous_outcome" if continuous else "binary_outcome"].to_numpy(float)
        for signature in signatures:
            x = frame[signature].to_numpy(float)
            if continuous:
                test = spearmanr(x, y)
                value, p, kind = test.statistic, test.pvalue, "spearman_rho"
            else:
                test = mannwhitneyu(x[y == 1], x[y == 0], alternative="two-sided")
                value = test.statistic / (np.sum(y == 1) * np.sum(y == 0))
                p = test.pvalue
                kind = "auroc"
            rows.append({"cohort": cohort, "signature": signature, "n": len(frame),
                         "effect_name": kind, "effect": value, "p": p,
                         "standardized_direction": value if continuous else 2 * (value - 0.5),
                         "source_cohort": sources[signature], "source_excluded": sources[signature] != cohort})
    metrics = pd.DataFrame(rows)
    metrics["q_within_cohort"] = metrics.groupby("cohort")["p"].transform(
        lambda x: pd.Series(np.minimum.accumulate(
            (np.sort(x.to_numpy()) * len(x) / np.arange(1, len(x) + 1))[::-1]
        )[::-1], index=x.sort_values().index).reindex(x.index).clip(upper=1).to_numpy())
    scores.to_csv(OUT / "multicohort_ssgsea_scores.tsv", sep="\t", index=False)
    metrics.to_csv(OUT / "multicohort_ssgsea_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(coverage_rows).to_csv(OUT / "multicohort_ssgsea_coverage.tsv", sep="\t", index=False)
    audit = {"method": "gseapy ssGSEA, rank normalization", "gseapy_version": gp.__version__,
             "cohort_n": scores.groupby("cohort").size().astype(int).to_dict(), "signatures": len(signatures)}
    (OUT / "multicohort_ssgsea_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    print(metrics.sort_values(["cohort", "p"])[["cohort", "signature", "effect_name", "effect", "p", "q_within_cohort", "source_excluded"]].to_string(index=False))


if __name__ == "__main__":
    main()
