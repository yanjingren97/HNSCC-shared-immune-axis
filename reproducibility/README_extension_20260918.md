# Biological specificity extension 20260918

This version archives code, frozen specifications, gene catalogs, matching diagnostics, software snapshots, source-file hashes and aggregate results for the HNSCC immune gene-set interpretation analysis. It contains no manuscript, expression matrices, patient-level scores or patient-level clinical records.

## Contents and reproduction

Run scripts from the repository root. `scripts/specificity_common.py` resolves that root and the existing results/data hierarchy. The public base repository supplies the original composition and incremental-analysis pipeline. Download source files using the acquisition scripts and manifests, and regenerate patient-level intermediate matrices locally. These inputs are deliberately not redistributed here. `reproducibility/input_file_versions.json` identifies the exact processed/input bytes used; historical source release labels may be less specific than these hashes.

The following are analysis entry points, not a claim that the repository is a one-command pipeline:

1. Original inputs: `acquire_independent_validation_data.py`, `analyze_independent_validation.py`, `build_signature_composition_pseudobulk.py`, `analyze_signature_composition.py` (follow the base repository README for source GEO files and annotation models).
2. Gene provenance: `build_specificity_provenance.py` and the fixed signature catalog.
3. Shared controls: `run_biological_nulls.py`, followed by the BH calculation in `finalize_specificity_outputs.py`; `run_equal_cell_specificity.py` and `complete_technical_audit.py` for sampling/selection.
4. Conditional protein associations: `run_specific_components.py`, then `check_specificity_sensitivity.py`.
5. TCGA abundance: `download_specificity_replication.py`, then `run_specificity_replication.py`.
6. Follow-up diagnostics: `diagnose_ifng6_20260918.py`; revised eight-set framework: `eight_set_framework_20260918.py`.
7. Independent feasibility: `download_independent_candidates.py`, `download_gencode33.py`, `independent_ifng6_20260918.py`. The eligibility gate failed; no independent protein association was calculated.

Core environment: Python 3.12.14, NumPy 2.5.3, pandas 3.0.5, SciPy 1.18.1, GSEApy 1.1.9. Annotation environment: Scanpy 1.12.4, CellTypist 1.7.1. Full recovered metadata is in `software_versions.json`; it is a recovered environment snapshot, not an immutable installation lockfile. Scrublet was invoked through Scanpy. Data-derived scientific plots use Matplotlib, not generative-image synthesis.

The `tmp/python_packages` search path in the scoring helper reflects the original local GSEApy installation; installing the same package in the active environment also permits import. Desktop export destinations in plotting helpers can be changed without affecting analysis definitions. Raw and derived patient-level files must be available at the documented project-relative locations before recomputation. Rebuilding all upstream annotation requires original input downloads and model files; reproducing aggregate-only inspection needs neither patient records nor controlled-access files.

## Statistical families

`null_summary.tsv` contains all 30 random-control tests with `q_all_30_tests`. Main-text immune-background/overlap-preserving examples are CPTAC/eight sets (p=0.000999, q=0.002141) and fixed-T/eight sets (p=0.003996, q=0.006309). Revised eight-catalog deletion tests are a separate exploratory family of three; `q_three_tests` is supplied in the revised deletion summary. Eight conditional score tests, five single-protein tests and five protein-deletion tests are separate BH families. The one-row panel diagnostic q is not the eight-score q.

## Analysis history

The initial nine-set analysis preceded this extension. On 17 September, the source audit identified Shared84 as a derivative of Dysf238, and eight non-derived sets became the main interpretation framework after the initial extension outputs had been inspected. Joint adjustment and subsequent stability diagnostics were exploratory follow-ups. The independent proteomic feasibility plan was fixed before looking at association outcomes. Seed values are integer identifiers, not execution dates: 20260917 (nulls), 20260918 (component analysis and separate follow-up diagnostic streams), 20260919 (TCGA), 20260920 (cell sampling), 20260921 (joint model), 20260922 (independent plan), 20260923 (eight-set recalculation). The primary IFNG6 interval and follow-up diagnostic interval differ slightly because their bootstrap streams differ; both are archived. Clinical analyses retain the original nine-set comparator and were not optimized again.

## Source versions

GO Biological Process 2025 was retrieved from Enrichr; exact source hash, keyword rule, 2,217-gene union and selected term names are archived. Cohort-specific measurable candidate backgrounds exclude the full original gene-set union. Matching uses percentile mean-expression rank, sample-SD rank and fraction of positive deposited values. The CPTAC last measure is a detection proxy, not an assay threshold. See the supplementary methods and `matching_diagnostics.tsv` for matching errors.

This archive is also supplied with the manuscript as `Supplementary_Analysis_Code_20260918.zip`. Package fingerprints are in `reproducibility/PACKAGE_SHA256.tsv`. No new DOI is claimed for this commit-based version.
