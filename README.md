# HNSCC shared immune-expression signal

Analysis code and aggregate results supporting *Distinct immune gene sets capture a shared transcriptional signal in HNSCC*.

## Scope

This research snapshot includes nine fixed gene sets, analysis configurations, patient-level analysis code, and aggregate outputs for T-cell composition controls, gene-overlap sensitivity, incremental response prediction, CPTAC biological references, and GSE284162 response discrimination. The manuscript, raw expression matrices, patient-level scores/outcomes, patient decoding tables, and private correspondence are not distributed.

## Main outputs

| Analysis | File |
|---|---|
| 22- and 18-patient composition sensitivity | `results/signature_revision_20260915/composition/redundancy_metrics.tsv` and `paired_structure_changes.tsv` |
| Absolute baseline and augmented LOPO MAE, differences and CIs | `results/signature_revision_20260915/incremental/incremental_MAE.tsv` |
| Three planned CPTAC associations, raw and BH-adjusted p values | `results/independent_validation_20260915/CPTAC_associations.tsv` |
| External unrestricted and platform-restricted permutations | `results/independent_validation_20260915/GSE284162_response_associations.tsv` |
| CPTAC gene-overlap sensitivity | `results/independent_validation_20260915/gene_disjoint/` |

MAE units are percentage points; improvement = baseline MAE minus augmented-model MAE. All 10 improvement confidence intervals include zero. The CD3-only partial correlation (0.750) is distinct from the exploratory dual-adjusted correlation (0.629). Shared score structure is not evidence of clinical predictive utility or a cell-intrinsic mechanism. Uniform ssGSEA does not reconstruct the original assay-specific predictors. Small samples and incompletely established cross-study patient independence limit inference. Configurations are analysis records, not prospective registration.

## Reproduction and data access

The aggregate tables can be inspected without downloading patient data. Full recomputation requires source expression matrices, patient annotations and outcomes at the paths specified in the scripts. These inputs are deliberately not bundled; this is not a self-contained patient-data release. Source studies and processed data:

- [GSE288199](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288199): 27 prediction patients; nested composition subsets of 22 and 18.
- [GSE286827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286827): 26 prediction patients.
- [GSE296954](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296954): 10 prediction patients; source cohort excluded for GZMK40 incremental evaluation.
- [CPTAC-HNSCC](https://www.linkedomics.org/data_download/CPTAC-HNSCC/): 108 matched patients, 83 with CD3 IHC; biological reference, not immunotherapy-response validation.
- [GSE284162](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284162): 12 pretreatment chemotherapy–immunotherapy samples, four author-labelled MPR.

Use a scientific Python environment with numpy, pandas, scipy, scikit-learn, threadpoolctl and GSEApy 1.1.9. Raw single-cell preprocessing additionally uses scanpy, anndata, celltypist, scrublet, h5py and rdata. These are dependency names, not a verified lockfile.

After obtaining and preparing the inputs, run from the repository root:

```sh
python scripts/build_signature_composition_pseudobulk.py
python scripts/analyze_signature_composition.py
python scripts/check_signature_composition_robustness.py
python scripts/evaluate_signature_incremental_value.py
python scripts/acquire_independent_validation_data.py
python scripts/analyze_independent_validation.py
python scripts/check_CPTAC_gene_disjoint_reference.py
```

The composition builder depends on previously generated cell annotations; see its input paths and the annotation scripts. Incremental prediction depends on the original analytical score/outcome table. Some historical scripts retain local input conventions and require configuration before use. Acquisition scripts download public files and may need working source servers. No analyses were refitted for this publication snapshot.

`SHA256SUMS.tsv` records the exact public files. This repository does not assign a license to third-party data or gene signatures.
