# S1 File: Reproducible analysis code

This archive accompanies the manuscript “Early longitudinal trajectories and rank persistence of published blood RNA signatures in sepsis: a multicohort transcriptomic study” by Xinzhuo Li, Zihao Yan, Xinyue Tu and Yi Gong. It contains the Python analysis workflow, frozen inputs used by the longitudinal modules, source-verified clinical mappings, the additional treatment-arm sensitivity and the code used to generate the manuscript figures.

The dated statistical analysis plan is distributed under the neutral filename, `SAP_LONGITUDINAL_V2_20260816.md`, because it was frozen for the preceding manuscript version on 16 August 2026. The plan date and scientific provisions are retained; only venue wording and paths are normalized in this public copy. The final revision did not change the fixed formulas, paired-analysis eligibility, or six-cohort primary synthesis. A final patient-level linkage audit clarified descriptive cohort counts by distinguishing unique case IDs from patients with an observed baseline sample. The added treatment-arm sensitivity is explicitly scoped to the separate GSE106878 directional-replication cohort.

## Analysis order

1. Reproduce the frozen longitudinal effects.
2. Estimate patient-rank persistence, diagnostic-signature change concordance and exact contribution recurrence.
3. Estimate the prespecified cohort-specific clinical contrasts.
4. Download and verify the two public GSE106878 inputs.
5. Reconstruct the eight fixed formulas in GSE106878.
6. Calculate arm-specific bootstrap estimates and HC3 treatment-indicator models for the four diagnostic signatures.
7. Calculate the secondary rank and treatment-adjusted mortality results.
8. Run the formula-preserving temporal-background benchmark.
9. Generate the manuscript and supplementary figures.
10. Export the final selection-transparency Figure 1 as `final/Figures/Fig1.tif`; this file supersedes the legacy Figure 1 generated in step 9.

The complete sequence can be run from this directory with:

```bash
python -m pip install -r requirements.txt
python run_all.py
```

The formula-preserving benchmark is the longest step. Deterministic SHA-256-derived random seeds are used throughout.

## Final patient-count clarification

The six primary datasets contain 302 unique case IDs and 299 patients with an observed baseline sample. GSE236713 contains 126 unique case IDs but 124 observed baseline patients; two case IDs have later samples only. GSE54514 contains 36 unique case IDs but 35 observed baseline patients; one nonsurvivor has an isolated Day-4 sample without a baseline link. Across the primary cohorts, 172 patients contribute a T24 pair, 146 contribute a T48 pair, and 264 unique patients contribute at least one of these two early paired analyses.

The earlier count of 277 refers to patients with records in at least two recorded time categories; four of those patients have baseline plus ICU-discharge records without forming a fixed T0-T4 analysis-window pair. The corresponding count with at least two fixed analysis windows is 273. These descriptive distinctions do not alter baseline standardization, paired ΔZ values, primary meta-analysis, rank-persistence analyses, or the GSE106878 replication.

The directory `_sources/repo_subset_v1.3.1/` retains the frozen historical scientific inputs, with neutralized paths and descriptive wording. Consequently, historical source files inside that directory may retain the earlier descriptive labels (for example, 126/36 as unique case-ID counts and a legacy Figure-1 script). The active figure workflow in `_work/` applies the audited 124/35 observed-baseline counts and displays 299 baseline patients. See `PATIENT_COUNT_CLARIFICATION_20260906.md` for the explicit provenance boundary.

## Public inputs

`fetch_public_inputs.py` downloads the GSE106878 series matrix and complete GPL10295 platform table from NCBI GEO. The files are accepted only when their SHA-256 values match the values fixed on 16 August 2026. Raw expression files from the six primary cohorts are not duplicated because the frozen derived inputs are included under `_sources/repo_subset_v1.3.1/`. Downloaded GEO metadata snapshots and complete PMC/PubMed XML records are also excluded; `PROVENANCE.csv` retains accessions, official source URLs, retrieval dates and recorded SHA-256 values. Project-generated mappings and the ENA/INSDC API tables retained during cohort screening remain included for provenance.

## Directory map

- `SAP_LONGITUDINAL_V2_20260816.md`: dated analysis plan and amendments, with neutral filename and unchanged scientific provisions.
- `MANIFEST_SHA256.txt`: integrity record for every file in this archive.
- `_work/`: executable analysis and figure scripts.
- `_work/run_gse106878_treatment_sensitivity.py`: arm-specific bootstrap and HC3/Holm sensitivity.
- `_work/build_selection_figure1.py`: portable Python source for the final Figure 1 layout.
- `_sources/repo_subset_v1.3.1/`: frozen formula registry and derived primary-cohort inputs.
- `_sources/cohort_audit/`: source-verified sample and clinical mappings.
- `_sources/revision_addition/GSE106878_patient_paired_changes.csv`: frozen 47-patient paired input for the treatment-arm sensitivity; regenerated by the reconstruction step when public inputs are fetched.
- `_sources/revision_addition/Figure_S6A_prespecified_pathway_correlations.csv` and `Figure_S6B_hpa_broad_lineage_contributions.csv`: exact source tables supplied with S2 Data, retained for a layout-only Figure S6 re-export. The legend has a dedicated row; no plotted values or inferential results were changed. The S6 plotting function prefers these frozen figure-source tables over legacy workbook sheets.
- `_work/analysis/`: generated machine-readable results after execution.
- `final/Figures/` and `final/Figure_Source_Data/`: generated figure files and source tables.

## Third-party source boundary

`LICENSE_CODE`, `LICENSE_DATA` and `LICENSE_NOTICE.md` are unchanged copies from the public project tag [v2.0.1](https://github.com/LiXinzhuo0425/sepsis-temporal-signatures/tree/v2.0.1), retrieved on 4 September 2026. Their historical release attribution is retained. The full-repository environment paths described in that notice are not all included in this smaller supporting-code archive. This public copy is included in version 2.1.0.

The project MIT and CC BY 4.0 licenses do not relicense third-party records. ENA API tables in `_sources/cohort_audit/` remain subject to ENA/INSDC public-record terms and are attributed in the repository-level `LICENSE_NOTICE.md`. Downloaded GEO quick/family records and PMC/PubMed full-text or abstract XML snapshots are not redistributed in this public archive. The official URLs and checksums retained in `PROVENANCE.csv` permit source-level verification without repackaging those records.

## Interpretation boundary

The clinical endpoints remain cohort-specific. The scripts do not pool incompatible mortality, organ-function and severity definitions. GSE106878 has one role as the separate directional-replication cohort. Its secondary mortality and rank analyses do not add another replication unit.

The GSE106878 treatment-arm sensitivity asks whether the four prespecified directions depend on one randomized arm. It is not an analysis of hydrocortisone efficacy. The reconstructed paired table is restored to original GEO sample order before the seeded bootstrap so the archived S26 arm intervals reproduce exactly.

The 5 September 2026 author-review revision clarifies the S26_DirectionalReplication notes in S1 Data: those rows retain the frozen descriptive strata; the formal between-arm sensitivity is reported in TreatmentSensitivity. Numerical results and the primary analyses are unchanged. In Appendix Table S13B, NA denotes that Holm correction is not applicable to intercepts; machine-readable coefficient outputs retain empty values for these non-applicable fields. The Figure 1 source uses Arial text at 8–12 pt and a three-row layout.

## Environment

The manuscript analyses were executed with Python 3.9.6. Exact package versions are listed in `requirements.txt`. Analysis-ready results, cohort registries and supporting analyses are supplied in S1 Data; exact source tables for the main and supplementary figures are supplied in S2 Data.
