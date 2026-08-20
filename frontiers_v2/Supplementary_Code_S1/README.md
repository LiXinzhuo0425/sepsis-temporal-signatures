# Supplementary Code S1

This directory contains the dated Frontiers V2 statistical analysis plan, the complete Python analysis scripts, the frozen legacy inputs needed by the new modules, source-verified clinical metadata mappings and the code used to generate the manuscript figures.

## Analysis order

1. Reproduce the frozen longitudinal effects.
2. Estimate patient-rank persistence, diagnostic-signature change concordance and exact contribution recurrence.
3. Estimate the prespecified cohort-specific clinical contrasts.
4. Download and verify the two public GSE106878 inputs.
5. Reconstruct the eight fixed formulas in GSE106878.
6. Calculate its secondary rank and treatment-adjusted mortality results.
7. Run the formula-preserving temporal-background benchmark.
8. Generate the manuscript and supplementary figures.

The complete sequence can be run from this directory with:

```bash
python -m pip install -r requirements.txt
python run_all.py
```

The formula-preserving benchmark is the longest step. Deterministic SHA-256-derived random seeds are used throughout.

## Public inputs

`fetch_public_inputs.py` downloads the GSE106878 series matrix and the complete GPL10295 platform table from NCBI GEO. The files are accepted only when their SHA-256 values match the values fixed on 16 August 2026. Raw expression files from the six primary cohorts are not duplicated here because the frozen, derived patient-level inputs used in the analysis are included under `_sources/repo_subset_v1.3.1/`. The public archive also excludes downloaded GEO metadata snapshots and complete PMC/PubMed XML records; `PROVENANCE.csv` retains their accessions, official source URLs, retrieval date and recorded SHA-256 values. Project-generated mappings and the ENA/INSDC API tables required by the feasibility audit remain included.

## Directory map

- `SAP_FRONTIERS_V2_20260816.md`: dated analysis plan and amendments.
- `MANIFEST_SHA256.txt`: integrity record for every file in this archive.
- `_work/`: executable analysis and figure scripts.
- `_sources/repo_subset_v1.3.1/`: frozen formula registry and derived primary-cohort inputs.
- `_sources/cohort_audit/`: source-verified sample and clinical mappings.
- `_work/analysis/`: generated machine-readable results after execution.
- `final/Figures/` and `final/Figure_Source_Data/`: generated figure files and source tables.

## Third-party source boundary

The project MIT and CC BY 4.0 licenses do not relicense third-party records. ENA API tables in `_sources/cohort_audit/` remain subject to ENA/INSDC public-record terms and are attributed in the repository-level `LICENSE_NOTICE.md`. Downloaded GEO quick/family records and PMC/PubMed full-text or abstract XML snapshots are not redistributed in this public archive. The official URLs and checksums retained in `PROVENANCE.csv` permit source-level verification without repackaging those records.

## Interpretation boundary

The clinical endpoints remain cohort-specific. The scripts do not pool incompatible mortality, organ-function and severity definitions. GSE106878 is counted once as an independent replication cohort. Its secondary mortality and rank analyses do not add another independent validation unit.

## Environment

The manuscript analyses were executed with Python 3.9.6. Exact package versions are listed in `requirements.txt`. Results tables in Data S1 and the figure source tables in Data S2 are the submission-facing records.
