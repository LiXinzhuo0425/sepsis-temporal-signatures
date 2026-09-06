# current patient-count clarification

Date: 6 September 2026

This note records the final patient-level linkage clarification used for the current submission package. It changes descriptive count labels only; the fixed formulas, baseline standardization, paired effects, primary random-effects synthesis, rank persistence, clinical contrasts, contribution decomposition, and GSE106878 analyses were independently reproduced and were not changed by this clarification.

## Primary-cohort count definitions

- Six primary datasets: 302 unique case IDs.
- Patients with an observed baseline sample: 299.
- T24 paired patients: 172.
- T48 paired patients: 146.
- Unique patients contributing a T24 or T48 primary pair: 264.
- Patients with records in at least two recorded time categories: 277.
- Patients with at least two fixed analysis-window records: 273.

### GSE236713

The source contains 126 unique case IDs. Two IDs do not have a baseline sample: patient 80 (Day 5 only in the audited linkage) and patient 84 (ICU-discharge sample only). Therefore the observed-baseline denominator is 124. Both later-only samples were part of the full-array upstream preprocessing/normalization context, but neither entered the baseline mean/SD, paired ΔZ, primary synthesis, or rank-persistence analyses. A duplicated Day-1 record for patient 105 was handled deterministically as defined by the locked workflow and did not duplicate that patient in baseline standardization.

### GSE54514

The source contains 36 unique case IDs. One record, NS_20 / GSM1317938, is an isolated Day-4 sample without a baseline link. Therefore the observed-baseline denominator is 35. Baseline mean/SD used the 35 observed baseline samples; primary T24 and T48 paired counts remain 31 and 28, respectively. The isolated Day-4 record did not enter the primary paired analyses or rank-persistence analyses.

## GSE95233 missingness

For the final submission package, the author-confirmed S22 / Table S9A / Table S9B treatment is retained unchanged. This count clarification does not modify that sensitivity analysis.

## Frozen-source boundary

`_sources/repo_subset_v1.3.1/` retains frozen historical scientific inputs, with venue wording and paths normalized in this public copy. Some historical fields or scripts inside that snapshot therefore retain earlier descriptive labels. The active current scripts under `_work/` apply the final audited baseline-count labels for generated Figure 1 and its figure-facing cohort metadata.
