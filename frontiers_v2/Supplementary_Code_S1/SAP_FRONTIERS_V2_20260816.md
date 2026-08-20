# Frontiers V2 statistical analysis plan amendment

Version: 1.1  
Frozen: 16 August 2026, Asia/Shanghai  
Target: *Frontiers in Medicine*, Intensive Care Medicine and Anesthesiology, Original Research  
Backend for all new statistical graphics and new quantitative analyses: Python

## 1. Purpose and status

This document freezes the analyses added for the Frontiers V2 manuscript before their model outputs are examined. It is an amendment to the corrected v1.2.1 numerical workflow and the journal-neutral v1.3.1 release. It is not described as a prospective registration. The six-cohort T24/T48 analysis, study-family eligibility rules, eight fixed signatures, formula implementations, scale definitions, primary effect estimates, missingness sensitivity, and exact contribution decompositions remain unchanged.

The new modules ask whether population-average temporal change is accompanied by preservation of patient rank, whether trajectories differ by independently defined clinical course, whether diagnostic signatures move together within patients, and whether leading formula contributors recur across cohorts. A new public cohort will be attempted as a source-defined external pressure test only when its subject-to-sample mapping and gene-level measurements can be reconstructed without outcome-driven decisions.

## 2. Locked legacy analysis

- Signatures: SIG001, SIG002, SIG003, SIG004, SIG022, SIG023, SIG033 and SIG034.
- Formulas and orientation: unchanged from v1.2.1.
- Strict landmark windows: T24, 12 to 36 h; T48, 36 to 60 h.
- Primary estimand: mean paired baseline-to-landmark change divided by the cohort-specific baseline standard deviation.
- Primary synthesis: signature-specific independent-cohort sets, REML random effects, standard Hartung-Knapp confidence intervals, prediction intervals where estimable, and the existing modified Hartung-Knapp sensitivity.
- T24 and T48 remain separate landmark populations. They are not treated as a formal patient-level T24-to-T48 contrast.
- GSE106878 remains outside the primary synthesis as a post-freeze directional-replication cohort.
- No signature is refitted, recalibrated, reselected or assigned a new cutoff.

## 3. Data and provenance gates

Every modeled record must have an unambiguous public sample identifier, patient identifier, source-defined time, signature value and, where relevant, a source-defined clinical label. Ambiguous mappings are recorded as unresolved and excluded. Clinical labels will not be inferred from expression patterns, plotted points or downstream outcomes.

Mortality, source-defined organ-function response and static baseline severity are separate constructs. Mortality cohorts will be synthesized only if the endpoint horizon and definition are compatible after source verification. Otherwise, cohort-specific estimates will be retained. GSE57065 SAPS II information, if accurately linked at patient level, is descriptive baseline-severity context and is not a longitudinal clinical-change endpoint.

## 4. New analysis families

### 4.1 Patient-rank persistence

Population: every primary cohort-by-signature-by-landmark combination with at least 10 complete baseline/follow-up pairs.

Estimands:

1. Spearman correlation between baseline and follow-up fixed-signature scores.
2. Median absolute displacement between a patient's baseline percentile rank and follow-up percentile rank, calculated within the same paired cohort.

Uncertainty: 2,000 patient-level bootstrap resamples with a deterministic seed derived from dataset, signature and landmark. Percentile confidence intervals will be reported. A bootstrap standard error on Fisher's z scale will be used for synthesis.

Synthesis: random-effects Fisher-z synthesis only when at least three eligible independent study families contribute. REML is used for between-study variance. Hartung-Knapp confidence intervals and prediction intervals are reported where estimable. Results are transformed back to Spearman-correlation units. No correlation cutoff defines stability.

Multiplicity: descriptive estimates are reported for all eight signatures. Holm adjustment is applied across the four diagnostic signatures within each landmark for tests of a zero pooled correlation, while interpretation rests on effect size and uncertainty.

### 4.2 Clinical-course anchoring

Primary signatures: SIG001, SIG002, SIG003 and SIG004. The other four signatures are secondary and placed in supplementary material.

Primary clinical endpoints:

- Source-verified mortality in GSE54514 and GSE95233. A common mortality synthesis is allowed only after endpoint compatibility is confirmed.
- Source-defined organ-function response in GSE110487, retained as a separate endpoint family.

Additional context:

- GSE236713 mortality is an external robustness analysis if its endpoint definition is verified.
- GSE57065 source-defined SAPS II category is descriptive baseline-severity context if exact patient mapping is available.

For a two-measurement cohort, the primary model is a linear mixed model of baseline-standardized oriented score with time, clinical group and their interaction, plus a patient random intercept. The time-by-group interaction is the target estimate. Because each included patient has two measurements, a group contrast in paired change is fitted as an algebraically aligned sensitivity analysis with heteroskedasticity-robust and patient-bootstrap uncertainty.

For a cohort with more than two verified measurements, time is categorical and the model includes categorical time-by-clinical-group terms and a patient random intercept. No linear trajectory is imposed. Source-defined days are retained and are not relabeled as strict T24 or T48 unless actual hours meet the locked landmark window.

Multiplicity: Holm correction across the four diagnostic signatures within each compatible endpoint family. Mortality and organ-function response remain separate families. All fitted interaction estimates, confidence intervals, raw P values and adjusted P values are reported. Model convergence, patient counts and observation counts accompany every estimate.

No cutoff search, monitoring AUC, decision-curve analysis, nomogram or new prognostic model is performed.

### 4.3 Cross-signature longitudinal concordance

Population: patients with paired change estimates for all four diagnostic signatures within a cohort and landmark.

Estimands:

- Six pairwise Spearman correlations among within-patient standardized changes.
- Pairwise sign agreement, with exact zeros reported separately.
- Proportion of patients for whom all four changes share the same non-zero direction.

Cohort estimates are synthesized on the Fisher-z scale only when at least three eligible independent study families contribute. All six pairs are retained. Holm correction is applied across six pooled pairwise-correlation tests within each landmark.

### 4.4 External longitudinal pressure test

PRJEB111201 is designated as a source-defined external 2026 pressure-test cohort. It does not enter the locked primary meta-analysis. Day 1 to Day 3 is the primary external contrast and is not called strict T48. Day 1 to Day 7 is secondary.

The module proceeds only if subject-to-run mapping is unambiguous, at least eight patients have valid paired Day 1 and Day 3 measurements, and a defensible gene-level expression matrix can be obtained or reconstructed with complete formula-gene coverage. Processing and exclusions must be independent of signature direction or clinical outcome. Every legally computable fixed signature is reported. If these gates fail, the attempt and reason are documented without substituting a different cohort after viewing results.

### 4.5 Contribution architecture extension

The exact v1.2.1 patient-level decomposition is retained. For each patient-score pair, the sum of component contributions must reproduce total score change with absolute error below 1e-8.

For each signature and landmark, the leading absolute contributor is identified within each eligible independent cohort. Leading-contributor recurrence is the proportion of cohorts sharing the modal leading gene; ties are reported. Dominance, cancellation and recurrence are displayed alongside drift and rank persistence descriptively. No correlation significance test is performed across only eight signatures.

### 4.6 Formula-preserving background benchmark

This module is optional. It enters the manuscript only if gene replacement preserves each signature's algebra and coefficients, candidate genes are matched using baseline-only mean and variance bins, at least 95% of 1,000 iterations succeed, and no follow-up or outcome information enters matching. A formula that cannot be randomized fairly is marked unavailable. Failure of this gate removes the module rather than changing the null design after viewing results.

## 5. Missing data and exclusions

New repeated-measurement analyses use complete pairs for the stated contrast. No outcome imputation is performed. Existing landmark missingness and location-shift sensitivity remain part of the evidence package. New exclusions are permitted only for unresolved identity/time mapping, failed prespecified measurement QC, missing formula genes, or a recorded model/data gate. Reasons are retained in an exclusion log.

## 6. Reporting and claim boundary

Positive, null and discordant results are retained. The manuscript may describe temporal behavior, rank persistence, clinical-course association, concordance and formula contribution. Associations with outcome or source-defined response do not establish treatment monitoring, prognostic utility, causality or transport of original assay thresholds. Mathematical contribution does not identify causal genes or cell types.

## 7. Completion gates

The V2 analysis is complete when:

- the legacy counts and manuscript-facing numerical outputs are reproduced from retained sources;
- every clinical record has traceable sample-to-patient-to-time-to-label mapping;
- rank persistence is generated for all eligible combinations with 2,000 bootstraps;
- at least two independently sourced clinical-anchor cohorts are analyzed;
- the four diagnostic signatures are fully reported in each primary endpoint family;
- GSE106878 remains an independent directional replication;
- PRJEB111201 is attempted under the locked gate and its success or failure is documented;
- contribution reconstruction passes the 1e-8 identity threshold;
- source tables, deterministic seeds, software versions and exclusions are archived;
- the article reports all material null and discordant findings within the stated claim boundary.

## 8. Independent statistical review amendment (version 1.1)

This quality-control amendment was added on 16 August 2026 after an initial implementation dry run and before the final analysis outputs were frozen. The changes respond to independently identified design issues and were not selected according to significance, direction or graphical appearance. All superseded dry-run files are replaced in the final package.

- Signature-specific primary-independent cohort sets are applied to every pooled rank, concordance and contribution result. For a cross-signature pair, a cohort is eligible only when it is eligible for both formulas.
- Patient percentile rank is defined as `(average rank - 1)/(n - 1) × 100`. Ranks and median displacement are recomputed within each bootstrap sample.
- Standard and modified Hartung-Knapp intervals are reported for new small-k syntheses. Modified inference uses `q*=max(q,1)`. Wording follows the more conservative modified result when procedures differ.
- For two-visit clinical cohorts, the adverse-minus-favorable paired-change contrast with HC3 covariance and a status-stratified patient bootstrap is the inferential estimate. The random-intercept time-by-status model is retained as an algebraically aligned consistency analysis. Baseline-adjusted follow-up ANCOVA is added as a sensitivity analysis.
- The combined GSE95233 D2/D3 contrast adjusts for follow-up day. Day-specific estimates are reported as a sensitivity.
- GSE54514 T48 is the main source-defined survival contrast, T24 is secondary, and later days are descriptive with group-specific numbers at risk. Its survival horizon remains unspecified. SIG002-SIG004 results are overlap-flagged exploratory because of possible study-program reuse.
- Clinical endpoints remain cohort-specific. No clinical meta-analysis is performed with fewer than three endpoint-compatible, independent cohorts.
- GSE57065 remains a static baseline-severity family rather than a recovery endpoint.
- A cohort's leading contributor is the gene with the largest absolute cohort-mean signed contribution. Numerical ties within an absolute tolerance of `1e-12` are retained as sets. Patient-level dominance and cancellation remain separate quantities.
- The treatment-adjusted GSE106878 28-day mortality analysis, if reported, is secondary and does not count as a second independent replication role.
