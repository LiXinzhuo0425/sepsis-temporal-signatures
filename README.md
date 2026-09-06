# Published blood RNA signatures in sepsis version 2.1.0

This code and derived data release accompanies *Early longitudinal trajectories and rank persistence of published blood RNA signatures in sepsis: a multicohort transcriptomic study*. The current study files are in `analysis_v2/`. They preserve eight published formulas without refitting and include the final treatment-arm sensitivity, patient-count clarification, figures and supplementary appendix.

The six primary cohorts contain 299 patients with observed baseline samples; 264 unique patients contribute a T24 or T48 pair (172 T24 pairs and 146 T48 pairs). The 302 unique case IDs include three patients without an observed baseline. GSE106878 remains one separate directional-replication cohort with 47 paired patients. Counts of patients with two recorded time categories and two fixed analysis windows are different descriptive quantities; see `analysis_v2/Supplementary_Code_S1/PATIENT_COUNT_CLARIFICATION_20260906.md`.

## Current files

- `analysis_v2/Supplementary_Data_S1.xlsx`: final analysis-ready results and supporting analyses. `TreatmentSensitivity` reports formal GSE106878 arm comparisons; `S26_DirectionalReplication` retains frozen descriptive strata.
- `analysis_v2/Supplementary_Data_S2.xlsx`: final exact figure-source tables.
- `analysis_v2/Figure_Source_Data/`: CSV exports named after the 28 source-workbook sheets.
- `analysis_v2/Figures/Fig1.tif` through `Fig6.tif`: byte-identical copies of the final supplied main figures.
- `analysis_v2/S1_Appendix.pdf`: final supplementary methods, Tables S1–S13, and Figures S1–S6.
- `analysis_v2/Supplementary_Code_S1/` and its ZIP: portable analysis code, dated plan, source mappings and frozen derived inputs.

The public workbook edits are limited to venue labels and the treatment-sensitivity tab name. All numerical cells and formulas were compared with the final supplied workbooks. Local submission correspondence and administrative forms are not part of this research archive. Venue wording and paths in public code copies were normalized; the plan date, scientific provisions and source numerical values are retained. Original input-file checksums and verification scope are recorded in `RELEASE_ALIGNMENT_v2.1.0.json`.

## Current reproduction workflow

Verify `MANIFEST_SHA256.txt` before executing code. In a separate working copy, change to `analysis_v2/Supplementary_Code_S1`, install `requirements.txt`, and run `python run_all.py`. The workflow reproduces the frozen baseline and longitudinal modules, verifies public input downloads, calculates the separate treatment-arm sensitivity, and builds the figures. The last figure step generates the current selection Figure 1. Final supplied TIFFs and the appendix are the authoritative presentation files; older PDF/SVG/TIFF alternatives are retained under `analysis_v2/Figures/historical_v2.0.1/`.

The scientific interpretation remains bounded: this study does not train a new classifier, validate clinical monitoring or treatment efficacy, transport original assay thresholds, infer cell proportions, or establish causal mechanisms. Clinical endpoints remain cohort-specific.

## Historical workflows and release history

The sections below describe retained earlier analyses and their verification records. Their figure numbering and historical reproduction status are not a claim that version 2.1.0 was rerun in full.

## Public inputs

Original expression data are available from NCBI GEO under GSE236713, GSE57065, GSE95233, GSE54514, GSE110487, GSE8121 and GSE106878. The first six cohorts constituted the primary analysis; GSE106878 was used as a separate replication cohort. PRJEB111201 was assessed for external-cohort feasibility using ENA public records but contributes no numerical manuscript result. Source-study terms govern reuse. No direct identifiers or newly linked clinical data are included. Patient and sample identifiers in `data/derived_patient_level` are the pseudonymous identifiers used in the public deposits.

The public v2 supplement excludes downloaded GEO quick/family records and PMC/PubMed XML snapshots. Project-generated mappings retain official source URLs, retrieval dates and recorded checksums. See `LICENSE_NOTICE.md` and `analysis_v2/Supplementary_Code_S1/README.md`.

## Reproduce Figures 1–6 and Supplementary Figure S1

Before running any code, verify the pristine extracted release from its root with `shasum -a 256 -c MANIFEST_SHA256.txt`. The manifest records the immutable published state. Run the reproduction workflow in a separate working copy because it intentionally regenerates `reproducibility_evidence/generated_v1.2.1` and refreshes the two `reproduction_result_v1.2.1.json` files.

1. Use Python 3.12.
2. Install `environment/requirements_release.txt`.
3. From the release root, run `python scripts/run_reproduction.py`.

The script uses only release-relative paths and regenerates Figures 1–6 plus Supplementary Figure S1 under `reproducibility_evidence/generated_v1.2.1`; it does not overwrite the retained historical project snapshots. Scientific content is verified through byte-identical comparison of all 12 main-figure source-data CSVs and the frozen Supplementary Figure S1 source CSV. Figure 2 and Supplementary Figure S1 are additionally required to match the decoded reference pixels only in the recorded Python/package, macOS arm64, and Arial environment. In another software, operating-system, or font-rendering environment, a clearly labelled `PASS_PLATFORM_RENDER_VARIATION` is permitted only when all scientific source files remain byte-identical and every relevant export check passes. For Figures 1, 3, 4, 5, and 6, the submission reference retains an author-reviewed presentation layout and the verifier records `PASS_AUTHORED_LAYOUT` after source-data comparison. The verifier also checks that TIFF files are RGB, 600 dpi and LZW-compressed, SVG files retain editable text, PDFs exist, and the retained PowerPoint sources remain editable. A successful run writes `reproducibility_evidence/reproduction_result_v1.2.1.json` with `PASS`.

## Reproduce the corrected analysis from GEO inputs

Portable Stage 3 and Stage 4 scripts are in `code/portable_analysis`. Download the six public inputs into the layout documented by `analysis_config.json`, set `SEPSIS_SIGNATURE_ANALYSIS_ROOT` to that analysis directory, install the recorded Python dependencies and the R/DESeq2 versions in `environment/R_packages_stage4.txt`, and run the Stage 3 and Stage 4 entry points below. The release does not redistribute the large raw GEO files.

The bundled DESeq2 1.50.2 library is a macOS arm64 environment snapshot licensed under LGPL-3.0-or-later; other platforms should install the recorded version and dependencies natively. See `LICENSE_NOTICE.md`.

- `03_00_09_environment_lock/run_stage3.py`
- `03_00_09_environment_lock/verify_stage3.py`
- `04_environment/reproducibility_rerun_stage4.py`
- `04_environment/verify_stage4.py`

The retained v1.2.1 numerical workflow includes the corrected SIG001 implementation and the v1.2.0 regression guard requiring every Stage 4 architecture component to match the exact signature-specific primary-independent cohort IDs, cohort count, and patient count from Stage 3. The complete formula registry, threshold-sensitivity evidence, architecture-set verification, Stage 4 log, and semantic rerun comparison are retained in `tables` and `reproducibility_evidence`.

The final workbook-assembly step is presentation-only and runs when `ARTIFACT_TOOL_MJS` points to a compatible workbook renderer. If that optional dependency is absent, the reproducibility runner skips workbook assembly after all scientific CSV/parquet outputs and verification checks have completed; the release workbooks are already included in `tables`.

## Reproduce the version 2 additions

From `analysis_v2/Supplementary_Code_S1`, install `requirements.txt` and run `python run_all.py`. The workflow uses release-relative paths, verifies downloaded GSE106878 inputs against fixed SHA-256 values and writes generated analyses and figures inside the supplement working tree. The formula-preserving benchmark is the longest step. The supplied Data S1 and Data S2 workbooks are the submission-facing result records.

## Reproduce the retained v1.3.0 additions

- In the frozen Stage 3 environment recorded by `code/portable_analysis/03_00_09_environment_lock/requirements_frozen.txt`, run `python code/analyse_common_cohort_sensitivity.py` to regenerate the S27 common-cohort evidence under `reproducibility_evidence/common_cohort_sensitivity`.
- Run `python code/build_presentation_figures.py` to regenerate the version 1.3.0 presentation figure references and their source extracts under `reproducibility_evidence/generated_presentation_v1.3.0`.

The figure builder uses the frozen release data and does not alter the v1.2.1 reference outputs. Figures 1 and 2 are scripted layouts. For Figure 3, the script produces a computational reference, while `reference_outputs/presentation_v1.3.0/Figure_3.tiff` is the author-approved final layout.

## Version 2.0.1 changes

- Removed venue-specific labels from current public paths, documentation, code comments, workbook provenance fields and release metadata.
- Corrected creator order to Xinzhuo Li, Zihao Yan, Xinyue Tu and Yi Gong.
- Preserved all analysis inputs, formulas, numerical outputs, figures and interpretation boundaries from version 2.0.0.
- Rebuilt the standalone code archive and SHA-256 manifests after the packaging-only normalization.

## Version 2.0.0 changes

- Added the dated version 2 statistical analysis plan, portable Python workflow, project-generated cohort mappings and audit checks.
- Added rank-persistence, diagnostic-change concordance, cohort-specific clinical-context, GSE106878, formula-preserving matched-background and gene-contribution recurrence modules.
- Added submission-facing Data S1 and Data S2, 28 figure-source CSV files, and Figures 1–6 plus Supplementary Figures S1–S6 in TIFF, PDF and SVG formats.
- Removed local workstation paths from the public supplement while preserving the frozen numerical source content in `_sources/repo_subset_v1.3.1`.
- Excluded downloaded GEO quick/family records and PMC/PubMed XML snapshots from the public supplement; official URLs, retrieval dates and recorded checksums remain available for verification.
- Updated version-specific metadata to DOI `10.5281/zenodo.22028159`; the concept DOI remains `10.5281/zenodo.21415496`.

## Version 1.3.1 changes

- Removed obsolete journal-specific labels from retained code comments, internal reproduction identifiers, audit labels, and two table-workbook provenance cells. Historical version meaning and all audit outcomes are preserved.
- Clarified in the Data S1 and Data S2 README worksheets that v1.3.1 is the current release, the corrected v1.2.1 workflow remains the numerical source, and v1.2.2 is retained as legacy documentation history.
- Updated release, citation, DOI, changelog, and manifest metadata to v1.3.1 and DOI `10.5281/zenodo.21862582`.
- Did not rerun or alter the frozen primary analysis. Cohorts, formulas, analysis populations, numeric values, statistical results, figures, and scientific interpretation are unchanged.

## Version 1.3.0 changes

- Adapted reader-facing terminology and file mapping for the version 1.3.0 presentation while keeping the public release venue-independent.
- Added Supplementary Data Table S27 and `code/analyse_common_cohort_sensitivity.py`, which derive the common-cohort sensitivity summary from retained inputs without replacing the primary analysis.
- Identified `tables/Figure_Source_Data.xlsx` as manuscript Data S2 and retained its traceable mapping to the figure-source CSVs.
- Updated the submission-facing Figure 3 typography and retained its final 600-dpi TIFF; plotted values, panel content, and numerical source data are unchanged.
- Updated repository, release, citation, DOI, changelog, and manifest metadata to v1.3.0.
- Preserved the corrected v1.2.1 numerical workflow without rerunning the frozen primary analysis.

## Version 1.2.2 changes

- Added modified Knapp–Hartung 95% confidence intervals and Holm-adjusted P values to the manuscript-facing Table 2 for the eight diagnostic signature–window contrasts.
- Removed the threshold-dependent `T48 pattern` column from Table 2.
- Identified four changes supported by both the standard and modified procedures and two estimates that are sensitive to the inferential method.
- Updated release documentation and citation metadata to v1.2.2.
- Preserved the v1.2.1 numerical workflow, cohorts, point estimates, standard Hartung–Knapp results, prediction intervals, heterogeneity estimates, figures, numerical source values, and reproducibility evidence without change.

## Version 1.2.1 changes

- Made `BASELINE_SD` an exact internal positive control by reusing the primary-analysis outputs and restored all 556 scaling-sensitivity records to the machine-readable workbook.
- Renamed healthy-control scaling as an available-cohort sensitivity and documented that the T48 comparison also changes cohort and patient composition.
- Removed superseded evidence-grade fields and replaced ambiguous label-stability codes with explicit threshold-sensitivity terms.
- Clarified GSE106878 as a post-freeze directional-replication cohort; treatment-arm estimates are descriptive, no between-arm interaction was tested, and all strata share the full-cohort baseline standard deviation.
- Revised Figure 6 to show separate baseline-anchored T24 and T48 landmark changes rather than a common three-time-point trajectory.
- Regenerated all six figures in one Python 3.12/Matplotlib 3.11.0 environment and synchronized the submission workbooks and source tables.
- Adopted the author-reviewed submission layouts for Figures 1, 4, 5, and 6 and retained their editable PowerPoint sources. Figure 4's supplied `RPGRIP1` label was corrected to `ZDHHC19` to match the frozen two-largest-contribution rule; the numerical figure-source files were not changed.
- Added the frozen Stage 2 eligibility criteria, 34-candidate registry, and A1/A2/B/X reproducibility grades that define the result-independent eight-signature selection boundary.
- Defined the S10 class-agreement denominator as all-cohort, strict-never-used, and prespecified non-pilot-only analyses. The last corresponds to the frozen internal code label `BLINDED_VALIDATION_ONLY` and does not imply formal masking; MAD scaling is retained separately for material sign reversal.
- Renamed the integrated pathway field to `T48_lowest_FDR_primary_pathway` and added a deterministic pathway-name tie-break for exact FDR ties without changing any selected row or numerical result.
- Added the 74-gene HPA mapping and complete source-file provenance to Supplementary Data Tables S14b and S16.
- Updated Supplementary Figure S2 and its generator to use the reader-facing labels `All cohorts`, `Primary independent`, `Strict never-used`, and `MAD scaling`.

## Version 1.2.0 changes

- Reran Stage 4 architecture summaries with the exact signature-specific primary-independent cohort set for every component.
- Corrected SIG002, SIG003, and SIG004 architecture denominators to four T48 study families and 118 patients; GSE54514 is excluded consistently from score and contribution summaries.
- Recomputed cohort-median dominance, cancellation, modal leading-gene agreement, material-gene count, direction consistency, baseline labels, and all four threshold-perturbation scenarios.
- Changed the descriptive SIG004 label from single-gene-dominant to consistent multigene drift. SIG003 remains single-gene-dominant and changes only in the +20% threshold scenario.
- Added `verify_architecture_primary_set.py` and a machine-readable PASS record that fail on cohort-ID, cohort-count, patient-count, or label drift.
- Regenerated Figure 4, Figure 6, Table 2, Supplementary Data Tables S11/S15/S17, and the matching figure-source workbook.
- Confirmed that Stage 3 effects, confidence intervals, prediction intervals, heterogeneity estimates, standard and modified Knapp–Hartung results, and multiplicity conclusions are unchanged.

## Earlier corrected-release changes

- Corrected the SIG001 Sepsis MetaScore coefficient; v1.0.1 is scientifically superseded.
- Regenerated corrected patient-level scores, paired changes, gene contributions, cohort estimates, meta-analyses and source data.
- Added the complete signature implementation registry and independent formula checks.
- Added architecture-threshold sensitivity evidence and marked threshold-sensitive categorical labels.
- Replaced Figure 5 panel-b hatching with distinct solid colours and thin white segment boundaries; the FDR marker remains independently visible in panel c.
- Replaced residual project-governance wording in Figure 1 with reader-facing study-design wording; no data or numerical output changed.
- Refreshed the outward-facing public-dataset search date to 20 July 2026 while preserving the dated underlying search records.
- Expanded Table 2 with T24/T48 cohort counts and explicit T48 heterogeneity and gene-contribution metrics.
- Locked the figure environment to Python 3.12 and Matplotlib 3.11.0; the clean run reproduces all six submission figures and 12 source-data files.
- Added modified Knapp–Hartung sensitivity using q*=max(1,qHK). Four of the six primary standard-Hartung–Knapp combinations retained multiplicity-adjusted evidence; the primary analysis itself was not replaced.
- Formally adjudicated four updated-search candidate cohorts. GSE106878 is reported as a standalone 47-patient post-freeze T24 directional-replication analysis; E-MEXP-3850, PRJEB111201 and phs003608 were excluded because a reproducible public expression input compatible with the frozen workflow was unavailable.
- Added exact Shapley/Owen attribution definitions, exact formula-source locations, and the public S25-S26 machine-readable tables.
- Corrected the public `data/source_tables` copies so that every Stage 3 and Stage 4 source table is byte-identical to the corrected reproduction project.
- Removed the redundant legacy Supplementary Figure S6 from archived engineering materials and normalized final release filenames without changing scientific content.
- Standardized reader-facing statistical notation to `ρ`, `α`, `β`, `γ`, `τ²`, `I²`, `ΔZ` and `φ`; machine-readable ASCII variable names and canonical MSigDB identifiers remain unchanged for reproducibility.

## Directory map

- `analysis_v2/Supplementary_Code_S1`: version 2 analysis plan, portable Python workflow, frozen v1.3.1 numerical subset and project-generated cohort-audit records.
- `analysis_v2/Supplementary_Code_S1.zip`: integrity-matched standalone copy of the version 2 code supplement.
- `analysis_v2/Supplementary_Data_S1.xlsx` and `analysis_v2/Supplementary_Data_S2.xlsx`: submission-facing result and figure-source workbooks.
- `analysis_v2/Figure_Source_Data`: 28 CSV files supporting Figures 1–6 and Supplementary Figures S1–S6.
- `analysis_v2/Figures`: TIFF, PDF and SVG exports for Figures 1–6 and Supplementary Figures S1–S6.
- `data/derived_patient_level`: corrected de-identified score, paired-change, gene-contribution and pathway-change matrices.
- `data/figure_source_data`: versioned CSV source data for Figures 1–6.
- `data/figure_source_data/presentation_v1.3.0`: source extracts for the three version 1.3.0 presentation figures.
- `tables`: Supplementary Data (Data S1), the figure-source workbook (Data S2), implementation registry and main Tables 1–2.
- `reference_outputs/main_figures_v1.2.1`: frozen main-figure references used by the numerical reproduction verifier.
- `reference_outputs/supplementary_figures_v1.2.1`: frozen Supplementary Figures S1–S5, including the corrected S2 terminology.
- `reference_outputs/presentation_v1.3.0`: TIFF files retained unchanged from the v1.3.0 presentation package.
- `editable_author_sources`: retained PowerPoint sources for the earlier author-reviewed layouts of Figures 1, 4, 5, and 6.
- `code/portable_analysis`: release-relative corrected analysis scripts.
- `reproducibility_evidence`: numerical comparisons, clean-run results, sensitivity evidence and the frozen Stage 2 signature-selection audit.
- `reproducibility_evidence/stage2_signature_selection_audit`: the prespecified eligibility criteria, complete 34-candidate registry and frozen A1/A2/B/X reproducibility grades.

## Version and citation

Release: `sepsis-temporal-signatures-v2.1.0`.

Repository: https://github.com/LiXinzhuo0425/sepsis-temporal-signatures

Concept DOI: https://doi.org/10.5281/zenodo.21415496

Version 2.1.0 DOI: https://doi.org/10.5281/zenodo.22413219

Code is licensed under MIT. Derived data, documentation and original figures are licensed under CC BY 4.0. Third-party records and dependencies retain their own terms. See `LICENSE_NOTICE.md`, `LICENSE_CODE` and `LICENSE_DATA`.
