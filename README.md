# Published blood RNA signatures in sepsis: version 2.0.0

This release accompanies *Longitudinal behavior and clinical associations of published blood RNA signatures in sepsis*. Version 2.0.0 adds the Frontiers V2 statistical analysis plan, portable analysis code, derived result workbooks, figure-source tables and publication figures. The supplied `_sources/repo_subset_v1.3.1` directory is retained byte-for-byte as immutable historical input.

The eight published formulas are not refitted, recalibrated or replaced. Version 2.0.0 extends the frozen longitudinal baseline with rank-persistence, diagnostic-change concordance, cohort-specific clinical-context, GSE106878, formula-preserving matched-background and contribution-recurrence analyses. Clinical endpoints remain cohort-specific, and GSE106878 is counted once as a separate replication cohort.

The public v1.0.1 archive is retained as superseded history. Version 1.1.2 corrected the SIG001 coefficient, version 1.2.0 reconciled Stage 4 architecture summaries with the signature-specific primary-independent cohort sets, version 1.2.1 completed the manuscript-facing data and documentation package, version 1.2.2 updated inference presentation in Table 2, and version 1.3.0 aligned the reader-facing package with the Biomedical Reports submission.

## Evidence boundary

The study evaluates published blood transcriptomic formulas under repeated sampling. It does not create a classifier, validate clinical monitoring or treatment response, transport original assay thresholds, deconvolve cell abundance, or establish biological mechanism.

## Public inputs

Original expression data are available from NCBI GEO under GSE236713, GSE57065, GSE95233, GSE54514, GSE110487, GSE8121 and GSE106878. The first six cohorts constituted the primary analysis; GSE106878 was used as a separate replication cohort. PRJEB111201 was assessed for external-cohort feasibility using ENA public records but contributes no numerical manuscript result. Source-study terms govern reuse. No direct identifiers or newly linked clinical data are included. Patient and sample identifiers in `data/derived_patient_level` are the pseudonymous identifiers used in the public deposits.

The public v2.0.0 supplement excludes downloaded GEO quick/family records and PMC/PubMed XML snapshots. Project-generated mappings retain official source URLs, retrieval dates and recorded checksums. See `LICENSE_NOTICE.md` and `frontiers_v2/Supplementary_Code_S1/README.md`.

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

## Reproduce the version 2.0.0 additions

From `frontiers_v2/Supplementary_Code_S1`, install `requirements.txt` and run `python run_all.py`. The workflow uses release-relative paths, verifies downloaded GSE106878 inputs against fixed SHA-256 values and writes generated analyses and figures inside the supplement working tree. The formula-preserving benchmark is the longest step. The supplied Data S1 and Data S2 workbooks are the submission-facing result records.

## Reproduce the retained v1.3.0 additions

- In the frozen Stage 3 environment recorded by `code/portable_analysis/03_00_09_environment_lock/requirements_frozen.txt`, run `python code/analyse_common_cohort_sensitivity.py` to regenerate the S27 common-cohort evidence under `reproducibility_evidence/common_cohort_sensitivity`.
- Run `python code/build_biomedical_reports_figures.py` to regenerate the Biomedical Reports figure references and their source extracts under `reproducibility_evidence/generated_biomedical_reports_v1.3.0`.

The figure builder uses the frozen release data and does not alter the v1.2.1 reference outputs. Figures 1 and 2 are scripted layouts. For Figure 3, the script produces a computational reference, while `reference_outputs/biomedical_reports_v1.3.0/Figure_3.tiff` is the author-approved final layout.

## Version 2.0.0 changes

- Added the dated Frontiers V2 statistical analysis plan, portable Python workflow, project-generated cohort mappings and audit checks.
- Added rank-persistence, diagnostic-change concordance, cohort-specific clinical-context, GSE106878, formula-preserving matched-background and gene-contribution recurrence modules.
- Added submission-facing Data S1 and Data S2, 28 figure-source CSV files, and Figures 1–6 plus Supplementary Figures S1–S6 in TIFF, PDF and SVG formats.
- Removed local workstation paths from the public supplement without modifying the frozen `_sources/repo_subset_v1.3.1` history.
- Excluded downloaded GEO quick/family records and PMC/PubMed XML snapshots from the public supplement; official URLs, retrieval dates and recorded checksums remain available for verification.
- Updated version-specific metadata to DOI `10.5281/zenodo.22028159`; the concept DOI remains `10.5281/zenodo.21415496`.

## Version 1.3.1 changes

- Removed obsolete journal-specific labels from retained code comments, internal reproduction identifiers, audit labels, and two table-workbook provenance cells. Historical version meaning and all audit outcomes are preserved.
- Clarified in the Data S1 and Data S2 README worksheets that v1.3.1 is the current release, the corrected v1.2.1 workflow remains the numerical source, and v1.2.2 is retained as legacy documentation history.
- Updated release, citation, DOI, changelog, and manifest metadata to v1.3.1 and DOI `10.5281/zenodo.21862582`.
- Did not rerun or alter the frozen primary analysis. Cohorts, formulas, analysis populations, numeric values, statistical results, figures, and scientific interpretation are unchanged.

## Version 1.3.0 changes

- Adapted reader-facing terminology and file mapping for the Biomedical Reports submission while keeping the public release journal-neutral.
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

- `frontiers_v2/Supplementary_Code_S1`: version 2.0.0 analysis plan, portable Python workflow, frozen v1.3.1 subset and project-generated cohort-audit records.
- `frontiers_v2/Supplementary_Code_S1.zip`: integrity-matched standalone copy of the version 2.0.0 code supplement.
- `frontiers_v2/Supplementary_Data_S1.xlsx` and `frontiers_v2/Supplementary_Data_S2.xlsx`: submission-facing result and figure-source workbooks.
- `frontiers_v2/Figure_Source_Data`: 28 CSV files supporting Figures 1–6 and Supplementary Figures S1–S6.
- `frontiers_v2/Figures`: TIFF, PDF and SVG exports for Figures 1–6 and Supplementary Figures S1–S6.
- `data/derived_patient_level`: corrected de-identified score, paired-change, gene-contribution and pathway-change matrices.
- `data/figure_source_data`: versioned CSV source data for Figures 1–6.
- `data/figure_source_data/biomedical_reports_v1.3.0`: source extracts for the three Biomedical Reports-facing main figures.
- `tables`: Supplementary Data (Data S1), the figure-source workbook (Data S2), implementation registry and main Tables 1–2.
- `reference_outputs/main_figures_v1.2.1`: frozen main-figure references used by the numerical reproduction verifier.
- `reference_outputs/supplementary_figures_v1.2.1`: frozen Supplementary Figures S1–S5, including the corrected S2 terminology.
- `reference_outputs/biomedical_reports_v1.3.0`: Biomedical Reports-facing TIFF files retained unchanged from the v1.3.0 presentation package.
- `editable_author_sources`: retained PowerPoint sources for the earlier author-reviewed layouts of Figures 1, 4, 5, and 6.
- `code/portable_analysis`: release-relative corrected analysis scripts.
- `reproducibility_evidence`: numerical comparisons, clean-run results, sensitivity evidence and the frozen Stage 2 signature-selection audit.
- `reproducibility_evidence/stage2_signature_selection_audit`: the prespecified eligibility criteria, complete 34-candidate registry and frozen A1/A2/B/X reproducibility grades.

## Version and citation

Release identifier: `sepsis-temporal-signatures-v2.0.0`.

Public repository: `https://github.com/LiXinzhuo0425/sepsis-temporal-signatures`.

Concept DOI for all versions: `10.5281/zenodo.21415496`.

Version-specific DOI for v2.0.0: `10.5281/zenodo.22028159` (https://doi.org/10.5281/zenodo.22028159).

Do not cite the superseded v1.0.1 version DOI as the corrected analysis.

## License

Original release code is available under the MIT License. Derived data, documentation and original figure/source-data materials are available under CC BY 4.0. Third-party source records and dependencies retain their original terms and are not relicensed by this repository. See `LICENSE_NOTICE.md`, `LICENSE_CODE` and `LICENSE_DATA`.
