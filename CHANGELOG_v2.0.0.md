# Version 2.0.0 (2026-08-20)

Version 2.0.0 adds the Frontiers V2 longitudinal-analysis materials while preserving the frozen v1.3.1 repository subset as immutable history. The eight published formulas and the corrected v1.2.1 primary numerical baseline are not refitted, recalibrated or replaced.

## Changes

- Added the dated Frontiers V2 statistical analysis plan, portable Python workflow, project-generated cohort mappings and audit checks.
- Added rank-persistence, diagnostic-change concordance, cohort-specific clinical-context, GSE106878, formula-preserving matched-background and gene-contribution recurrence modules.
- Added submission-facing Data S1 and Data S2, 28 figure-source CSV files, and Figures 1–6 plus Supplementary Figures S1–S6 in TIFF, PDF and SVG formats.
- Made the cohort-audit validation paths release-relative and removed local workstation paths from provenance records.
- Excluded downloaded GEO quick/family records and PMC/PubMed XML snapshots from the public supplement because the repository licenses do not grant uniform redistribution rights; official URLs, retrieval dates and recorded checksums remain in `PROVENANCE.csv`.
- Updated release, citation, DOI, changelog, license notice and manifest metadata to version 2.0.0 and DOI `10.5281/zenodo.22028159`.

## Unchanged scientific boundary

- The 46-file `_sources/repo_subset_v1.3.1` history is byte-identical to the supplied frozen subset.
- Published signature gene membership and algebra are unchanged.
- No assay threshold is transported, and no classifier is refitted or recalibrated.
- The version 2.0.0 additions do not replace the frozen v1.2.1 primary numerical workflow.

Version 1.3.1 and all earlier releases remain part of the immutable version history.
