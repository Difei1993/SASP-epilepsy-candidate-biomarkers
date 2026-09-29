# SASP-related candidate biomarkers in epilepsy

Analysis scripts supporting the manuscript **Identification and preliminary validation of senescence-associated secretory phenotype-related candidate biomarkers in epilepsy**.

## Public datasets

- **GSE256068** — discovery cohort.
- **GSE190451** — independent human cortical validation dataset (3 TLE and 3 non-epileptic controls).

Expression data are not redistributed here; they can be downloaded from GEO.

## Revised reproducibility scripts

Run from the repository root:

1. `scripts/00_config.R`
2. `scripts/01_prepare_GSE190451.R`
3. `scripts/02_candidate_validation_GSE190451.R`
4. `scripts/03_clinical_association_GSE256068.R`
5. `scripts/99_session_info.R`

The original feature-prioritization workflow produced five genes for cross-dataset assessment: **SERPINE1, CCL2, IGFBP4, C3, and SPX**. The GSE190451 script evaluates these genes without hard-coding the final supported set.

### Statistical note for GSE190451

The GEO record provides processed **TPM** expression values for six samples. The revised validation code analyzes `log2(TPM + 1)` values using **limma** moderated linear modeling and exports both nominal and Benjamini-Hochberg-adjusted p values. This avoids relying on a Wilcoxon significance rule in a 3-versus-3 dataset.

## Legacy analysis

`legacy_original_scripts/README.md` documents the original code archive. The original archive contained machine-specific absolute paths and an earlier GSE134697 validation workflow; those are not used for the revised GSE190451 validation scripts above.

## Manual database/web steps

STRING, miRNet, GeneMANIA, Cytoscape, PubChem and CB-Dock2 steps are described in `docs/manual_database_steps.md`.

## Contact

Difei Wang  
Department of Neurosurgery, Children's Hospital of Chongqing Medical University
