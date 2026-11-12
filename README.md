# SASP-related candidate biomarkers in epilepsy

Analysis materials for *Identification and preliminary validation of senescence-associated secretory phenotype-related candidate biomarkers associated with epilepsy*.

The final manuscript discusses SERPINE1, IGFBP4, C3 and SPX as four candidate biomarkers. CCL2 was part of an earlier five-gene cross-dataset assessment and is not part of the final four-gene panel.

## Data and code scope

The study uses publicly available GEO datasets, including GSE256068 and GSE190451. Source expression data are not redistributed by this repository. The raw-count matrix, sample-group table and reproducibility script used for the GSE190451 Figure 2e check are included in the manuscript's Supplementary Code 1 submission. During peer review, contact the corresponding authors for that supplementary archive if it is not yet publicly posted by the journal.

The repository contains the verified GSE190451 raw-count DESeq2 reproduction, the archived GSE256068 cohort/Figure 3 script, and distinct historical or supplementary scripts. It is not a complete end-to-end reproduction of every manuscript figure or the corrected Supplementary Figure 1.

## GSE190451 analysis used for Figure 2e

The manuscript's Figure 2e annotations use an epilepsy-versus-control comparison of the GSE190451 **raw count** matrix with DESeq2, a `~ Type` design, control as reference, three samples per group, and Benjamini-Hochberg adjustment. Nominal and adjusted P values have distinct meanings. `scripts/GSE190451_DESEQ2_REPRODUCIBLE.R` implements this comparison, with run instructions in `docs/GSE190451_DESeq2_README.md` and a five-gene check in `docs/GSE190451_FIVE_GENE_VERIFIED_RESULTS.csv`.

An independent run using the archived count and group tables agrees with the reported five-gene expression directions and nominal/adjusted significance calls. It does **not** exactly reproduce every historical coefficient and P value; the original complete software/preprocessing provenance was unavailable. The results table preserves the differences instead of replacing the historical table.

## GSE256068 cohort and Figure 3

The archived script for the verified 133-sample cohort and current Figure 3 is supplied as `scripts/GSE256068_FIGURE3_VERIFIED.R`. Run instructions, limits and the archived software inventory are in `docs/GSE256068_Figure3_README.md` and `docs/GSE256068_Figure3_sessionInfo.txt`. These files were published from the verified archive; no analysis or figure was rerun during the publication update.

## Historical TPM/limma analysis

The existing `scripts/01_prepare_GSE190451.R` and `scripts/02_candidate_validation_GSE190451.R` analyze processed `log2(TPM + 1)` expression values with limma-trend (`eBayes(trend = TRUE)`). They remain available as a distinct historical or supplementary analysis. They are not the statistical source of Figure 2e's raw-count DESeq2 annotations. Keep the `legacy_original_scripts/` directory as historical provenance.

The existing GSE256068 clinical-association script should likewise not be treated as an exact reproduction of the corrected four-gene Supplementary Figure 1 without a separate comparison against the final manuscript's 121-case, two-sided Spearman and Wilcoxon specifications.

## Contact

Difei Wang, Department of Neurosurgery, Children's Hospital of Chongqing Medical University.

## Publication update on 12 November 2026

This update publishes previously verified code and outputs and preserves all existing scripts. The final four-gene panel, frozen figures and reported statistics were not changed. Historical and independent DESeq2 numerical results are labelled separately in the five-gene CSV. See `docs/REPRODUCIBILITY_UPDATE_20261112.md`.
