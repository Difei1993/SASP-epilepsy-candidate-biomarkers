# Verified GSE256068 sample reconstruction and Figure 3

`scripts/GSE256068_FIGURE3_VERIFIED.R` is the archived script used to construct the current Figure 3, including the 133-sample reconstruction. Its archived TIFF is pixel-identical to the frozen submission Figure 3. The script has been published without rerunning the model, calibration or figure export during this final editorial update.

## Inputs and command

Obtain the NCBI-generated GSE256068 TPM matrix and GEO series-matrix metadata from the [GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256068). Source matrices and individual-level clinical records are not redistributed here. Pass the two local source files explicitly:

```text
Rscript scripts/GSE256068_FIGURE3_VERIFIED.R <TPM.tsv.gz> <series_matrix.txt.gz> <output_directory> 1000
```

The TPM file must have `GeneID` in its first column and GSM sample accessions as the remaining columns. This uses the NCBI-generated matrix and the historical title-prefix selection, not the full submitter-processed matrix. The script matches public sample identifiers, excludes the ControlHippocampus prefix and samples absent from the NCBI matrix, orders TSC, FCD2a, FCD2b, TLE-HS and ControlCortex samples, and asserts 121 epilepsy samples and 12 controls. The cohort spans brain regions and pathological subtypes; it is not a homogeneous temporal-cortex cohort.

The model uses log2(TPM + 1) for SERPINE1, IGFBP4, C3 and SPX, epilepsy as 1 and control as 0, inverse-logit probabilities, seed 3 and 1,000 bootstrap repetitions. Calibration and decision-curve analyses use that same discovery cohort. These are exploratory internal assessments and do not establish clinical utility or external calibration. See `docs/Figure3_MODEL_EXPLANATION.txt` and `docs/GSE256068_Figure3_sessionInfo.txt` for the archived interpretation and software inventory.

The older `scripts/03_prepare_GSE256068_metadata.R` is a general metadata helper. It does not implement this complete final 133-sample selection. Likewise, `scripts/04_clinical_association_GSE256068.R` must not be treated as the exact source of the corrected Supplementary Figure 1, whose final reporting uses 121 matched epilepsy cases and the prespecified two-sided Spearman/Wilcoxon comparisons.
