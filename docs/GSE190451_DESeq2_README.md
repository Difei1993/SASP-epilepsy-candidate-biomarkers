# GSE190451 raw-count DESeq2

`scripts/GSE190451_DESEQ2_REPRODUCIBLE.R` is byte-identical to the R script in the manuscript's Supplementary Code 1. It analyzes 39,374 symbol-level raw-count rows from six GSE190451 samples, with three controls and three epilepsy samples. The model is `~ Type`; the contrast is epilepsy versus control; P values are adjusted by the Benjamini-Hochberg method. The optional historical table is used only for comparison and is not part of model fitting.

After obtaining and extracting Supplementary Code 1, run from the repository root:

```text
Rscript scripts/GSE190451_DESEQ2_REPRODUCIBLE.R counts.csv groups.csv output_dir
Rscript scripts/GSE190451_DESEQ2_REPRODUCIBLE.R counts.csv groups.csv historical.csv output_dir_with_history
```

The second command is optional. It produces a historical comparison; the first produces the same current five-gene DESeq2 values. Inputs must have matching sample IDs in the same order, a `Type` group column with `control` and `EP`, and non-negative integer counts. The historical CSV, if supplied, must include `GeneSymbol`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, and `padj`. On Windows with a non-UTF-8 R locale, use an ASCII-only working path.

The script writes all-gene results, a five-gene table, a run summary and sessionInfo; optional historical comparison outputs are generated only by the second command. The verification run used R 4.6.1 and DESeq2 1.52.0; package details are in `docs/GSE190451_R_sessionInfo.txt`.

The five-gene directions and nominal/adjusted significance calls agree with the historical table, but exact historical numerical reproduction was not achieved. SERPINE1, IGFBP4 and C3 have BH-adjusted P < 0.05; SPX has nominal P < 0.05 but BH-adjusted P about 0.06175; CCL2 does not have nominal P < 0.05. This check does not change the final four-gene panel or constitute new feature selection.
