#!/usr/bin/env Rscript
# Article: Identification and preliminary validation of senescence-associated secretory phenotype-related candidate biomarkers associated with epilepsy
# Journal: Neurochemical Research
# Authors: Difei Wang, Hui Gan, Yuxin Wu, Baohui Yang, Han Xiao, Yue Qin,
# Yujie Zhu, Haotian Tang, Kaiyi Kang, Yudong Zhou, Bin Zou, Xuan Zhai, Wu Yan
# Correspondence: Xuan Zhai, Department of Neurosurgery, Liangjiang Campus,
# Children's Hospital Affiliated to Chongqing Medical University, Chongqing,
# China; zhaixuan@163.com. Wu Yan, Department of Neurosurgery, The First
# People's Hospital of Urumqi, Urumqi, Xinjiang, China; yanwu1974@163.com
# Reproduce the six-sample GSE190451 EP-versus-control DESeq2 result.
# Usage: Rscript Supplementary_Code_1_GSE190451_DESeq2.R counts.csv groups.csv output_dir
# Optional legacy comparison: append historical.csv before output_dir.
# The historical table is used only for comparison; it does not select genes or enter the model.
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L)) {
  stop("Expected counts.csv groups.csv output_dir, or counts.csv groups.csv historical.csv output_dir")
}
count_file <- normalizePath(args[[1]], mustWork = TRUE)
group_file <- normalizePath(args[[2]], mustWork = TRUE)
history_file <- if (length(args) == 4L) normalizePath(args[[3]], mustWork = TRUE) else NULL
out_dir <- args[[length(args)]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
groups <- read.csv(group_file, row.names = 1, check.names = FALSE)
historical <- if (!is.null(history_file)) read.csv(history_file, check.names = FALSE) else NULL
stopifnot(identical(colnames(counts), rownames(groups)))
stopifnot(identical(colnames(groups), "Type"))
stopifnot(!anyNA(counts), !anyNA(groups$Type), !anyDuplicated(rownames(counts)))
stopifnot(ncol(counts) == 6L, nrow(counts) == 39374L)
count_matrix <- as.matrix(counts)
stopifnot(all(is.finite(count_matrix)), all(count_matrix >= 0),
          all(count_matrix == round(count_matrix)))
storage.mode(count_matrix) <- "integer"
# Preserve all provided symbol-level integer count rows. DESeq2 handles all-zero rows.
groups$Type <- factor(groups$Type, levels = c("control", "EP"))
stopifnot(!anyNA(groups$Type), identical(as.integer(table(groups$Type)), c(3L, 3L)))

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = groups, design = ~ Type)
dds <- DESeq(dds, fitType = "parametric", quiet = TRUE)
res <- results(dds, contrast = c("Type", "EP", "control"),
               pAdjustMethod = "BH", alpha = 0.1,
               independentFiltering = TRUE, cooksCutoff = TRUE)
result <- as.data.frame(res)
result$GeneSymbol <- rownames(result)
result <- result[, c("GeneSymbol", "baseMean", "log2FoldChange", "lfcSE",
                     "stat", "pvalue", "padj")]
write.csv(result, file.path(out_dir, "GSE190451_DESEQ2_all_reproduced.csv"), row.names = FALSE, na = "")

five <- c("SERPINE1", "IGFBP4", "C3", "SPX", "CCL2")
stopifnot(all(five %in% result$GeneSymbol))
verified <- result[match(five, result$GeneSymbol), ]
fields <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
verified$nominal_p_lt_0_05 <- verified$pvalue < 0.05
verified$BH_FDR_lt_0_05 <- verified$padj < 0.05
verified$expression_direction <- ifelse(verified$log2FoldChange > 0, "Up in EP", "Down in EP")
verified$comparison <- "EP versus control"
verified$input <- "01.exprcount_GSE190451.csv; 01.group_GSE190451.csv"
if (!is.null(historical)) {
  stopifnot(all(c("GeneSymbol", fields) %in% colnames(historical)),
            all(five %in% historical$GeneSymbol))
  old <- historical[match(five, historical$GeneSymbol), ]
  for (field in fields) {
    verified[[paste0("historical_", field)]] <- old[[field]]
    verified[[paste0("abs_diff_", field)]] <- abs(verified[[field]] - old[[field]])
  }
  verified$historical_exact_match <- vapply(seq_len(nrow(verified)), function(i) {
    all(vapply(fields, function(field) {
      isTRUE(all.equal(verified[[field]][i], old[[field]][i], tolerance = 1e-8))
    }, logical(1)))
  }, logical(1))
}
columns <- c("GeneSymbol", "comparison", fields, "expression_direction",
             "nominal_p_lt_0_05", "BH_FDR_lt_0_05")
if (!is.null(historical)) {
  columns <- c(columns, "historical_exact_match", paste0("historical_", fields),
               paste0("abs_diff_", fields))
}
verified <- verified[, c(columns, "input")]
write.csv(verified, file.path(out_dir, "GSE190451_FIVE_GENE_VERIFIED_RESULTS.csv"),
          row.names = FALSE, na = "")

if (!is.null(historical)) {
  comparison <- merge(result, historical[, c("GeneSymbol", fields)],
                      by = "GeneSymbol", suffixes = c("_reproduced", "_historical"),
                      sort = FALSE)
  comparison$reproduced_p_lt_0_05 <- comparison$pvalue_reproduced < 0.05
  comparison$historical_p_lt_0_05 <- comparison$pvalue_historical < 0.05
  comparison$reproduced_FDR_lt_0_05 <- comparison$padj_reproduced < 0.05
  comparison$historical_FDR_lt_0_05 <- comparison$padj_historical < 0.05
  write.csv(comparison, file.path(out_dir, "GSE190451_DESEQ2_HISTORY_COMPARISON.csv"),
            row.names = FALSE, na = "")
}
summary_lines <- c(
  sprintf("R version: %s", R.version.string),
  sprintf("DESeq2 version: %s", as.character(packageVersion("DESeq2"))),
  sprintf("Input dimensions: %s genes x %s samples", nrow(count_matrix), ncol(count_matrix)),
  sprintf("Groups: %s", paste(sprintf("%s=%s", names(table(groups$Type)), table(groups$Type)), collapse = ", ")),
  sprintf("Current result rows with non-missing P and BH P: %s", sum(complete.cases(result))),
  sprintf("Current P NA: %s; current BH P NA: %s", sum(is.na(result$pvalue)), sum(is.na(result$padj))),
  sprintf("Reproduced nominal P<0.05: %s", sum(result$pvalue < 0.05, na.rm = TRUE)),
  sprintf("Reproduced BH FDR<0.05: %s", sum(result$padj < 0.05, na.rm = TRUE))
)
if (!is.null(historical)) {
  summary_lines <- c(summary_lines,
    sprintf("Historical result rows: %s", nrow(historical)),
    sprintf("Matched result rows: %s", nrow(comparison)),
    sprintf("Historical rows absent from input: %s", sum(!historical$GeneSymbol %in% rownames(count_matrix))),
    sprintf("Historical rows with current BH P NA: %s", sum(is.na(comparison$padj_reproduced))),
    sprintf("Historical nominal P<0.05: %s", sum(historical$pvalue < 0.05, na.rm = TRUE)),
    sprintf("Historical BH FDR<0.05: %s", sum(historical$padj < 0.05, na.rm = TRUE)),
    sprintf("Five-gene nominal calls agree: %s", all(verified$nominal_p_lt_0_05 == (old$pvalue < 0.05))),
    sprintf("Five-gene FDR calls agree: %s", all(verified$BH_FDR_lt_0_05 == (old$padj < 0.05))),
    "Exact numerical historical reproduction: NO (historical software/preprocessing parameters were not supplied).")
  for (field in fields) {
    left <- comparison[[paste0(field, "_reproduced")]]
    right <- comparison[[paste0(field, "_historical")]]
    summary_lines <- c(summary_lines,
      sprintf("%s: max absolute difference=%.12g; differing NA status=%s",
              field, max(abs(left - right), na.rm = TRUE), sum(is.na(left) != is.na(right))))
  }
} else {
  summary_lines <- c(summary_lines, "Historical comparison: not requested; core DESeq2 results are independent of historical input.")
}
writeLines(summary_lines, file.path(out_dir, "GSE190451_reproduction_summary.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "GSE190451_R_sessionInfo.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
