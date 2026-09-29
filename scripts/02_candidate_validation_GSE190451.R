# 02_candidate_validation_GSE190451.R
# Cross-dataset assessment of the five feature-prioritized candidate genes.
#
# GSE190451 contains 3 TLE and 3 non-epileptic control samples and provides TPM.
# We use limma moderated linear modeling on log2(TPM + 1), and report both
# nominal and BH-adjusted p values.

source("scripts/00_config.R")

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Install Bioconductor package 'limma'.")
}

expr_file <- file.path(PROCESSED_DIR, "GSE190451_log2TPM.csv")
group_file <- file.path(PROCESSED_DIR, "GSE190451_groups.csv")

if (!file.exists(expr_file) || !file.exists(group_file)) {
  stop("Run scripts/01_prepare_GSE190451.R first.")
}

expr <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
meta <- read.csv(group_file, stringsAsFactors = FALSE)

stopifnot(identical(colnames(expr), meta$sample))
group <- factor(meta$group, levels = c("control", "EP"))

missing <- setdiff(INITIAL_CANDIDATES, rownames(expr))
if (length(missing) > 0L) {
  stop("Candidate genes missing from GSE190451: ",
       paste(missing, collapse = ", "))
}

candidate_expr <- expr[INITIAL_CANDIDATES, , drop = FALSE]

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contrast <- limma::makeContrasts(EP - control, levels = design)

fit <- limma::lmFit(candidate_expr, design)
fit <- limma::contrasts.fit(fit, contrast)
fit <- limma::eBayes(fit, trend = TRUE)

tab <- limma::topTable(fit, number = Inf, sort.by = "none")
tab$gene <- rownames(tab)
rownames(tab) <- NULL

ep_median <- apply(
  candidate_expr[, group == "EP", drop = FALSE],
  1, median, na.rm = TRUE
)
ctrl_median <- apply(
  candidate_expr[, group == "control", drop = FALSE],
  1, median, na.rm = TRUE
)
median_log2FC <- ep_median - ctrl_median

tab$median_log2FC <- median_log2FC[tab$gene]
tab$direction <- ifelse(tab$median_log2FC > 0, "up", "down")
tab$nominal_p_lt_0.05 <- tab$P.Value < 0.05
tab$BH_FDR_lt_0.05 <- tab$adj.P.Val < 0.05

out <- tab[, c(
  "gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val",
  "median_log2FC", "direction",
  "nominal_p_lt_0.05", "BH_FDR_lt_0.05"
)]

write.csv(
  out,
  file.path(RESULT_DIR, "GSE190451_candidate_validation.csv"),
  row.names = FALSE, quote = FALSE
)

supported_nominal <- out$gene[out$nominal_p_lt_0.05]
write.csv(
  data.frame(gene = supported_nominal),
  file.path(RESULT_DIR, "GSE190451_supported_candidates_nominal.csv"),
  row.names = FALSE, quote = FALSE
)

print(out)
message("Nominal p < 0.05: ",
        paste(supported_nominal, collapse = ", "))
