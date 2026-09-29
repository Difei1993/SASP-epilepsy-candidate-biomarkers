# 01_prepare_GSE190451.R
# Download and prepare the processed TPM matrix for GSE190451.
# GEO design: 3 non-epileptic controls and 3 TLE samples.

source("scripts/00_config.R")

url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190451/",
  "suppl/GSE190451_TLE_samples.txt.gz"
)
dest <- file.path(RAW_DIR, "GSE190451_TLE_samples.txt.gz")

if (!file.exists(dest)) {
  message("Downloading GSE190451 processed TPM matrix from GEO...")
  download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
}

x <- read.delim(gzfile(dest), header = TRUE, check.names = FALSE,
                stringsAsFactors = FALSE)

sample_cols <- grep("^(Normal|TLE)[_\\.-]?[123]$",
                    colnames(x), value = TRUE, ignore.case = TRUE)
if (length(sample_cols) != 6L) {
  sample_cols <- grep("Normal|TLE", colnames(x), value = TRUE, ignore.case = TRUE)
}
if (length(sample_cols) != 6L) {
  stop("Could not identify the six GSE190451 sample columns. Inspect colnames(x).")
}

candidate_gene_cols <- c(
  "gene", "Gene", "gene_name", "GeneName", "symbol", "Symbol",
  "gene_symbol", "GeneSymbol"
)
gene_col <- candidate_gene_cols[candidate_gene_cols %in% colnames(x)][1]

if (is.na(gene_col)) {
  non_sample <- setdiff(colnames(x), sample_cols)
  char_cols <- non_sample[vapply(
    x[non_sample],
    function(z) is.character(z) || is.factor(z),
    logical(1)
  )]
  gene_col <- if (length(char_cols)) char_cols[1] else non_sample[1]
}

expr <- x[, sample_cols, drop = FALSE]
expr[] <- lapply(expr, function(z) suppressWarnings(as.numeric(z)))
gene <- trimws(as.character(x[[gene_col]]))

keep <- !is.na(gene) & nzchar(gene)
expr <- expr[keep, , drop = FALSE]
gene <- gene[keep]

if (mean(grepl("^ENSG", gene)) > 0.5) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install AnnotationDbi and org.Hs.eg.db for Ensembl-to-symbol mapping.")
  }
  gene0 <- sub("\\..*$", "", gene)
  map <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = gene0, keytype = "ENSEMBL", column = "SYMBOL",
    multiVals = "first"
  )
  gene <- unname(map[gene0])
}

if (mean(grepl("^[0-9]+$", gene), na.rm = TRUE) > 0.5) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install AnnotationDbi and org.Hs.eg.db for Entrez-to-symbol mapping.")
  }
  map <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = gene, keytype = "ENTREZID", column = "SYMBOL",
    multiVals = "first"
  )
  gene <- unname(map[gene])
}

keep <- !is.na(gene) & nzchar(gene)
expr <- expr[keep, , drop = FALSE]
gene <- gene[keep]

expr$gene_symbol <- gene
expr <- aggregate(. ~ gene_symbol, data = expr, FUN = median, na.rm = TRUE)
rownames(expr) <- expr$gene_symbol
expr$gene_symbol <- NULL

expr_log2 <- log2(as.matrix(expr) + 1)
colnames(expr_log2) <- gsub("\\.", "_", colnames(expr_log2))
groups <- ifelse(
  grepl("^TLE", colnames(expr_log2), ignore.case = TRUE),
  "EP", "control"
)

write.csv(
  expr_log2,
  file.path(PROCESSED_DIR, "GSE190451_log2TPM.csv"),
  quote = FALSE
)
write.csv(
  data.frame(sample = colnames(expr_log2), group = groups),
  file.path(PROCESSED_DIR, "GSE190451_groups.csv"),
  row.names = FALSE, quote = FALSE
)

message("Prepared GSE190451: ", nrow(expr_log2),
        " genes x ", ncol(expr_log2), " samples.")
print(table(groups))
