# 03_clinical_association_GSE256068.R
# Exploratory association of candidate-gene expression with available clinical
# variables in epilepsy cases from GSE256068.
#
# Required local inputs:
#   data/processed/GSE256068_log2TPM.csv
#   data/processed/GSE256068_metadata.csv
#
# Metadata columns:
#   sample, group, sex, age, seizure_frequency

source("scripts/00_config.R")

expr_file <- file.path(PROCESSED_DIR, "GSE256068_log2TPM.csv")
meta_file <- file.path(PROCESSED_DIR, "GSE256068_metadata.csv")

if (!file.exists(expr_file) || !file.exists(meta_file)) {
  stop(
    "Missing GSE256068 processed expression/metadata files. ",
    "Prepare them from the original project files before running this script."
  )
}

expr <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)

required <- c("sample", "group", "sex", "age", "seizure_frequency")
missing_cols <- setdiff(required, colnames(meta))
if (length(missing_cols) > 0L) {
  stop("Missing metadata columns: ",
       paste(missing_cols, collapse = ", "))
}

ep_meta <- meta[
  tolower(meta$group) %in% c("ep", "epilepsy", "tle", "case"),
  , drop = FALSE
]
common_samples <- intersect(ep_meta$sample, colnames(expr))
ep_meta <- ep_meta[match(common_samples, ep_meta$sample), , drop = FALSE]
ep_expr <- expr[, common_samples, drop = FALSE]

candidate_file <- file.path(
  RESULT_DIR, "GSE190451_supported_candidates_nominal.csv"
)
if (file.exists(candidate_file)) {
  candidates <- read.csv(candidate_file, stringsAsFactors = FALSE)$gene
} else {
  candidates <- INITIAL_CANDIDATES
  warning("Validation candidate file not found; using five initial candidates.")
}
candidates <- intersect(candidates, rownames(ep_expr))

sex_results <- do.call(rbind, lapply(candidates, function(g) {
  d <- data.frame(
    expression = as.numeric(ep_expr[g, ]),
    sex = ep_meta$sex
  )
  d <- d[complete.cases(d) & nzchar(d$sex), ]
  if (length(unique(d$sex)) != 2L) {
    return(data.frame(
      gene = g, variable = "sex",
      statistic = NA, p_value = NA
    ))
  }
  wt <- wilcox.test(expression ~ sex, data = d, exact = FALSE)
  data.frame(
    gene = g, variable = "sex",
    statistic = unname(wt$statistic),
    p_value = wt$p.value
  )
}))

corr_one <- function(g, variable) {
  x <- suppressWarnings(as.numeric(ep_meta[[variable]]))
  y <- as.numeric(ep_expr[g, ])
  keep <- complete.cases(x, y)
  if (sum(keep) < 4L) {
    return(data.frame(
      gene = g, variable = variable,
      rho = NA, p_value = NA, n = sum(keep)
    ))
  }
  ct <- suppressWarnings(
    cor.test(x[keep], y[keep], method = "spearman", exact = FALSE)
  )
  data.frame(
    gene = g, variable = variable,
    rho = unname(ct$estimate),
    p_value = ct$p.value,
    n = sum(keep)
  )
}

continuous_results <- do.call(
  rbind,
  c(
    lapply(candidates, corr_one, variable = "age"),
    lapply(candidates, corr_one, variable = "seizure_frequency")
  )
)

write.csv(
  sex_results,
  file.path(RESULT_DIR, "clinical_association_sex.csv"),
  row.names = FALSE, quote = FALSE
)
write.csv(
  continuous_results,
  file.path(RESULT_DIR, "clinical_association_continuous.csv"),
  row.names = FALSE, quote = FALSE
)

print(sex_results)
print(continuous_results)
