# 00_config.R
# Portable project configuration for the revised analysis.
# Run from the repository root.

options(stringsAsFactors = FALSE)

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
DATA_DIR <- file.path(PROJECT_ROOT, "data")
RAW_DIR <- file.path(DATA_DIR, "raw")
PROCESSED_DIR <- file.path(DATA_DIR, "processed")
RESULT_DIR <- file.path(PROJECT_ROOT, "results")

for (d in c(DATA_DIR, RAW_DIR, PROCESSED_DIR, RESULT_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Five genes entering cross-dataset assessment after feature prioritization.
INITIAL_CANDIDATES <- c("SERPINE1", "CCL2", "IGFBP4", "C3", "SPX")

SEED_LASSO <- 8L
SEED_SVMRFE <- 15L
SEED_GENERAL <- 20260929L

message("PROJECT_ROOT: ", PROJECT_ROOT)
message("INITIAL_CANDIDATES: ", paste(INITIAL_CANDIDATES, collapse = ", "))
