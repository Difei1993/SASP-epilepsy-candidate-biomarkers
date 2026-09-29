# 03_prepare_GSE256068_metadata.R
# Retrieve public sample metadata for GSE256068 directly from GEO and extract
# sex, age, disease state, and seizure frequency when available.

source("scripts/00_config.R")

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("Install Bioconductor package 'GEOquery'.")
}

gse <- GEOquery::getGEO("GSE256068", GSEMatrix = TRUE, getGPL = FALSE)
if (is.list(gse)) gse <- gse[[1]]
pd <- Biobase::pData(gse)

char_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
if (!length(char_cols)) {
  stop("No characteristics_ch1 columns were found in GSE256068 metadata.")
}

char_text <- apply(pd[, char_cols, drop = FALSE], 1, function(z) {
  paste(z[!is.na(z)], collapse = "; ")
})

extract_field <- function(x, label_regex) {
  m <- regexpr(
    paste0("(?i)(?:", label_regex, ")\\s*:\\s*([^;]+)"),
    x, perl = TRUE
  )
  out <- rep(NA_character_, length(x))
  hit <- m > 0
  if (any(hit)) {
    tmp <- regmatches(x[hit], m[hit])
    out[hit] <- sub(
      paste0("(?i)^(?:", label_regex, ")\\s*:\\s*"),
      "", tmp, perl = TRUE
    )
  }
  trimws(out)
}

sex <- extract_field(char_text, "sex|gender")
age <- extract_field(char_text, "age")
seizure_frequency <- extract_field(
  char_text, "seizure frequency|seisure frequency"
)
disease_state <- extract_field(
  char_text, "disease state|diagnosis"
)
tissue <- extract_field(char_text, "tissue")

age_num <- suppressWarnings(as.numeric(gsub("[^0-9.]+", "", age)))
freq_num <- suppressWarnings(as.numeric(gsub("[^0-9.]+", "", seizure_frequency)))

group <- ifelse(
  grepl("control|non-epileptic|normal", disease_state, ignore.case = TRUE),
  "control", "EP"
)

meta <- data.frame(
  sample = rownames(pd),
  title = if ("title" %in% colnames(pd)) pd$title else rownames(pd),
  group = group,
  disease_state = disease_state,
  tissue = tissue,
  sex = sex,
  age = age_num,
  seizure_frequency = freq_num,
  stringsAsFactors = FALSE
)

write.csv(
  meta,
  file.path(PROCESSED_DIR, "GSE256068_metadata.csv"),
  row.names = FALSE, quote = FALSE
)

message("Wrote GSE256068 metadata: ", nrow(meta), " samples.")
print(table(meta$group, useNA = "ifany"))
