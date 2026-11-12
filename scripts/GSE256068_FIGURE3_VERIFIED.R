options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript build_figure3_code_verified.R <TPM.tsv.gz> <series_matrix.txt.gz> <out_dir> <boot_reps>")
}

tpm_file <- normalizePath(args[[1]], mustWork = TRUE)
series_file <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- args[[3]]
boot_reps <- as.integer(args[[4]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(rms)
  library(pROC)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

figure_contract <- c(
  "Core conclusion: A four-gene logistic model estimates the probability of epilepsy in the reconstructed GSE256068 discovery cohort.",
  "Archetype: quantitative grid with a full-width nomogram and two supporting validation panels.",
  "Backend: R only.",
  "Outcome: EP=1 and Control=0 in every panel.",
  "Prediction scale: inverse-logit probability of EP, never the untransformed linear predictor.",
  "Caveat: apparent/training-cohort performance; 121 EP and 12 controls; exploratory model only."
)
writeLines(figure_contract, file.path(out_dir, "Figure3_figure_contract.txt"))

parse_geo_series <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  parsed <- strsplit(lines[grepl("^!Sample_", lines)], "\t")
  clean <- function(x) sub('^"(.*)"$', '\\1', x)
  parsed <- lapply(parsed, clean)
  values <- function(key, occurrence = 1L) {
    hits <- parsed[vapply(parsed, function(x) identical(x[[1]], key), logical(1))]
    if (length(hits) < occurrence) stop("Missing GEO field: ", key)
    hits[[occurrence]][-1]
  }
  chars <- parsed[vapply(parsed, function(x) identical(x[[1]], "!Sample_characteristics_ch1"), logical(1))]
  out <- data.frame(
    GSM = values("!Sample_geo_accession"),
    sample_name = values("!Sample_title"),
    original_name = values("!Sample_description"),
    brain_region = values("!Sample_source_name_ch1"),
    stringsAsFactors = FALSE
  )
  for (block in chars) {
    vals <- block[-1]
    field <- tolower(trimws(sub(":.*$", "", vals[[1]])))
    out[[field]] <- trimws(sub("^[^:]+:", "", vals))
  }
  out
}

meta <- parse_geo_series(series_file)
tpm <- data.table::fread(tpm_file, data.table = FALSE, check.names = FALSE)
stopifnot(identical(names(tpm)[[1]], "GeneID"))
matrix_gsms <- names(tpm)[-1]

meta$prefix <- sub("_[0-9]+$", "", meta$sample_name)
meta$in_ncbi_tpm <- meta$GSM %in% matrix_gsms
meta$included_133 <- meta$in_ncbi_tpm & meta$prefix != "ControlHippocampus"
meta$analysis_group <- ifelse(meta$prefix == "ControlCortex", "Control",
                              ifelse(grepl("^Control", meta$prefix), "Excluded control", "EP"))

prefix_order <- c("TSC", "FCD2a", "FCD2b", "TLE-HS", "ControlCortex")
selected <- meta[meta$included_133, , drop = FALSE]
selected$prefix_rank <- match(selected$prefix, prefix_order)
selected$matrix_rank <- match(selected$GSM, matrix_gsms)
selected <- selected[order(selected$prefix_rank, selected$matrix_rank), , drop = FALSE]
stopifnot(nrow(selected) == 133L,
          sum(selected$analysis_group == "EP") == 121L,
          sum(selected$analysis_group == "Control") == 12L)

genes <- c("SERPINE1", "IGFBP4", "C3", "SPX")
mapping <- AnnotationDbi::select(org.Hs.eg.db, keys = genes,
                                 keytype = "SYMBOL", columns = "ENTREZID")
mapping <- mapping[!duplicated(mapping$SYMBOL), ]
mapping <- mapping[match(genes, mapping$SYMBOL), ]
row_idx <- match(mapping$ENTREZID, as.character(tpm$GeneID))
if (anyNA(row_idx)) stop("One or more model genes are absent from the official TPM matrix")

expr <- as.matrix(tpm[row_idx, matrix_gsms, drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- mapping$SYMBOL
expr <- log2(expr + 1)

model_data <- as.data.frame(t(expr[, selected$GSM, drop = FALSE]))
model_data <- model_data[, genes, drop = FALSE]
model_data$GSM <- rownames(model_data)
model_data$sample_name <- selected$sample_name[match(model_data$GSM, selected$GSM)]
model_data$group <- selected$analysis_group[match(model_data$GSM, selected$GSM)]
model_data$EP <- as.integer(model_data$group == "EP")
model_data <- model_data[, c("GSM", "sample_name", "group", "EP", genes)]
stopifnot(!anyNA(model_data))

fit_glm <- glm(EP ~ SERPINE1 + IGFBP4 + C3 + SPX,
               data = model_data, family = binomial(link = "logit"))
model_data$linear_predictor <- predict(fit_glm, type = "link")
model_data$predicted_probability_EP <- predict(fit_glm, type = "response")

dd <- rms::datadist(model_data[, c("EP", genes)])
options(datadist = "dd")
fit_lrm <- rms::lrm(EP ~ SERPINE1 + IGFBP4 + C3 + SPX,
                    data = model_data, x = TRUE, y = TRUE)

set.seed(3)
cal <- rms::calibrate(fit_lrm, method = "boot", B = boot_reps)
cal_columns <- c("predy", "calibrated.orig", "calibrated.corrected",
                 "index.orig", "training", "test", "optimism", "index.corrected", "n")
available_cal_columns <- intersect(cal_columns, colnames(cal))
cal_df <- data.frame(lapply(available_cal_columns, function(nm) as.numeric(cal[, nm])),
                     check.names = FALSE)
names(cal_df) <- available_cal_columns

net_benefit <- function(y, p, threshold) {
  positive <- p >= threshold
  tp <- sum(positive & y == 1L)
  fp <- sum(positive & y == 0L)
  tp / length(y) - fp / length(y) * threshold / (1 - threshold)
}

thresholds <- seq(0.01, 0.99, by = 0.01)
probabilities <- list(`Four-gene model` = model_data$predicted_probability_EP)
for (gene in genes) {
  fit_one <- glm(reformulate(gene, response = "EP"), data = model_data,
                 family = binomial(link = "logit"))
  probabilities[[gene]] <- predict(fit_one, type = "response")
}

dca_df <- do.call(rbind, lapply(names(probabilities), function(model_name) {
  data.frame(
    threshold = thresholds,
    net_benefit = vapply(thresholds, function(t) {
      net_benefit(model_data$EP, probabilities[[model_name]], t)
    }, numeric(1)),
    curve = model_name,
    stringsAsFactors = FALSE
  )
}))
prevalence <- mean(model_data$EP)
dca_df <- rbind(
  dca_df,
  data.frame(threshold = thresholds,
             net_benefit = prevalence - (1 - prevalence) * thresholds / (1 - thresholds),
             curve = "Treat all"),
  data.frame(threshold = thresholds, net_benefit = 0, curve = "Treat none")
)

roc_obj <- pROC::roc(model_data$EP, model_data$predicted_probability_EP,
                     levels = c(0, 1), direction = "<", quiet = TRUE)
auc_value <- as.numeric(pROC::auc(roc_obj))
brier <- mean((model_data$EP - model_data$predicted_probability_EP)^2)

coef_tab <- data.frame(
  term = rownames(summary(fit_glm)$coefficients),
  summary(fit_glm)$coefficients,
  row.names = NULL,
  check.names = FALSE
)
metrics <- data.frame(
  metric = c("n", "EP_n", "Control_n", "event_coding", "model_formula",
             "prediction_scale", "AUC_apparent", "Brier_score_apparent",
             "bootstrap_repetitions", "TPM_transform", "model_status"),
  value = c(nrow(model_data), sum(model_data$EP), sum(model_data$EP == 0L),
            "EP=1; Control=0", "EP ~ SERPINE1 + IGFBP4 + C3 + SPX",
            "Probability of EP (inverse-logit of linear predictor)",
            signif(auc_value, 8), signif(brier, 8), boot_reps,
            "log2(TPM + 1)", "Exploratory; apparent/training-cohort performance")
)

write.csv(model_data, file.path(out_dir, "Figure3_source_data.csv"), row.names = FALSE, quote = TRUE)
write.csv(coef_tab, file.path(out_dir, "Figure3_model_coefficients.csv"), row.names = FALSE, quote = TRUE)
write.csv(cal_df, file.path(out_dir, "Figure3_calibration_source_data.csv"), row.names = FALSE, quote = TRUE)
write.csv(dca_df, file.path(out_dir, "Figure3_DCA_source_data.csv"), row.names = FALSE, quote = TRUE)
write.csv(metrics, file.path(out_dir, "Figure3_model_metrics.csv"), row.names = FALSE, quote = TRUE)
capture.output(summary(fit_glm), file = file.path(out_dir, "Figure3_model_summary.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "Figure3_R_sessionInfo.txt"))

palette <- c(
  `Four-gene model` = "#1F4E79",
  SERPINE1 = "#C23B32",
  IGFBP4 = "#8C564B",
  C3 = "#D89C00",
  SPX = "#5B5FC7",
  `Treat all` = "#9A9A9A",
  `Treat none` = "#222222"
)

theme_pub <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5),
      plot.margin = margin(4, 5, 4, 5),
      panel.grid = element_blank()
    )
}

beta <- stats::coef(fit_glm)
gene_beta <- beta[genes]
gene_ranges <- lapply(genes, function(gene) range(model_data[[gene]], finite = TRUE))
names(gene_ranges) <- genes
contribution_ranges <- mapply(function(b, r) range(b * r), gene_beta, gene_ranges,
                              SIMPLIFY = FALSE)
point_scale <- max(vapply(contribution_ranges, function(x) diff(x), numeric(1))) / 100
gene_max_points <- vapply(contribution_ranges, function(x) diff(x) / point_scale, numeric(1))
total_max_points <- sum(gene_max_points)
base_lp <- unname(beta[["(Intercept)"]]) +
  sum(vapply(contribution_ranges, min, numeric(1)))

nom_rows <- data.frame(
  label = c("Points", genes, "Total points", "Predicted probability of epilepsy"),
  y = c(7, 6, 5, 4, 3, 2, 1),
  face = c("plain", rep("plain", length(genes)), "bold", "plain"),
  stringsAsFactors = FALSE
)

gene_ticks <- do.call(rbind, lapply(seq_along(genes), function(i) {
  gene <- genes[[i]]
  values <- pretty(gene_ranges[[gene]], n = 5)
  values <- values[values >= gene_ranges[[gene]][1] & values <= gene_ranges[[gene]][2]]
  contribution <- gene_beta[[gene]] * values
  min_contribution <- min(contribution_ranges[[gene]])
  points <- (contribution - min_contribution) / point_scale
  data.frame(gene = gene, value = values, points = points,
             x = points, y = 7 - i, stringsAsFactors = FALSE)
}))

gene_segments <- do.call(rbind, lapply(seq_along(genes), function(i) {
  gene <- genes[[i]]
  data.frame(gene = gene, x = 0, xend = gene_max_points[[gene]], y = 7 - i)
}))

point_ticks <- seq(0, 100, by = 10)
point_tick_df <- data.frame(x = point_ticks, label = point_ticks)
total_ticks <- pretty(c(0, total_max_points), n = 7)
total_ticks <- total_ticks[total_ticks >= 0 & total_ticks <= total_max_points]
total_x <- total_ticks / total_max_points * 100
total_tick_df <- data.frame(x = total_x, label = total_ticks)

prob_ticks <- c(0.01, 0.20, 0.80, 0.99)
prob_total <- (qlogis(prob_ticks) - base_lp) / point_scale
keep_prob <- is.finite(prob_total) & prob_total >= 0 & prob_total <= total_max_points
prob_ticks <- prob_ticks[keep_prob]
prob_x <- prob_total[keep_prob] / total_max_points * 100
prob_tick_df <- data.frame(x = prob_x,
                           label = format(prob_ticks, trim = TRUE, scientific = FALSE))

p_nom <- ggplot() +
  geom_segment(data = gene_segments,
               aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 0.42, colour = "black") +
  geom_segment(data = gene_ticks,
               aes(x = x, xend = x, y = y, yend = y - 0.10),
               linewidth = 0.35, colour = "black") +
  geom_text(data = gene_ticks,
            aes(x = x, y = y - 0.22, label = format(value, digits = 3, trim = TRUE)),
            size = 2.1, family = "Arial", vjust = 1) +
  geom_segment(aes(x = 0, xend = 100, y = 7, yend = 7), linewidth = 0.42) +
  geom_segment(data = point_tick_df,
               aes(x = x, xend = x, y = 7, yend = 6.90), linewidth = 0.35) +
  geom_text(data = point_tick_df, aes(x = x, y = 6.80, label = label),
            size = 2.05, family = "Arial", vjust = 1) +
  geom_segment(aes(x = 0, xend = 100, y = 2, yend = 2), linewidth = 0.42) +
  geom_segment(data = total_tick_df,
               aes(x = x, xend = x, y = 2, yend = 1.90), linewidth = 0.35) +
  geom_text(data = total_tick_df, aes(x = x, y = 1.80, label = label),
            size = 2.05, family = "Arial", vjust = 1) +
  geom_segment(aes(x = min(prob_x), xend = max(prob_x), y = 1, yend = 1), linewidth = 0.42) +
  geom_segment(data = prob_tick_df,
               aes(x = x, xend = x, y = 1, yend = 0.90), linewidth = 0.35) +
  geom_text(data = prob_tick_df, aes(x = x, y = 0.80, label = label),
            size = 2.05, family = "Arial", vjust = 1) +
  geom_text(data = nom_rows, aes(x = -4, y = y, label = label, fontface = face),
            hjust = 1, size = 2.5, family = "Arial") +
  coord_cartesian(xlim = c(-29, 103), ylim = c(0.45, 7.35), clip = "off", expand = FALSE) +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(4, 8, 3, 12))

p_cal <- ggplot(cal_df, aes(x = predy)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.45,
              linetype = "dashed", colour = "#333333") +
  geom_line(aes(y = calibrated.orig, colour = "Apparent"), linewidth = 0.65) +
  geom_line(aes(y = calibrated.corrected, colour = "Bias-corrected"), linewidth = 0.65) +
  scale_colour_manual(values = c(Apparent = "#1F4E79", `Bias-corrected` = "#D06B32")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "Predicted probability of epilepsy",
       y = "Observed probability of epilepsy") +
  annotate("text", x = 0.04, y = 0.96,
           label = paste0("n = 133; bootstrap B = ", boot_reps),
           hjust = 0, vjust = 1, size = 2.25, family = "Arial") +
  theme_pub() +
  theme(legend.position = c(0.74, 0.16),
        legend.background = element_rect(fill = "white", colour = NA))

dca_df$curve <- factor(dca_df$curve,
                       levels = c("Four-gene model", genes, "Treat all", "Treat none"))
p_dca <- ggplot(dca_df, aes(threshold, net_benefit, colour = curve, linetype = curve)) +
  geom_line(linewidth = 0.62) +
  scale_colour_manual(values = palette, drop = FALSE) +
  scale_linetype_manual(values = c(
    `Four-gene model` = "solid", SERPINE1 = "solid", IGFBP4 = "solid",
    C3 = "solid", SPX = "solid", `Treat all` = "longdash", `Treat none` = "dotted"
  ), drop = FALSE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.05, 1.0), expand = FALSE) +
  labs(x = "Threshold probability for epilepsy", y = "Net benefit") +
  theme_pub() +
  theme(legend.position = c(0.23, 0.38),
        legend.background = element_rect(fill = "white", colour = "#D9D9D9", linewidth = 0.25),
        legend.key.height = unit(3.5, "mm"))

figure <- (p_nom / (p_cal | p_dca)) +
  plot_layout(heights = c(1.18, 1), widths = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold", family = "Arial"))

# Static export contract for source preflight: width_mm = 183
width_mm <- 183
height_mm <- 170
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(file.path(out_dir, "Figure3_CODE_VERIFIED.svg"),
                 width = width_in, height = height_in)
print(figure)
dev.off()

grDevices::cairo_pdf(file.path(out_dir, "Figure3_CODE_VERIFIED.pdf"),
                     width = width_in, height = height_in, family = "Arial")
print(figure)
dev.off()

ragg::agg_tiff(file.path(out_dir, "Figure3_CODE_VERIFIED.tiff"),
               width = width_in, height = height_in, units = "in", res = 600,
               compression = "lzw", background = "white")
print(figure)
dev.off()

ragg::agg_png(file.path(out_dir, "Figure3_CODE_VERIFIED.png"),
              width = width_in, height = height_in, units = "in", res = 300,
              background = "white")
print(figure)
dev.off()

message("Figure 3 rebuilt from ", nrow(model_data), " samples (EP=",
        sum(model_data$EP), ", Control=", sum(model_data$EP == 0L), ").")
message("Apparent AUC=", signif(auc_value, 5), "; Brier=", signif(brier, 5), ".")
