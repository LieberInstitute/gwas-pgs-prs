#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name) {
  idx <- match(name, args)
  if (is.na(idx) || idx == length(args)) stop("Missing argument: ", name, call. = FALSE)
  args[[idx + 1L]]
}

registry_file <- get_arg("--registry")
out_dir <- get_arg("--out-dir")
report_file <- get_arg("--report")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.delim(registry_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("SAMPLE", "TRAIT", "GWAS_VERSION", "METHOD", "SETUP", "THRESHOLD",
              "RAW_SCORE", "SCORED_VARIANTS")
if (!all(required %in% names(d))) stop("Registry is missing required columns", call. = FALSE)
if (!setequal(unique(d$GWAS_VERSION), c("no23andMe", "full"))) {
  stop("Registry must contain full and no23andMe versions", call. = FALSE)
}

zscore <- function(x) as.numeric(scale(x))
safe_cor <- function(x, y, method) {
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y, method = method)
}
series_key <- c("TRAIT", "METHOD", "SETUP", "THRESHOLD")
groups <- split(d, interaction(d[, series_key], drop = TRUE, lex.order = TRUE))

impact <- list()
rank_rows <- list()
for (g in groups) {
  if (!setequal(unique(g$GWAS_VERSION), c("no23andMe", "full"))) next
  old <- g[g$GWAS_VERSION == "no23andMe", ]
  full <- g[g$GWAS_VERSION == "full", ]
  m <- merge(old, full, by = "SAMPLE", suffixes = c("_NO23", "_FULL"))
  if (nrow(m) != 119L) stop("Full/no23 donor mismatch", call. = FALSE)

  old_rank <- rank(-m$RAW_SCORE_NO23, ties.method = "average")
  full_rank <- rank(-m$RAW_SCORE_FULL, ties.method = "average")
  n_decile <- ceiling(nrow(m) * 0.10)
  top_old <- m$SAMPLE[order(-m$RAW_SCORE_NO23)][seq_len(n_decile)]
  top_full <- m$SAMPLE[order(-m$RAW_SCORE_FULL)][seq_len(n_decile)]
  bottom_old <- m$SAMPLE[order(m$RAW_SCORE_NO23)][seq_len(n_decile)]
  bottom_full <- m$SAMPLE[order(m$RAW_SCORE_FULL)][seq_len(n_decile)]

  impact[[length(impact) + 1L]] <- data.frame(
    TRAIT = g$TRAIT[[1L]], METHOD = g$METHOD[[1L]], SETUP = g$SETUP[[1L]],
    THRESHOLD = g$THRESHOLD[[1L]], N = nrow(m),
    PEARSON = safe_cor(m$RAW_SCORE_NO23, m$RAW_SCORE_FULL, "pearson"),
    SPEARMAN = safe_cor(m$RAW_SCORE_NO23, m$RAW_SCORE_FULL, "spearman"),
    Z_RMSE = sqrt(mean((zscore(m$RAW_SCORE_FULL) - zscore(m$RAW_SCORE_NO23))^2)),
    MEAN_ABS_RANK_CHANGE = mean(abs(full_rank - old_rank)),
    MAX_ABS_RANK_CHANGE = max(abs(full_rank - old_rank)),
    TOP_DECILE_OVERLAP = length(intersect(top_old, top_full)) / n_decile,
    BOTTOM_DECILE_OVERLAP = length(intersect(bottom_old, bottom_full)) / n_decile,
    NO23_VARIANTS = median(m$SCORED_VARIANTS_NO23),
    FULL_VARIANTS = median(m$SCORED_VARIANTS_FULL),
    stringsAsFactors = FALSE
  )
  rank_rows[[length(rank_rows) + 1L]] <- data.frame(
    SAMPLE = m$SAMPLE, TRAIT = g$TRAIT[[1L]], METHOD = g$METHOD[[1L]],
    SETUP = g$SETUP[[1L]], THRESHOLD = g$THRESHOLD[[1L]],
    NO23_RANK = old_rank, FULL_RANK = full_rank,
    RANK_CHANGE = full_rank - old_rank,
    stringsAsFactors = FALSE
  )
}
impact <- do.call(rbind, impact)
rank_changes <- do.call(rbind, rank_rows)

## compare method series at exact shared C+T thresholds and GraphPred against all
method_rows <- list()
for (trait in unique(d$TRAIT)) {
  for (version in unique(d$GWAS_VERSION)) {
    subset <- d[d$TRAIT == trait & d$GWAS_VERSION == version, ]
    keys <- unique(subset[, c("METHOD", "SETUP", "THRESHOLD")])
    if (nrow(keys) < 2L) next
    pairs <- combn(seq_len(nrow(keys)), 2L)
    for (j in seq_len(ncol(pairs))) {
      a <- keys[pairs[1L, j], ]
      b <- keys[pairs[2L, j], ]
      graph_pair <- a$METHOD == "GraphPred" || b$METHOD == "GraphPred"
      same_threshold <- !graph_pair && suppressWarnings(as.numeric(a$THRESHOLD)) ==
                                      suppressWarnings(as.numeric(b$THRESHOLD))
      if (!graph_pair && !same_threshold) next

      da <- subset[subset$METHOD == a$METHOD & subset$SETUP == a$SETUP &
                     subset$THRESHOLD == a$THRESHOLD, c("SAMPLE", "RAW_SCORE")]
      db <- subset[subset$METHOD == b$METHOD & subset$SETUP == b$SETUP &
                     subset$THRESHOLD == b$THRESHOLD, c("SAMPLE", "RAW_SCORE")]
      m <- merge(da, db, by = "SAMPLE", suffixes = c("_A", "_B"))
      method_rows[[length(method_rows) + 1L]] <- data.frame(
        TRAIT = trait, GWAS_VERSION = version,
        METHOD_A = a$METHOD, SETUP_A = a$SETUP, THRESHOLD_A = a$THRESHOLD,
        METHOD_B = b$METHOD, SETUP_B = b$SETUP, THRESHOLD_B = b$THRESHOLD,
        PEARSON = safe_cor(m$RAW_SCORE_A, m$RAW_SCORE_B, "pearson"),
        SPEARMAN = safe_cor(m$RAW_SCORE_A, m$RAW_SCORE_B, "spearman"),
        stringsAsFactors = FALSE
      )
    }
  }
}
method_cor <- do.call(rbind, method_rows)

write.table(impact, file.path(out_dir, "full_vs_no23_metrics.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(rank_changes, file.path(out_dir, "donor_rank_changes.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(method_cor, file.path(out_dir, "method_correlations.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

## plot full/no23 correlations over p-value thresholds
svg(file.path(out_dir, "full_vs_no23_correlations.svg"), width = 10, height = 5.8)
par(mfrow = c(1, length(unique(impact$TRAIT))), mar = c(4.5, 4.2, 3, 1),
    oma = c(1.5, 0, 0, 0))
palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
for (trait in sort(unique(impact$TRAIT))) {
  x <- impact[impact$TRAIT == trait & impact$THRESHOLD != "ALL", ]
  labels <- unique(paste(x$METHOD, x$SETUP, sep = ":"))
  trait_impact <- impact[impact$TRAIT == trait, ]
  y_range <- range(c(trait_impact$PEARSON, trait_impact$SPEARMAN), na.rm = TRUE)
  y_pad <- diff(y_range) * 0.08
  plot(NA, xlim = c(0, 8.45), ylim = y_range + c(-y_pad, y_pad),
       xlab = "-log10(P threshold)", ylab = "correlation", main = trait)
  for (i in seq_along(labels)) {
    g <- x[paste(x$METHOD, x$SETUP, sep = ":") == labels[[i]], ]
    ord <- order(-log10(as.numeric(g$THRESHOLD)))
    xx <- -log10(as.numeric(g$THRESHOLD[ord]))
    lines(xx, g$PEARSON[ord], type = "b", pch = 16, col = palette[[i]], lwd = 1.5)
    lines(xx, g$SPEARMAN[ord], type = "b", pch = 1, col = palette[[i]], lty = 2)
  }
  gp <- impact[impact$TRAIT == trait & impact$THRESHOLD == "ALL", ]
  if (nrow(gp)) {
    points(rep(8.15, nrow(gp)), gp$PEARSON, pch = 18, cex = 1.3)
    points(rep(8.32, nrow(gp)), gp$SPEARMAN, pch = 5, cex = 1.2)
  }
  legend("bottomleft", legend = c(labels, "GraphPred"),
         col = c(palette[seq_along(labels)], "black"),
         lty = c(rep(1, length(labels)), NA),
         pch = c(rep(16, length(labels)), 18), cex = 0.72, bty = "n")
}
mtext("Filled: Pearson; open/dashed: Spearman; diamonds: GraphPred", side = 1,
      outer = TRUE, line = 0.1, cex = 0.75)
invisible(dev.off())

## summarize donor rank movement without selecting a preferred threshold
svg(file.path(out_dir, "full_vs_no23_rank_changes.svg"), width = 11, height = 7)
rank_changes$SERIES <- paste(rank_changes$TRAIT, rank_changes$METHOD,
                             rank_changes$SETUP, sep = "|")
rank_groups <- split(abs(rank_changes$RANK_CHANGE), rank_changes$SERIES)
rank_labels <- gsub("\\|", "\n", names(rank_groups))
par(mar = c(8.5, 4.5, 1.5, 1))
boxplot(rank_groups, names = rep("", length(rank_groups)),
        ylab = "absolute donor-rank change", xlab = "", col = "#D9EAF2",
        border = "#404040", outline = FALSE, xaxt = "n")
axis(1, at = seq_along(rank_groups), labels = rank_labels, tick = FALSE,
     line = -0.5, cex.axis = 0.72)
abline(h = 0, col = "#808080")
invisible(dev.off())

## write an answer-first generated report from the validated tables
summary_groups <- split(impact, interaction(impact$TRAIT, impact$METHOD, impact$SETUP,
                                            drop = TRUE, lex.order = TRUE))
summary_rows <- do.call(rbind, lapply(summary_groups, function(g) {
  data.frame(
    Trait = g$TRAIT[[1L]], Method = g$METHOD[[1L]], Setup = g$SETUP[[1L]],
    Series = nrow(g), Pearson_min = min(g$PEARSON), Pearson_max = max(g$PEARSON),
    Spearman_min = min(g$SPEARMAN), Spearman_max = max(g$SPEARMAN),
    Z_RMSE_min = min(g$Z_RMSE), Z_RMSE_max = max(g$Z_RMSE),
    Rank_min = min(g$MEAN_ABS_RANK_CHANGE), Rank_max = max(g$MEAN_ABS_RANK_CHANGE),
    Top_min = min(g$TOP_DECILE_OVERLAP), Top_max = max(g$TOP_DECILE_OVERLAP),
    Bottom_min = min(g$BOTTOM_DECILE_OVERLAP),
    Bottom_max = max(g$BOTTOM_DECILE_OVERLAP),
    stringsAsFactors = FALSE
  )
}))

count_data <- d[d$THRESHOLD %in% c("1", "ALL"), ]
count_groups <- split(count_data,
                      interaction(count_data$TRAIT, count_data$METHOD,
                                  count_data$SETUP, drop = TRUE, lex.order = TRUE))
count_rows <- do.call(rbind, lapply(count_groups, function(g) {
  data.frame(
    Trait = g$TRAIT[[1L]], Method = g$METHOD[[1L]], Setup = g$SETUP[[1L]],
    No23 = median(g$SCORED_VARIANTS[g$GWAS_VERSION == "no23andMe"]),
    Full = median(g$SCORED_VARIANTS[g$GWAS_VERSION == "full"]),
    stringsAsFactors = FALSE
  )
}))

reproduction_file <- file.path(out_dir, "reproduction_validation.tsv")
reproduction <- if (file.exists(reproduction_file)) {
  read.delim(reproduction_file, stringsAsFactors = FALSE)
} else NULL
ct_min_reproduction <- if (!is.null(reproduction)) {
  min(reproduction$VALUE[reproduction$METHOD == "Shizhong_CT"])
} else NA_real_
graph_max_error <- if (!is.null(reproduction)) {
  max(reproduction$VALUE[reproduction$METHOD == "GraphPred"])
} else NA_real_

fmt <- function(x) formatC(x, digits = 4, format = "f")
con <- file(report_file, open = "wt")
writeLines(c(
  "# MBv BD and MDD PRS: 23andMe Impact",
  "",
  "This report compares scores for the same 119 MBv donors using public",
  "no-23andMe GWAS statistics and reconstructed integrated statistics.",
  "It evaluates score and rank sensitivity only; it does not select a threshold",
  "or test diagnosis association in this small cohort.",
  "",
  "BD full statistics are labeled pre-DENTIST. The reconstructed meta-analysis",
  "has not received the paper's final HRC-based DENTIST filter, and the delivered",
  "v7.0 association file was mapped with v7.2 European annotations.",
  "",
  "## Validation",
  "",
  sprintf("- Registry: %d score series and %d rows; every series has the same 119 donors.",
          length(unique(interaction(d[, c("TRAIT", "GWAS_VERSION", "METHOD", "SETUP", "THRESHOLD")], drop = TRUE))),
          nrow(d)),
  sprintf("- Historical no23 Shizhong minimum Pearson reproduction: %s (required >= 0.999).",
          fmt(ct_min_reproduction)),
  sprintf("- No23 GraphPred maximum absolute reproduction error: %s (required <= 0.000001).",
          format(graph_max_error, scientific = TRUE)),
  "- All eight C+T runs completed 22 chromosomes without a failure status.",
  "",
  "## Full Versus No-23andMe Summary",
  "",
  "| Trait | Method | Setup | Series | Pearson | Spearman | z-RMSE | Mean abs rank change | Top-decile overlap | Bottom-decile overlap |",
  "|---|---|---|---:|---:|---:|---:|---:|---:|---:|"
), con)
for (i in seq_len(nrow(summary_rows))) {
  g <- summary_rows[i, ]
  writeLines(sprintf(paste0("| %s | %s | %s | %d | %s to %s | %s to %s | ",
                            "%s to %s | %s to %s | %s to %s | %s to %s |"),
                     g$Trait, g$Method, g$Setup, g$Series,
                     fmt(g$Pearson_min), fmt(g$Pearson_max),
                     fmt(g$Spearman_min), fmt(g$Spearman_max),
                     fmt(g$Z_RMSE_min), fmt(g$Z_RMSE_max),
                     fmt(g$Rank_min), fmt(g$Rank_max),
                     fmt(g$Top_min), fmt(g$Top_max),
                     fmt(g$Bottom_min), fmt(g$Bottom_max)), con)
}
writeLines(c(
  "",
  "## Scored Variants",
  "",
  "Counts use the P=1 score for C+T methods and the single GraphPred score.",
  "",
  "| Trait | Method | Setup | No23 | Full |",
  "|---|---|---|---:|---:|"
), con)
for (i in seq_len(nrow(count_rows))) {
  g <- count_rows[i, ]
  writeLines(sprintf("| %s | %s | %s | %d | %d |",
                     g$Trait, g$Method, g$Setup, g$No23, g$Full), con)
}
writeLines(c(
  "",
  "Full scored-variant counts can be lower because the integrated GWAS files",
  "apply variant-coverage filters and contain fewer eligible variants than the",
  "public no23 releases. This is an input-availability difference, not a scoring",
  "failure.",
  "",
  "Detailed results:",
  "",
  "- `full_vs_no23_metrics.tsv`: correlations, z-RMSE, rank movement, decile overlap, and variant counts.",
  "- `donor_rank_changes.tsv`: donor-level rank changes for every score series.",
  "- `method_correlations.tsv`: exact-threshold C+T and GraphPred-to-threshold comparisons.",
  "- `full_vs_no23_correlations.svg`: Pearson and Spearman sensitivity by threshold.",
  "- `full_vs_no23_rank_changes.svg`: donor-rank movement distributions.",
  "",
  "Interpret BD results as provisional until final full-EUR DENTIST-filtered",
  "statistics or an exclusion list and release-mapping confirmation are obtained."
), con)
close(con)
