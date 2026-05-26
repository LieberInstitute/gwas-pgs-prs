#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: compare_prs_methods.R <score_manifest.tsv> <out_dir>")
}

manifest <- read.table(args[[1]], header = TRUE, sep = "\t", quote = "",
                       comment.char = "", check.names = FALSE)
out_dir <- args[[2]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required <- c("method", "label", "path", "sample_col", "score_col", "threshold")
missing <- setdiff(required, names(manifest))
if (length(missing)) stop("manifest missing columns: ", paste(missing, collapse = ","))

read_any <- function(path) {
  if (grepl("[.]csv$", path, ignore.case = TRUE)) read.csv(path, check.names = FALSE) else
    read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
               check.names = FALSE)
}

fmt_thr <- function(x) ifelse(is.na(x) | x == ".", "NA", as.character(x))
safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  as.numeric(stats::cor(x[ok], y[ok], method = method))
}
zscore <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

records <- list()
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  x <- read_any(row$path)
  records[[i]] <- data.frame(
    method = row$method,
    label = row$label,
    threshold = suppressWarnings(as.numeric(row$threshold)),
    threshold_label = fmt_thr(row$threshold),
    score_col = paste(row$method, row$label, fmt_thr(row$threshold), row$score_col, sep = "__"),
    SAMPLE = as.character(x[[row$sample_col]]),
    raw_score = as.numeric(x[[row$score_col]]),
    stringsAsFactors = FALSE
  )
}

long <- do.call(rbind, records)
long$z_score <- ave(long$raw_score, long$score_col, FUN = zscore)

samples <- unique(long$SAMPLE)
score_cols <- unique(long$score_col)

inventory <- aggregate(SAMPLE ~ method + label, long, function(x) length(unique(x)))
names(inventory)[names(inventory) == "SAMPLE"] <- "sample_count"
write.table(inventory, file.path(out_dir, "score_inventory.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

availability <- unique(long[c("method", "label", "threshold_label", "score_col")])
write.table(availability, file.path(out_dir, "threshold_availability.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

pair_rows <- list()
idx <- 0L
for (label in unique(long$label)) {
  cols <- unique(long$score_col[long$label == label])
  if (length(cols) < 2) next
  for (i in seq_len(length(cols) - 1L)) {
    for (j in seq.int(i + 1L, length(cols))) {
      x <- long[long$score_col == cols[i], ]
      y <- long[long$score_col == cols[j], ]
      x <- x[match(samples, x$SAMPLE), ]
      y <- y[match(samples, y$SAMPLE), ]
      idx <- idx + 1L
      pair_rows[[idx]] <- data.frame(
        label = label,
        score_x = cols[i],
        method_x = unique(x$method),
        threshold_x = unique(x$threshold_label),
        score_y = cols[j],
        method_y = unique(y$method),
        threshold_y = unique(y$threshold_label),
        n = sum(is.finite(x$raw_score) & is.finite(y$raw_score)),
        pearson = safe_cor(x$raw_score, y$raw_score, "pearson"),
        spearman = safe_cor(x$raw_score, y$raw_score, "spearman"),
        stringsAsFactors = FALSE
      )
    }
  }
}

pairwise <- if (length(pair_rows)) do.call(rbind, pair_rows) else data.frame()
write.table(pairwise, file.path(out_dir, "all_pairwise_correlations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

wide <- data.frame(SAMPLE = samples, stringsAsFactors = FALSE)
for (col in score_cols) {
  sub <- long[long$score_col == col, ]
  wide[[col]] <- sub$raw_score[match(samples, sub$SAMPLE)]
}
write.table(wide, file.path(out_dir, "combined_scores_raw.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

wide_z <- data.frame(SAMPLE = samples, stringsAsFactors = FALSE)
for (col in score_cols) {
  sub <- long[long$score_col == col, ]
  wide_z[[col]] <- sub$z_score[match(samples, sub$SAMPLE)]
}
write.table(wide_z, file.path(out_dir, "combined_scores_z.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")
