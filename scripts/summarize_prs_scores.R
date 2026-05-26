#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: summarize_prs_scores.R <score_manifest.tsv> <out_dir>")
}

manifest <- read.table(args[[1]], header = TRUE, sep = "\t", quote = "",
                       comment.char = "", check.names = FALSE)
out_dir <- args[[2]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required <- c("method", "label", "path", "sample_col", "score_col", "threshold")
missing <- setdiff(required, names(manifest))
if (length(missing)) stop("manifest missing columns: ", paste(missing, collapse = ","))

read_any <- function(path) {
  if (grepl("[.]csv$", path, ignore.case = TRUE)) {
    read.csv(path, check.names = FALSE)
  } else {
    read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
               check.names = FALSE)
  }
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

records <- list()
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  x <- read_any(row$path)
  if (!row$sample_col %in% names(x)) stop("missing sample_col in ", row$path)
  if (!row$score_col %in% names(x)) stop("missing score_col in ", row$path)
  records[[i]] <- data.frame(
    method = row$method,
    label = row$label,
    threshold = row$threshold,
    score_col = paste(row$method, row$label, row$threshold, row$score_col, sep = "__"),
    SAMPLE = as.character(x[[row$sample_col]]),
    raw_score = as.numeric(x[[row$score_col]]),
    stringsAsFactors = FALSE
  )
}

long <- do.call(rbind, records)
long$z_score <- ave(long$raw_score, long$score_col, FUN = zscore)

samples <- unique(long$SAMPLE)
wide_from <- function(value_col) {
  out <- data.frame(SAMPLE = samples, stringsAsFactors = FALSE)
  for (col in unique(long$score_col)) {
    sub <- long[long$score_col == col, ]
    out[[col]] <- sub[[value_col]][match(samples, sub$SAMPLE)]
  }
  out
}

write.table(long, file.path(out_dir, "scores_long.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")
write.table(wide_from("raw_score"), file.path(out_dir, "scores_raw.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")
write.table(wide_from("z_score"), file.path(out_dir, "scores_z.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")
