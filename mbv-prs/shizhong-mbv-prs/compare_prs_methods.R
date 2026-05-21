#!/usr/bin/env Rscript

## compare Shizhong C+T, PRSice C+T, and bcftools GraphPred PRS outputs

options(stringsAsFactors = FALSE)

out_dir <- "prs_method_comparison"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

input_files <- list(
  demog = "../MBv_demographics_n119.tab",
  graphpred = "../scores_bcftools/mbv_graphpred_scores.tsv",
  prsice = "../prsice/out/mbv_prsice_scores.tsv",
  shizhong_bpd = "PRS_Bipolar.csv",
  shizhong_mdd = "PRS_MDD.csv",
  shizhong_sczd = "PRS_SCZ.csv"
)

for (f in unlist(input_files)) {
  if (!file.exists(f)) stop("Missing required input file: ", f)
}

read_tsv <- function(path) {
  read.table(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             check.names = FALSE)
}

read_csv <- function(path) {
  read.csv(path, header = TRUE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA")
}

fmt_thr <- function(x) {
  ## stable labels used for matching and readable output
  vapply(x, function(v) {
    if (is.na(v)) return("NA")
    if (v < 1e-3) return(format(v, scientific = TRUE, digits = 1))
    format(v, scientific = FALSE, trim = TRUE)
  }, character(1))
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  as.numeric(cor(x[ok], y[ok], method = method))
}

rmse <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) == 0) return(NA_real_)
  sqrt(mean((x[ok] - y[ok]) ^ 2))
}

sample_check <- function(reference, observed, label) {
  missing <- setdiff(reference, observed)
  extra <- setdiff(observed, reference)
  if (length(missing) || length(extra)) {
    stop(label, " sample mismatch. Missing: ", paste(missing, collapse = ","),
         " Extra: ", paste(extra, collapse = ","))
  }
}

demog <- read_tsv(input_files$demog)
if (!"BrNum" %in% names(demog)) stop("Demographics file lacks BrNum column")
names(demog)[names(demog) == "BrNum"] <- "SAMPLE"

graphpred <- read_tsv(input_files$graphpred)
prsice <- read_tsv(input_files$prsice)
shizhong_raw <- list(
  BPD = read_csv(input_files$shizhong_bpd),
  MDD = read_csv(input_files$shizhong_mdd),
  SCZD = read_csv(input_files$shizhong_sczd)
)

for (d in names(shizhong_raw)) {
  names(shizhong_raw[[d]])[1] <- "SAMPLE"
}

reference_samples <- demog$SAMPLE
sample_check(reference_samples, graphpred$SAMPLE, "GraphPred")
sample_check(reference_samples, prsice$SAMPLE, "PRSice")
for (d in names(shizhong_raw)) {
  sample_check(reference_samples, shizhong_raw[[d]]$SAMPLE, paste("Shizhong", d))
}

## keep all data in demographic order
graphpred <- graphpred[match(reference_samples, graphpred$SAMPLE), , drop = FALSE]
prsice <- prsice[match(reference_samples, prsice$SAMPLE), , drop = FALSE]
for (d in names(shizhong_raw)) {
  shizhong_raw[[d]] <- shizhong_raw[[d]][match(reference_samples, shizhong_raw[[d]]$SAMPLE), , drop = FALSE]
}

score_records <- list()
record_idx <- 0

add_record <- function(score_records, source_file, method, disorder, threshold,
                       score_col, sample, score) {
  ## each record contains one named score vector plus enough metadata for joins
  data.frame(
    source_file = source_file,
    method = method,
    disorder = disorder,
    threshold = threshold,
    threshold_label = fmt_thr(threshold),
    score_col = score_col,
    SAMPLE = sample,
    raw_score = as.numeric(score),
    z_score = zscore(as.numeric(score)),
    stringsAsFactors = FALSE
  )
}

## Shizhong thresholded score records
shizhong_files <- c(BPD = input_files$shizhong_bpd,
                   MDD = input_files$shizhong_mdd,
                   SCZD = input_files$shizhong_sczd)
for (d in names(shizhong_raw)) {
  x <- shizhong_raw[[d]]
  for (col in setdiff(names(x), "SAMPLE")) {
    threshold <- as.numeric(sub("^p cutoff=", "", col))
    record_idx <- record_idx + 1
    score_records[[record_idx]] <- add_record(
      score_records, shizhong_files[[d]], "Shizhong_C+T", d, threshold,
      paste0("Shizhong_", d, "_Pt_", fmt_thr(threshold)), x$SAMPLE, x[[col]]
    )
  }
}

## PRSice thresholded score records
prsice_cols <- setdiff(names(prsice), "SAMPLE")
for (col in prsice_cols) {
  m <- regexec("^PRSice_(BPD|MDD|SCZD)_Pt_(.+)$", col)
  parts <- regmatches(col, m)[[1]]
  if (length(parts) != 3) stop("Cannot parse PRSice column: ", col)
  d <- parts[2]
  threshold <- as.numeric(parts[3])
  record_idx <- record_idx + 1
  score_records[[record_idx]] <- add_record(
    score_records, input_files$prsice, "PRSice_C+T", d, threshold,
    col, prsice$SAMPLE, prsice[[col]]
  )
}

## GraphPred has one LD-aware score per disorder and no p-value threshold
graphpred_map <- c(
  BPD = "BPD_BIP_2024.EUR_pgs_a0.5_b5e-08",
  MDD = "MDD_MDD_2025_pgs_a0.5_b2e-08",
  SCZD = "SCZD_SCZ_2022.EUR_pgs_a0.5_b2e-07"
)
for (d in names(graphpred_map)) {
  col <- graphpred_map[[d]]
  if (!col %in% names(graphpred)) stop("Missing GraphPred column: ", col)
  record_idx <- record_idx + 1
  score_records[[record_idx]] <- add_record(
    score_records, input_files$graphpred, "GraphPred_bcftools", d, NA_real_,
    paste0("GraphPred_", d), graphpred$SAMPLE, graphpred[[col]]
  )
}

long_scores <- do.call(rbind, score_records)

## wide tables preserve one row per sample for easy downstream plotting
score_wide <- function(value_col) {
  pieces <- split(long_scores, long_scores$score_col)
  out <- data.frame(SAMPLE = reference_samples, stringsAsFactors = FALSE)
  for (nm in names(pieces)) {
    p <- pieces[[nm]]
    out[[nm]] <- p[[value_col]][match(reference_samples, p$SAMPLE)]
  }
  out
}

combined_raw <- merge(demog, score_wide("raw_score"), by = "SAMPLE", sort = FALSE)
combined_raw <- combined_raw[match(reference_samples, combined_raw$SAMPLE), , drop = FALSE]
combined_z <- merge(demog, score_wide("z_score"), by = "SAMPLE", sort = FALSE)
combined_z <- combined_z[match(reference_samples, combined_z$SAMPLE), , drop = FALSE]

write_tsv(combined_raw, file.path(out_dir, "combined_scores_raw.tsv"))
write_tsv(combined_z, file.path(out_dir, "combined_scores_z.tsv"))

## input inventory at the method/disorder level
inventory_rows <- list()
inv_idx <- 0
for (method in unique(long_scores$method)) {
  for (d in unique(long_scores$disorder[long_scores$method == method])) {
    sub <- long_scores[long_scores$method == method & long_scores$disorder == d, ]
    cols <- unique(sub$score_col)
    inv_idx <- inv_idx + 1
    inventory_rows[[inv_idx]] <- data.frame(
      source_file = paste(unique(sub$source_file), collapse = ";"),
      method = method,
      disorder = d,
      threshold_count = length(unique(sub$threshold_label)),
      sample_count = length(unique(sub$SAMPLE)),
      score_columns = paste(cols, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}
score_inventory <- do.call(rbind, inventory_rows)
write_tsv(score_inventory, file.path(out_dir, "score_inventory.tsv"))

## threshold availability table focuses on C+T methods
availability_rows <- list()
avail_idx <- 0
for (d in c("BPD", "MDD", "SCZD")) {
  sh <- sort(unique(long_scores$threshold[long_scores$method == "Shizhong_C+T" &
                                            long_scores$disorder == d]))
  pr <- sort(unique(long_scores$threshold[long_scores$method == "PRSice_C+T" &
                                            long_scores$disorder == d]))
  all_thr <- sort(unique(c(sh, pr)))
  for (thr in all_thr) {
    avail_idx <- avail_idx + 1
    availability_rows[[avail_idx]] <- data.frame(
      disorder = d,
      threshold = thr,
      threshold_label = fmt_thr(thr),
      shizhong_available = thr %in% sh,
      prsice_available = thr %in% pr,
      exact_match_available = (thr %in% sh) && (thr %in% pr),
      stringsAsFactors = FALSE
    )
  }
}
threshold_availability <- do.call(rbind, availability_rows)
write_tsv(threshold_availability, file.path(out_dir, "threshold_availability.tsv"))

get_vector <- function(method, disorder, threshold = NA_real_) {
  sub <- long_scores[long_scores$method == method & long_scores$disorder == disorder, ]
  if (is.na(threshold)) {
    sub <- sub[is.na(sub$threshold), ]
  } else {
    sub <- sub[!is.na(sub$threshold) & sub$threshold == threshold, ]
  }
  if (!nrow(sub)) return(NULL)
  sub[match(reference_samples, sub$SAMPLE), ]
}

cor_row <- function(x, y, label_x, label_y, disorder, threshold = NA_real_) {
  data.frame(
    disorder = disorder,
    threshold = threshold,
    threshold_label = fmt_thr(threshold),
    score_x = label_x,
    score_y = label_y,
    n = sum(is.finite(x$raw_score) & is.finite(y$raw_score)),
    pearson = safe_cor(x$raw_score, y$raw_score, "pearson"),
    spearman = safe_cor(x$raw_score, y$raw_score, "spearman"),
    z_rmse = rmse(x$z_score, y$z_score),
    stringsAsFactors = FALSE
  )
}

## exact-threshold Shizhong vs PRSice correlations
matched_rows <- list()
match_idx <- 0
for (d in c("BPD", "MDD", "SCZD")) {
  shared <- threshold_availability$threshold[
    threshold_availability$disorder == d & threshold_availability$exact_match_available
  ]
  for (thr in shared) {
    sh <- get_vector("Shizhong_C+T", d, thr)
    pr <- get_vector("PRSice_C+T", d, thr)
    match_idx <- match_idx + 1
    matched_rows[[match_idx]] <- cor_row(
      pr, sh, paste0("PRSice_", d, "_Pt_", fmt_thr(thr)),
      paste0("Shizhong_", d, "_Pt_", fmt_thr(thr)), d, thr
    )
  }
}
matched_threshold_correlations <- do.call(rbind, matched_rows)
write_tsv(matched_threshold_correlations,
          file.path(out_dir, "matched_threshold_correlations.tsv"))

## GraphPred vs every C+T threshold
graph_rows <- list()
graph_idx <- 0
for (d in c("BPD", "MDD", "SCZD")) {
  gp <- get_vector("GraphPred_bcftools", d, NA_real_)
  for (method in c("PRSice_C+T", "Shizhong_C+T")) {
    thrs <- sort(unique(long_scores$threshold[long_scores$method == method &
                                                long_scores$disorder == d]))
    for (thr in thrs) {
      ct <- get_vector(method, d, thr)
      graph_idx <- graph_idx + 1
      graph_rows[[graph_idx]] <- cbind(
        method_compared_to_graphpred = method,
        cor_row(gp, ct, paste0("GraphPred_", d),
                paste0(method, "_", d, "_Pt_", fmt_thr(thr)), d, thr),
        stringsAsFactors = FALSE
      )
    }
  }
}
graphpred_threshold_correlations <- do.call(rbind, graph_rows)
graphpred_threshold_correlations$best_abs_pearson_for_method <- FALSE
for (key in unique(paste(graphpred_threshold_correlations$disorder,
                         graphpred_threshold_correlations$method_compared_to_graphpred))) {
  idx <- which(paste(graphpred_threshold_correlations$disorder,
                     graphpred_threshold_correlations$method_compared_to_graphpred) == key)
  vals <- abs(graphpred_threshold_correlations$pearson[idx])
  if (all(is.na(vals))) next
  graphpred_threshold_correlations$best_abs_pearson_for_method[idx[which.max(vals)]] <- TRUE
}
write_tsv(graphpred_threshold_correlations,
          file.path(out_dir, "graphpred_threshold_correlations.tsv"))

## full within-disorder pairwise audit grid
pair_rows <- list()
pair_idx <- 0
for (d in c("BPD", "MDD", "SCZD")) {
  score_cols <- unique(long_scores$score_col[long_scores$disorder == d])
  for (i in seq_len(length(score_cols) - 1)) {
    for (j in seq((i + 1), length(score_cols))) {
      x <- long_scores[long_scores$score_col == score_cols[i], ]
      y <- long_scores[long_scores$score_col == score_cols[j], ]
      x <- x[match(reference_samples, x$SAMPLE), ]
      y <- y[match(reference_samples, y$SAMPLE), ]
      pair_idx <- pair_idx + 1
      pair_rows[[pair_idx]] <- data.frame(
        disorder = d,
        score_x = score_cols[i],
        method_x = unique(x$method),
        threshold_x = unique(x$threshold_label),
        score_y = score_cols[j],
        method_y = unique(y$method),
        threshold_y = unique(y$threshold_label),
        n = sum(is.finite(x$raw_score) & is.finite(y$raw_score)),
        pearson = safe_cor(x$raw_score, y$raw_score, "pearson"),
        spearman = safe_cor(x$raw_score, y$raw_score, "spearman"),
        z_rmse = rmse(x$z_score, y$z_score),
        stringsAsFactors = FALSE
      )
    }
  }
}
all_pairwise_correlations <- do.call(rbind, pair_rows)
write_tsv(all_pairwise_correlations,
          file.path(out_dir, "all_pairwise_correlations.tsv"))

## group summaries use standardized scores so methods are on the same scale
long_with_dx <- merge(long_scores, demog[, c("SAMPLE", "PrimaryDx")],
                      by = "SAMPLE", sort = FALSE)
summary_rows <- list()
summary_idx <- 0
group_keys <- unique(long_with_dx[, c("method", "disorder", "threshold_label",
                                      "score_col", "PrimaryDx")])
for (i in seq_len(nrow(group_keys))) {
  key <- group_keys[i, ]
  sub <- long_with_dx[
    long_with_dx$method == key$method &
      long_with_dx$disorder == key$disorder &
      long_with_dx$threshold_label == key$threshold_label &
      long_with_dx$score_col == key$score_col &
      long_with_dx$PrimaryDx == key$PrimaryDx,
  ]
  vals <- sub$z_score
  summary_idx <- summary_idx + 1
  summary_rows[[summary_idx]] <- data.frame(
    method = key$method,
    disorder = key$disorder,
    threshold_label = key$threshold_label,
    score_col = key$score_col,
    PrimaryDx = key$PrimaryDx,
    n = sum(is.finite(vals)),
    mean_z = mean(vals, na.rm = TRUE),
    sd_z = sd(vals, na.rm = TRUE),
    median_z = median(vals, na.rm = TRUE),
    iqr_z = IQR(vals, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
diagnosis_group_summary <- do.call(rbind, summary_rows)
write_tsv(diagnosis_group_summary,
          file.path(out_dir, "diagnosis_group_summary.tsv"))

readme_lines <- c(
  "# PRS Method Comparison Outputs",
  "",
  "This directory compares three PRS score sources for the same 119 MBv samples:",
  "",
  "- `Shizhong_C+T`: top-level `PRS_Bipolar.csv`, `PRS_MDD.csv`, and `PRS_SCZ.csv`.",
  "- `PRSice_C+T`: `../prsice/out/mbv_prsice_scores.tsv`.",
  "- `GraphPred_bcftools`: `../scores_bcftools/mbv_graphpred_scores.tsv`.",
  "",
  "## Threshold handling",
  "",
  "PRSice and Shizhong both output clumping-and-thresholding scores, but their threshold grids are not identical. Exact shared thresholds are compared one-to-one: `1e-06`, `1e-04`, `0.001`, `0.01`, `0.05`, `0.1`, `0.5`, and `1`. PRSice-only and Shizhong-only thresholds are listed in `threshold_availability.tsv`. GraphPred has one LD-aware score per disorder, so it is correlated against every PRSice and Shizhong threshold.",
  "",
  "## Tables",
  "",
  "- `score_inventory.tsv`: input source, method, disorder, score column list, threshold count, and sample count.",
  "- `threshold_availability.tsv`: threshold grid availability in PRSice and Shizhong by disorder.",
  "- `combined_scores_raw.tsv`: demographics plus all raw score columns.",
  "- `combined_scores_z.tsv`: demographics plus all score columns standardized across the 119 samples.",
  "- `matched_threshold_correlations.tsv`: PRSice vs Shizhong at exact shared thresholds.",
  "- `graphpred_threshold_correlations.tsv`: GraphPred vs each PRSice/Shizhong threshold, with best absolute Pearson flags.",
  "- `all_pairwise_correlations.tsv`: all within-disorder pairwise method/threshold correlations.",
  "- `diagnosis_group_summary.tsv`: z-score summaries by method, disorder, threshold, and PrimaryDx.",
  "",
  "## Notes",
  "",
  "Correlation signs and magnitudes reflect score conventions as stored in the source files. Z-score RMSE is computed after standardizing each score column over all 119 samples."
)
writeLines(readme_lines, file.path(out_dir, "README.md"))

message("Wrote PRS method comparison outputs to ", out_dir)
