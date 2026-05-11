options(stringsAsFactors = FALSE)

outdir <- "prs_diagnosis_assoc"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

demog_path <- "MBv_demographics_n119.tab"
scores_path <- "mbv_prs_merged_scores.tsv"
pca_path <- "qc/mbv.pca.eigenvec"
overlap_counts_path <- "gwas_sig_overlaps/gwas_significant_overlap_counts.csv"
overlap_pairwise_path <- "gwas_sig_overlaps/gwas_significant_pairwise_overlaps.csv"

## read local phenotype, score, and ancestry-PC inputs
demog <- read.delim(demog_path, check.names = FALSE, quote = "", comment.char = "")
demog <- demog[names(demog) != ""]
scores <- read.delim(scores_path, check.names = FALSE)
pca <- read.delim(pca_path, check.names = FALSE)

fail <- function(msg) stop(msg, call. = FALSE)

if (anyDuplicated(demog$BrNum)) fail("duplicate BrNum values in demographics")
if (anyDuplicated(scores$SAMPLE)) fail("duplicate SAMPLE values in PRS score table")
if (anyDuplicated(pca$IID)) fail("duplicate IID values in PCA table")
if (!all(demog$BrNum %in% scores$SAMPLE)) fail("some demographic samples are missing from PRS table")
if (!all(demog$BrNum %in% pca$IID)) fail("some demographic samples are missing from PCA table")

## join by sample IDs used by the previous PRS workflow
x <- merge(demog, scores, by.x = "BrNum", by.y = "SAMPLE", all = FALSE)
pca$BrNum <- pca$IID
x <- merge(x, pca, by = "BrNum", all = FALSE)

if (nrow(x) != 119) fail("joined table does not contain 119 samples")
dx_counts <- table(x$PrimaryDx)
expected_dx <- c(Bipolar = 40, Control = 40, MDD = 39)
if (!all(expected_dx[names(expected_dx)] == dx_counts[names(expected_dx)])) {
  fail("diagnosis counts do not match expected 40 Bipolar, 40 Control, 39 MDD")
}

score_cols <- setdiff(names(scores), "SAMPLE")
pc_cols <- paste0("PC", 1:5)
if (!all(pc_cols %in% names(x))) fail("PC1-PC5 are not available")

## define score metadata from stable column-name patterns
score_meta <- data.frame(score = score_cols)
score_meta$method <- ifelse(grepl("^PRSice_", score_meta$score), "PRSice", "GraphPred_score")
score_meta$score_disorder <- sub("^PRSice_([^_]+)_.*$", "\\1", score_meta$score)
score_meta$score_disorder[score_meta$method == "GraphPred_score"] <-
  sub("^([^_]+)_.*$", "\\1", score_meta$score[score_meta$method == "GraphPred_score"])
score_meta$threshold <- ifelse(
  score_meta$method == "PRSice",
  sub("^PRSice_[^_]+_Pt_", "", score_meta$score),
  "GraphPred"
)
score_meta$threshold_numeric <- suppressWarnings(as.numeric(score_meta$threshold))

if (any(!score_meta$score_disorder %in% c("BPD", "MDD", "SCZD"))) {
  fail("unexpected score disorder in PRS columns")
}

for (col in c(score_cols, pc_cols)) {
  if (any(!is.finite(x[[col]]))) fail(paste("non-finite values found in", col))
}

contrasts <- list(
  MDD_vs_Control = list(case = "MDD", control = "Control", case_values = "MDD"),
  BPD_vs_Control = list(case = "Bipolar", control = "Control", case_values = "Bipolar"),
  CaseUnion_vs_Control = list(case = "MDD+Bipolar", control = "Control", case_values = c("MDD", "Bipolar")),
  BPD_vs_MDD = list(case = "Bipolar", control = "MDD", case_values = "Bipolar", control_values = "MDD")
)

auc_rank <- function(y, s) {
  ## rank-based AUC equals probability a random case scores above a random control
  r <- rank(s, ties.method = "average")
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cohen_d <- function(case_vals, control_vals) {
  n1 <- length(case_vals)
  n0 <- length(control_vals)
  sp <- sqrt(((n1 - 1) * var(case_vals) + (n0 - 1) * var(control_vals)) / (n1 + n0 - 2))
  (mean(case_vals) - mean(control_vals)) / sp
}

fit_one <- function(dat, y, score_z, adjusted) {
  ## fit logistic model and return PRS coefficient statistics
  df <- data.frame(y = y, PRS_z = score_z, Sex = factor(dat$Sex), dat[, pc_cols, drop = FALSE])
  form <- if (adjusted) {
    as.formula(paste("y ~ PRS_z + Sex +", paste(pc_cols, collapse = " + ")))
  } else {
    y ~ PRS_z
  }
  fit <- glm(form, family = binomial(), data = df)
  s <- coef(summary(fit))
  if (!"PRS_z" %in% rownames(s)) fail("PRS_z coefficient missing from model")
  beta <- unname(coef(fit)["PRS_z"])
  se <- unname(s["PRS_z", "Std. Error"])
  data.frame(
    beta = beta,
    se = se,
    z = unname(s["PRS_z", "z value"]),
    p = unname(s["PRS_z", "Pr(>|z|)"]),
    OR = exp(beta),
    OR_low = exp(beta - 1.96 * se),
    OR_high = exp(beta + 1.96 * se)
  )
}

assoc_rows <- list()
group_rows <- list()
k <- 1
gk <- 1

for (contrast_name in names(contrasts)) {
  cc <- contrasts[[contrast_name]]
  control_values <- if (!is.null(cc$control_values)) cc$control_values else cc$control
  keep_dx <- x$PrimaryDx %in% c(cc$case_values, control_values)
  dat <- x[keep_dx, ]
  y <- ifelse(dat$PrimaryDx %in% cc$case_values, 1, 0)

  for (score in score_cols) {
    score_z <- as.numeric(scale(dat[[score]], center = mean(x[[score]]), scale = sd(x[[score]])))
    case_vals <- dat[[score]][y == 1]
    control_vals <- dat[[score]][y == 0]

    meta <- score_meta[score_meta$score == score, ]
    group_rows[[gk]] <- data.frame(
      contrast = contrast_name,
      score = score,
      method = meta$method,
      score_disorder = meta$score_disorder,
      threshold = meta$threshold,
      n_case = sum(y == 1),
      n_control = sum(y == 0),
      mean_case = mean(case_vals),
      mean_control = mean(control_vals),
      mean_diff = mean(case_vals) - mean(control_vals),
      cohen_d = cohen_d(case_vals, control_vals)
    )
    gk <- gk + 1

    for (adjusted in c(FALSE, TRUE)) {
      fit_stats <- fit_one(dat, y, score_z, adjusted)
      assoc_rows[[k]] <- cbind(
        data.frame(
          contrast = contrast_name,
          case_label = cc$case,
          control_label = cc$control,
          model = ifelse(adjusted, "adjusted_sex_PC1_PC5", "unadjusted"),
          score = score,
          method = meta$method,
          score_disorder = meta$score_disorder,
          threshold = meta$threshold,
          threshold_numeric = meta$threshold_numeric,
          n = length(y),
          n_case = sum(y == 1),
          n_control = sum(y == 0),
          AUC = auc_rank(y, score_z),
          mean_case = mean(case_vals),
          mean_control = mean(control_vals),
          mean_diff = mean(case_vals) - mean(control_vals),
          cohen_d = cohen_d(case_vals, control_vals)
        ),
        fit_stats
      )
      k <- k + 1
    }
  }
}

assoc <- do.call(rbind, assoc_rows)
group_summary <- do.call(rbind, group_rows)
assoc$q_fdr_all <- p.adjust(assoc$p, method = "BH")
assoc$q_fdr_within_contrast_model <- ave(
  assoc$p,
  assoc$contrast,
  assoc$model,
  FUN = function(z) p.adjust(z, method = "BH")
)

## flag internally best PRSice thresholds by contrast, model, and score family
assoc$exploratory_best_internal <- FALSE
prsice_idx <- which(assoc$method == "PRSice")
best_rows <- c()
for (cn in unique(assoc$contrast)) {
  for (md in unique(assoc$model)) {
    for (sd in unique(assoc$score_disorder)) {
      idx <- prsice_idx[assoc$contrast[prsice_idx] == cn &
                          assoc$model[prsice_idx] == md &
                          assoc$score_disorder[prsice_idx] == sd]
      if (length(idx) > 0) {
        best_rows <- c(best_rows, idx[which.max(assoc$AUC[idx])])
      }
    }
  }
}
assoc$exploratory_best_internal[best_rows] <- TRUE
best <- assoc[assoc$exploratory_best_internal, ]
best <- best[order(best$contrast, best$model, best$score_disorder), ]

write.csv(assoc, file.path(outdir, "prs_association_all_scores.csv"), row.names = FALSE, quote = TRUE)
write.csv(best, file.path(outdir, "prs_best_thresholds_exploratory.csv"), row.names = FALSE, quote = TRUE)
write.csv(group_summary, file.path(outdir, "prs_group_summary.csv"), row.names = FALSE, quote = TRUE)

## prepare overlap context from the previous significant-variant analysis
counts <- read.csv(overlap_counts_path)
pairwise <- read.csv(overlap_pairwise_path)
exact_counts <- counts[counts$universe == "exact_variant_CHROM_POS_REF_ALT", ]
get_count <- function(region) exact_counts$count[match(region, exact_counts$region)]
mdd_bpd <- pairwise[pairwise$universe == "exact_variant_CHROM_POS_REF_ALT" &
                      pairwise$pair == "BPD_MDD", ]
overlap_note <- data.frame(
  metric = c(
    "BPD_significant_variants",
    "MDD_significant_variants",
    "BPD_MDD_shared_significant_variants",
    "BPD_MDD_jaccard",
    "fraction_BPD_significant_shared_with_MDD",
    "fraction_MDD_significant_shared_with_BPD"
  ),
  value = c(
    get_count("BPD_total"),
    get_count("MDD_total"),
    mdd_bpd$intersection,
    mdd_bpd$jaccard,
    mdd_bpd$intersection / get_count("BPD_total"),
    mdd_bpd$intersection / get_count("MDD_total")
  )
)
write.csv(overlap_note, file.path(outdir, "prs_overlap_interpretation.csv"), row.names = FALSE, quote = TRUE)

if (!requireNamespace("ggplot2", quietly = TRUE)) fail("ggplot2 is required for plots")
library(ggplot2)

plot_dat <- assoc[assoc$model == "unadjusted", ]
plot_dat$score_disorder <- factor(plot_dat$score_disorder, levels = c("BPD", "MDD", "SCZD"))
p_auc <- ggplot(plot_dat, aes(x = threshold_numeric, y = AUC, color = score_disorder, group = score_disorder)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey55") +
  geom_line(data = plot_dat[plot_dat$method == "PRSice", ], linewidth = 0.5) +
  geom_point(data = plot_dat[plot_dat$method == "PRSice", ], size = 1.8) +
  geom_point(data = plot_dat[plot_dat$method == "GraphPred_score", ], aes(x = 1.2), shape = 17, size = 2.5) +
  scale_x_log10(na.value = 1.2) +
  facet_wrap(~ contrast, ncol = 2) +
  labs(x = "PRSice p-value threshold; triangle at right = GraphPred", y = "AUC", color = "Score") +
  theme_bw(base_size = 10)
ggsave(file.path(outdir, "prs_auc_by_threshold.svg"), p_auc, width = 9, height = 6)

top_scores <- unique(best$score[best$model == "unadjusted"])
dist_rows <- list()
dk <- 1
for (contrast_name in names(contrasts)) {
  cc <- contrasts[[contrast_name]]
  control_values <- if (!is.null(cc$control_values)) cc$control_values else cc$control
  keep_dx <- x$PrimaryDx %in% c(cc$case_values, control_values)
  dat <- x[keep_dx, ]
  group <- ifelse(dat$PrimaryDx %in% cc$case_values, cc$case, cc$control)
  for (score in top_scores) {
    meta <- score_meta[score_meta$score == score, ]
    dist_rows[[dk]] <- data.frame(
      contrast = contrast_name,
      group = group,
      score = score,
      score_label = paste(meta$score_disorder, meta$threshold, sep = " Pt="),
      PRS_z = as.numeric(scale(dat[[score]], center = mean(x[[score]]), scale = sd(x[[score]])))
    )
    dk <- dk + 1
  }
}
dist_dat <- do.call(rbind, dist_rows)
p_dist <- ggplot(dist_dat, aes(x = group, y = PRS_z, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, height = 0, size = 0.9, alpha = 0.7) +
  facet_grid(score_label ~ contrast, scales = "free_x") +
  labs(x = NULL, y = "PRS z-score") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(outdir, "prs_selected_distributions.svg"), p_dist, width = 11, height = 9)

top_unadj <- assoc[assoc$model == "unadjusted", ]
top_unadj <- top_unadj[order(top_unadj$contrast, -top_unadj$AUC), ]
top_unadj <- do.call(rbind, lapply(split(top_unadj, top_unadj$contrast), head, 5))

md <- file.path(outdir, "prs_diagnosis_association_report.md")
sink(md)
cat("# PRS Diagnosis Association Report\n\n")
cat("Samples: 119 total; 39 MDD, 40 Bipolar/BPD, 40 Controls. No SCZ diagnosis group is present, so SCZD PRS is evaluated only as a cross-disorder score.\n\n")
cat("Models: unadjusted logistic regression and sensitivity logistic regression adjusted for Sex plus PC1-PC5. Each PRS was standardized over all 119 samples. PRSice best thresholds are internal exploratory flags, not confirmatory choices.\n\n")
cat("## Best Internal PRSice Thresholds\n\n")
print(best[, c("contrast", "model", "score_disorder", "threshold", "n_case", "n_control", "AUC", "OR", "OR_low", "OR_high", "p", "q_fdr_within_contrast_model", "cohen_d")], row.names = FALSE)
cat("\n## Top Unadjusted Scores by AUC\n\n")
print(top_unadj[, c("contrast", "score", "method", "score_disorder", "threshold", "AUC", "OR", "p", "cohen_d")], row.names = FALSE)
cat("\n## GWAS Significant Variant Overlap Context\n\n")
print(overlap_note, row.names = FALSE)
cat("\nMDD-BPD genome-wide significant exact-variant overlap is modest by Jaccard, but asymmetric by disorder because BPD has fewer significant variants. This overlap is not equivalent to PRS overlap because PRSice thresholds include sub-significant variants and GraphPred uses LD-aware genome-wide weights.\n\n")
cat("## Threshold Selection Guidance\n\n")
cat("Do not select a PRSice threshold for confirmatory downstream analysis solely because it separates MBv cases and controls best. This is overfit at n=119. Prefer a pre-specified score, an external validation cohort, cross-validation if sample size allows, or PRSice permutation/empirical p-values for association. Report all thresholds here and label internally best thresholds as exploratory.\n\n")
cat("## Output Files\n\n")
cat("- `prs_association_all_scores.csv`\n")
cat("- `prs_best_thresholds_exploratory.csv`\n")
cat("- `prs_group_summary.csv`\n")
cat("- `prs_overlap_interpretation.csv`\n")
cat("- `prs_auc_by_threshold.svg`\n")
cat("- `prs_selected_distributions.svg`\n")
sink()

cat("wrote outputs to", outdir, "\n")
