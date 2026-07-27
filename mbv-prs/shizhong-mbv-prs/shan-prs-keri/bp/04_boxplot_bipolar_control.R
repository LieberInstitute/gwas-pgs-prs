 ## ---- setup ----
library(data.table)
library(stringr)
library(ggplot2)
library(ggforce)   # for facet_wrap_paginate()

# Paths
setwd("out")
prs_dir <- "."  # directory containing score.*.sscore files
phe_path <- "/dcs05/lieber/liebercentral/libdGenotype_LIBD001/BrainGenotyping/genotyped_brnum_pheno_3125.tab"

## ---- read all PRS score files ----
files <- list.files(prs_dir, pattern = "^score\\..+\\.sscore$", full.names = TRUE)
files <- files[!grepl("\\.log$", files, ignore.case = TRUE)]
stopifnot(length(files) > 0)

read_one <- function(f) {
  dt <- fread(f)
  setnames(dt, old = names(dt)[1], new = "IID")
  base <- basename(f)
  thr_label <- sub("^score\\.", "", base)
  thr_label <- sub("\\.sscore$", "", thr_label)
  thr_value <- suppressWarnings(as.numeric(thr_label))

 data.table(
    IID = dt$IID,
    score = dt$SCORE1_AVG,
    threshold_value = thr_value
  )
}

scores <- rbindlist(lapply(files, read_one), use.names = TRUE, fill = TRUE)

# Order thresholds by numeric value; fall back to lexicographic if NA
setorder(scores, threshold_value)
scores[, threshold_label := factor(
  paste0("p cutoff=", threshold_value),
  levels = paste0("p cutoff=", sort(unique(threshold_value)))
)]

## ---- WIDE PRS table (one column per threshold, ordered) ----
prs_wide <- dcast(
  scores,
  IID ~ threshold_label,
  value.var = "score",
  fun.aggregate = mean
)
setcolorder(prs_wide, c("IID", levels(scores$threshold_label)))
fwrite(prs_wide, file = "PRS_scores_by_threshold.csv")

## ---- read phenotype and harmonize ----
phe <- fread(phe_path, sep = "\t", header = TRUE)
setnames(phe, "BrNum", "IID")

# Keep only Control vs SCZD (map common schizophrenia labels to "SCZD")
phe[, Dx2 := fifelse(Dx %in% c("Control"), "Control",
              fifelse(Dx %in% c("Bipolar"), "Bipolar", NA_character_))]
phe <- phe[!is.na(Dx2)]
phe <- phe[, .(IID, Group = Dx2, Race)]

## ---- merge PRS with phenotype ----
dat <- merge(scores, phe, by = "IID")
dat <- dat[Race %in% c("CAUC")]
dat[, Group := factor(Group, levels = c("Control","Bipolar"))]

## ---- function to draw paginated boxplots with p-values ----
plot_prs_by_race <- function(df, race) {
  dfr <- df[Race == race]
  # Compute Wilcoxon p-values per threshold
  pvals <- dfr[, {
    if (length(unique(na.omit(Group))) < 2L) {
      .(p = NA_real_, ymax = max(score, na.rm=TRUE))
    } else {
      pv <- tryCatch(t.test(score ~ Group)$p.value, error = function(e) NA_real_)
      .(p = pv, ymax = max(score, na.rm=TRUE))
    }
  }, by = threshold_label]

  # Nicely formatted labels
  pvals[, p_label := ifelse(is.na(p), "p=NA", paste0("p=", formatC(p, format="e", digits=2)))]
  # Put text roughly above the boxes
  rng <- range(dfr$score, na.rm = TRUE)
  bump <- 0.06 * (rng[2] - rng[1] + 1e-12)
  pvals[, x := 1.5]                      # halfway between Control (1) and SCZD (2)
  pvals[, y := ymax + bump]

  # Number of pages (9 facets per page)
  n_thr <- nlevels(dfr$threshold_label)
  n_pages <- ceiling(n_thr / 9)

  # Create one PDF per race with multiple pages
  pdf(file = "bipolar_control_bipolarPRS_boxplot.pdf", width = 12, height = 12)
  for (pg in seq_len(n_pages)) {
    g <- ggplot(dfr, aes(x = Group, y = score, fill = Group)) +
      geom_boxplot(width = 0.7, outlier.shape = NA) +
	  geom_jitter(
		color = "black",
		width = 0.18,
		size = 1.4,
		alpha = 0.6
	  ) +
      geom_text(data = pvals, aes(x = x, y = y, label = p_label),
                inherit.aes = FALSE, size = 4) +
      ggforce::facet_wrap_paginate(~ threshold_label, ncol = 3, nrow = 3, page = pg, scales = "free_y") +
      labs(title = "Bipolar PRS",
           x = NULL, y = "PRS") +
      scale_fill_manual(values = c("Control" = "#4C78A8", "Bipolar" = "#F58518")) +
      theme_bw(base_size = 16) +
      theme(legend.position = "top",
            strip.background = element_rect(fill = "grey90"),
            axis.text.x = element_text(angle = 0, hjust = 0.5))
    print(g)
  }
  dev.off()
}

## ---- run plots ----
plot_prs_by_race(dat, "CAUC")

