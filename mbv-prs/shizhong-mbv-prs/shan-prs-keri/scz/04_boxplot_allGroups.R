## ---- setup ----
library(data.table)
library(stringr)
library(ggplot2)
library(ggforce)   # facet_wrap_paginate()

setwd("out")
prs_dir  <- "."
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

setorder(scores, threshold_value)
scores[, threshold_label := factor(
  paste0("p cutoff=", threshold_value),
  levels = paste0("p cutoff=", sort(unique(threshold_value)))
)]

## ---- WIDE PRS table ----
prs_wide <- dcast(scores, IID ~ threshold_label, value.var = "score", fun.aggregate = mean)
setcolorder(prs_wide, c("IID", levels(scores$threshold_label)))
fwrite(prs_wide, file = "PRS_scores_by_threshold.csv")

## ---- read phenotype and harmonize ----
phe <- fread(phe_path, sep = "\t", header = TRUE)
setnames(phe, "BrNum", "IID")

# Keep Control, MDD, Bipolar
phe[, Dx2 := fifelse(Dx %in% c("Control"), "Control",
              fifelse(Dx %in% c("MDD"), "MDD",
              fifelse(Dx %in% c("Bipolar"), "Bipolar", NA_character_)))]
phe <- phe[!is.na(Dx2)]
phe <- phe[, .(IID, Group = Dx2, Race)]

## ---- merge PRS with phenotype ----
dat <- merge(scores, phe, by = "IID")
dat <- dat[Race %in% c("CAUC")]
dat[, Group := factor(Group, levels = c("Control", "MDD", "Bipolar"))]

## ---- function to draw paginated boxplots with p-values ----
plot_prs_by_race <- function(df, race) {
  dfr <- df[Race == race]

  # Compute p-values + per-threshold y-range
  pvals <- dfr[, {
    sc <- score[is.finite(score)]
    ymin <- if (length(sc)) min(sc) else NA_real_
    ymax <- if (length(sc)) max(sc) else NA_real_
    bump <- 0.06 * ((ymax - ymin) + 1e-12)

    pv_mdd <- tryCatch({
      tmp <- .SD[Group %in% c("Control", "MDD")]
      if (length(unique(tmp$Group)) < 2L) NA_real_
      else t.test(score ~ Group, data = tmp)$p.value
    }, error = function(e) NA_real_)

    pv_bip <- tryCatch({
      tmp <- .SD[Group %in% c("Control", "Bipolar")]
      if (length(unique(tmp$Group)) < 2L) NA_real_
      else t.test(score ~ Group, data = tmp)$p.value
    }, error = function(e) NA_real_)

    .(
      p_mdd = pv_mdd,
      p_bip = pv_bip,
      ymin = ymin,
      ymax = ymax,
      bump = bump
    )
  }, by = threshold_label]

  # Format labels
  pvals[, `:=`(
    lab_mdd = ifelse(is.na(p_mdd), "p=NA",
                     paste0("p=", formatC(p_mdd, format = "e", digits = 2))),
    lab_bip = ifelse(is.na(p_bip), "p=NA",
                     paste0("p=", formatC(p_bip, format = "e", digits = 2)))
  )]

  # Bracket geometry (x positions correspond to factor levels)
  # Control=1, MDD=2, Bipolar=3
  pvals[, `:=`(
  x1_mdd = 1, x2_mdd = 2,
  x1_bip = 1, x2_bip = 3,

  # bracket heights
  y_mdd = ymax + 1.2 * bump,
  y_bip = ymax + 4 * bump,

  # text heights (above brackets)
  ytxt_mdd = ymax + 1.8 * bump,
  ytxt_bip = ymax + 4.6 * bump
  )]

  # Pages
  n_thr <- nlevels(dfr$threshold_label)
  n_pages <- ceiling(n_thr / 9)

  pdf(file = "PRS_SCZ.pdf", width = 12, height = 12)
  for (pg in seq_len(n_pages)) {

    g <- ggplot(dfr, aes(x = Group, y = score, fill = Group)) +
      geom_boxplot(width = 0.7, outlier.shape = NA) +
      geom_jitter(color = "black", width = 0.18, size = 1.4, alpha = 0.6) +

      ## ---- MDD vs Control bracket ----
      geom_segment(
        data = pvals,
        aes(x = x1_mdd, xend = x2_mdd, y = y_mdd, yend = y_mdd),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = pvals,
        aes(x = (x1_mdd + x2_mdd) / 2, y = ytxt_mdd, label = lab_mdd),
        inherit.aes = FALSE,
        size = 4
      ) +

      ## ---- Bipolar vs Control bracket ----
      geom_segment(
        data = pvals,
        aes(x = x1_bip, xend = x2_bip, y = y_bip, yend = y_bip),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = pvals,
        aes(x = (x1_bip + x2_bip) / 2, y = ytxt_bip, label = lab_bip),
        inherit.aes = FALSE,
        size = 4
      ) +

      ggforce::facet_wrap_paginate(
        ~ threshold_label,
        ncol = 3, nrow = 3, page = pg,
        scales = "free_y"
      ) +
      labs(title = "", x = NULL, y = "PRS") +
      scale_fill_manual(
        values = c("Control" = "#4C78A8", "MDD" = "#F58518", "Bipolar" = "#54A24B")
      ) +
      theme_bw(base_size = 16) +
      theme(
        legend.position = "top",
        strip.background = element_rect(fill = "grey90"),
        axis.text.x = element_text(angle = 0, hjust = 0.5)
      )

    print(g)
  }
  dev.off()
}

## ---- run plots ----
plot_prs_by_race(dat, "CAUC")
