options(stringsAsFactors = FALSE)

bcf <- "/opt/sw/bin/bcftools"
threshold_p <- 5e-8
threshold_lp <- -log10(threshold_p)

files <- c(
  BPD = "/home/gpertea/work/ref/GWAS/BPD/bip2024_eur_no23andMe.hg38.bcf",
  MDD = "/home/gpertea/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf",
  SCZD = "/home/gpertea/work/ref/GWAS/SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf"
)

outdir <- "gwas_sig_overlaps"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_sig <- function(label, path) {
  ## query GWAS-VCF summary fields and filter in R for transparent thresholding
  cmd <- sprintf(
    "%s query -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT[\\t%%LP]\\n' %s",
    shQuote(bcf), shQuote(path)
  )
  x <- read.table(pipe(cmd), sep = "\t", header = FALSE, quote = "",
                  comment.char = "", fill = TRUE,
                  col.names = c("CHROM", "POS", "ID", "REF", "ALT", "LP"))
  x$LP <- as.numeric(x$LP)
  x <- x[is.finite(x$LP) & x$LP >= threshold_lp, ]
  x$disorder <- label
  x$variant_id <- paste(x$CHROM, x$POS, x$REF, x$ALT, sep = ":")
  x$position_id <- paste(x$CHROM, x$POS, sep = ":")
  x$P <- 10^(-x$LP)
  x
}

sig <- do.call(rbind, Map(read_sig, names(files), files))
write.csv(sig, file.path(outdir, "gwas_significant_variants_long.csv"),
          row.names = FALSE, quote = TRUE)

make_membership <- function(sig, id_col, prefix) {
  ids <- sort(unique(sig[[id_col]]))
  out <- data.frame(id = ids)
  for (d in names(files)) {
    y <- sig[sig$disorder == d, ]
    hit <- match(ids, y[[id_col]])
    out[[d]] <- !is.na(hit)
    out[[paste0(d, "_LP")]] <- NA_real_
    out[[paste0(d, "_P")]] <- NA_real_
    out[[paste0(d, "_ID")]] <- NA_character_
    out[[paste0(d, "_LP")]][!is.na(hit)] <- y$LP[hit[!is.na(hit)]]
    out[[paste0(d, "_P")]][!is.na(hit)] <- y$P[hit[!is.na(hit)]]
    out[[paste0(d, "_ID")]][!is.na(hit)] <- y$ID[hit[!is.na(hit)]]
  }
  out$n_disorders <- rowSums(out[names(files)])
  out$membership <- apply(out[names(files)], 1, function(z) paste(names(files)[z], collapse = "&"))
  write.csv(out, file.path(outdir, paste0(prefix, "_membership.csv")),
            row.names = FALSE, quote = TRUE)
  out
}

membership_exact <- make_membership(sig, "variant_id", "exact_variant")
membership_pos <- make_membership(sig, "position_id", "position")

region_counts <- function(m, universe_label) {
  b <- m$BPD
  d <- m$MDD
  s <- m$SCZD
  data.frame(
    universe = universe_label,
    region = c("BPD_only", "MDD_only", "SCZD_only",
               "BPD_MDD_only", "BPD_SCZD_only", "MDD_SCZD_only",
               "BPD_MDD_SCZD", "BPD_total", "MDD_total", "SCZD_total",
               "union_total"),
    count = c(
      sum(b & !d & !s),
      sum(!b & d & !s),
      sum(!b & !d & s),
      sum(b & d & !s),
      sum(b & !d & s),
      sum(!b & d & s),
      sum(b & d & s),
      sum(b),
      sum(d),
      sum(s),
      nrow(m)
    )
  )
}

counts <- rbind(
  region_counts(membership_exact, "exact_variant_CHROM_POS_REF_ALT"),
  region_counts(membership_pos, "position_CHROM_POS")
)
write.csv(counts, file.path(outdir, "gwas_significant_overlap_counts.csv"),
          row.names = FALSE, quote = TRUE)

pairwise <- function(m, universe_label) {
  data.frame(
    universe = universe_label,
    pair = c("BPD_MDD", "BPD_SCZD", "MDD_SCZD"),
    intersection = c(sum(m$BPD & m$MDD), sum(m$BPD & m$SCZD), sum(m$MDD & m$SCZD)),
    union = c(sum(m$BPD | m$MDD), sum(m$BPD | m$SCZD), sum(m$MDD | m$SCZD)),
    jaccard = c(
      sum(m$BPD & m$MDD) / sum(m$BPD | m$MDD),
      sum(m$BPD & m$SCZD) / sum(m$BPD | m$SCZD),
      sum(m$MDD & m$SCZD) / sum(m$MDD | m$SCZD)
    )
  )
}

pairwise_counts <- rbind(
  pairwise(membership_exact, "exact_variant_CHROM_POS_REF_ALT"),
  pairwise(membership_pos, "position_CHROM_POS")
)
write.csv(pairwise_counts, file.path(outdir, "gwas_significant_pairwise_overlaps.csv"),
          row.names = FALSE, quote = TRUE)

draw_venn <- function(counts_one, title, path) {
  ## fixed-circle schematic, not area-proportional
  vals <- setNames(counts_one$count, counts_one$region)
  svg(path, width = 7, height = 5)
  par(mar = c(0, 0, 2, 0))
  plot.new()
  plot.window(xlim = c(0, 10), ylim = c(0, 7))
  symbols(4, 4.1, circles = 2.25, inches = FALSE, add = TRUE,
          bg = adjustcolor("#D95F02", 0.30), fg = "#D95F02", lwd = 2)
  symbols(6, 4.1, circles = 2.25, inches = FALSE, add = TRUE,
          bg = adjustcolor("#1B9E77", 0.30), fg = "#1B9E77", lwd = 2)
  symbols(5, 2.45, circles = 2.25, inches = FALSE, add = TRUE,
          bg = adjustcolor("#7570B3", 0.30), fg = "#7570B3", lwd = 2)
  text(2.6, 6.35, sprintf("BPD\nn=%s", vals["BPD_total"]), cex = 1.0)
  text(7.4, 6.35, sprintf("MDD\nn=%s", vals["MDD_total"]), cex = 1.0)
  text(5.0, 0.25, sprintf("SCZD\nn=%s", vals["SCZD_total"]), cex = 1.0)
  text(3.1, 4.35, vals["BPD_only"], cex = 1.1)
  text(6.9, 4.35, vals["MDD_only"], cex = 1.1)
  text(5.0, 1.45, vals["SCZD_only"], cex = 1.1)
  text(5.0, 4.85, vals["BPD_MDD_only"], cex = 1.1)
  text(4.05, 3.0, vals["BPD_SCZD_only"], cex = 1.1)
  text(5.95, 3.0, vals["MDD_SCZD_only"], cex = 1.1)
  text(5.0, 3.55, vals["BPD_MDD_SCZD"], cex = 1.1, font = 2)
  title(main = title)
  text(8.5, 0.7, sprintf("Union: %s", vals["union_total"]), cex = 0.9)
  dev.off()
}

draw_venn(
  counts[counts$universe == "exact_variant_CHROM_POS_REF_ALT", ],
  sprintf("Genome-wide significant overlap, exact variant, P <= %.0e", threshold_p),
  file.path(outdir, "gwas_significant_exact_variant_venn.svg")
)
draw_venn(
  counts[counts$universe == "position_CHROM_POS", ],
  sprintf("Genome-wide significant overlap, position, P <= %.0e", threshold_p),
  file.path(outdir, "gwas_significant_position_venn.svg")
)

md <- file.path(outdir, "gwas_significant_overlap_report.md")
sink(md)
cat("# GWAS Significant Variant Overlap Report\n\n")
cat(sprintf("Threshold: `P <= %.0e`, equivalent to `LP >= %.5f`.\n\n", threshold_p, threshold_lp))
cat("Primary overlap unit: exact `CHROM:POS:REF:ALT` variant ID.\n")
cat("Secondary overlap unit: `CHROM:POS` position ID.\n\n")
cat("## Exact Variant Counts\n\n")
print(counts[counts$universe == "exact_variant_CHROM_POS_REF_ALT", c("region", "count")], row.names = FALSE)
cat("\n## Position Counts\n\n")
print(counts[counts$universe == "position_CHROM_POS", c("region", "count")], row.names = FALSE)
cat("\n## Pairwise Overlaps\n\n")
print(pairwise_counts, row.names = FALSE)
cat("\n## Files\n\n")
cat("- `gwas_significant_variants_long.csv`\n")
cat("- `exact_variant_membership.csv`\n")
cat("- `position_membership.csv`\n")
cat("- `gwas_significant_overlap_counts.csv`\n")
cat("- `gwas_significant_pairwise_overlaps.csv`\n")
cat("- `gwas_significant_exact_variant_venn.svg`\n")
cat("- `gwas_significant_position_venn.svg`\n")
sink()
