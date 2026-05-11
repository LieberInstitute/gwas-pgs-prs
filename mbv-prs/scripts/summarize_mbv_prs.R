options(stringsAsFactors = FALSE)

zscore <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

read_table <- function(path) {
  read.table(path, header = TRUE, sep = "", check.names = FALSE, quote = "",
             comment.char = "", fill = TRUE)
}

read_graphpred <- function(label) {
  path <- file.path("scores_bcftools", paste0(label, ".graphpred.tsv"))
  x <- read_table(path)
  if (!"SAMPLE" %in% names(x)) names(x)[1] <- "SAMPLE"
  names(x)[names(x) != "SAMPLE"] <- paste(label, names(x)[names(x) != "SAMPLE"], sep = "_")
  x
}

read_prsice <- function(label) {
  path <- file.path("prsice", "out", paste0(label, ".all_score"))
  if (!file.exists(path)) return(NULL)
  x <- read_table(path)
  sample_col <- if ("IID" %in% names(x)) "IID" else names(x)[1]
  x$SAMPLE <- x[[sample_col]]
  x$SAMPLE <- sub("^([^_]+)_\\1$", "\\1", x$SAMPLE)
  keep <- setdiff(names(x), c("FID", "IID", "SAMPLE", "In_Regression"))
  x <- x[c("SAMPLE", keep)]
  names(x)[names(x) != "SAMPLE"] <- paste("PRSice", label, names(x)[names(x) != "SAMPLE"], sep = "_")
  x
}

merge_all <- function(xs) {
  xs <- Filter(Negate(is.null), xs)
  Reduce(function(a, b) merge(a, b, by = "SAMPLE", all = TRUE), xs)
}

write_z <- function(raw_path, z_path) {
  x <- read.table(raw_path, header = TRUE, sep = "\t", check.names = FALSE)
  z <- x
  for (nm in setdiff(names(z), "SAMPLE")) {
    if (is.numeric(z[[nm]]) || suppressWarnings(all(is.na(as.numeric(z[[nm]])) == is.na(z[[nm]])))) {
      z[[nm]] <- zscore(z[[nm]])
    }
  }
  write.table(z, z_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

graph <- merge_all(lapply(c("BPD", "MDD", "SCZD"), read_graphpred))
write.table(graph, "scores_bcftools/mbv_graphpred_scores.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write_z("scores_bcftools/mbv_graphpred_scores.tsv",
        "scores_bcftools/mbv_graphpred_scores_z.tsv")

prsice <- merge_all(lapply(c("BPD", "MDD", "SCZD"), read_prsice))
if (!is.null(prsice)) {
  write.table(prsice, "prsice/out/mbv_prsice_scores.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write_z("prsice/out/mbv_prsice_scores.tsv", "prsice/out/mbv_prsice_scores_z.tsv")
  merged <- merge(graph, prsice, by = "SAMPLE", all = TRUE)
} else {
  merged <- graph
}
write.table(merged, "mbv_prs_merged_scores.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

count_lines <- function(path, header = TRUE) {
  if (!file.exists(path)) return(NA_integer_)
  n <- length(readLines(path, warn = FALSE))
  if (header) n <- max(0L, n - 1L)
  n
}

variant_counts <- if (file.exists("qc/variant_counts.tsv")) {
  read.table("qc/variant_counts.tsv", header = TRUE, sep = "\t", check.names = FALSE)
} else {
  data.frame(stage = character(), variants = integer())
}
get_variant_count <- function(stage) {
  val <- variant_counts$variants[variant_counts$stage == stage]
  if (length(val) == 0) return(NA_integer_)
  val[1]
}

report <- data.frame(
  metric = c("input_samples", "final_samples", "input_variants", "final_variants",
             "ds_out_of_range_sites", "related_pairs_ge_0.0884", "pca_samples"),
  value = c(
    count_lines("qc/vcf.samples", header = FALSE),
    count_lines("qc/mbv.qc.psam", header = TRUE),
    get_variant_count("input"),
    get_variant_count("final_qc"),
    count_lines("qc/ds_out_of_range_sites.tsv", header = FALSE),
    count_lines("qc/mbv.king.kin0", header = TRUE),
    count_lines("qc/mbv.pca.eigenvec", header = TRUE)
  )
)

if (file.exists("qc/scoring_overlap.tsv")) {
  overlap <- read.table("qc/scoring_overlap.tsv", header = TRUE, sep = "\t",
                        check.names = FALSE)
  report <- rbind(report, data.frame(metric = paste0("overlap_", overlap$disorder),
                                     value = overlap$variants))
}

write.table(report, "qc/target_qc_report.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
