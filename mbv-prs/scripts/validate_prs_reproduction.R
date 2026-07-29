#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name) {
  idx <- match(name, args)
  if (is.na(idx) || idx == length(args)) stop("Missing argument: ", name, call. = FALSE)
  args[[idx + 1L]]
}

registry_file <- get_arg("--registry")
archive_root <- get_arg("--archive-root")
output <- get_arg("--output")
d <- read.delim(registry_file, check.names = FALSE, stringsAsFactors = FALSE)

results <- list()
for (trait in c("BD", "MDD")) {
  archive_name <- if (trait == "BD") "PRS_Bipolar.csv" else "PRS_MDD.csv"
  old <- read.csv(file.path(archive_root, "shizhong-mbv-prs", archive_name),
                  check.names = FALSE, stringsAsFactors = FALSE)
  new <- d[d$TRAIT == trait & d$GWAS_VERSION == "no23andMe" &
             d$METHOD == "Shizhong_CT" & d$SETUP == "historical", ]
  for (threshold in unique(new$THRESHOLD)) {
    column <- paste0("p cutoff=", threshold)
    if (!column %in% names(old)) stop("Archived threshold missing: ", column, call. = FALSE)
    score <- new[new$THRESHOLD == threshold, c("SAMPLE", "RAW_SCORE")]
    merged <- merge(old[, c("IID", column)], score, by.x = "IID", by.y = "SAMPLE")
    correlation <- cor(merged[[column]], merged$RAW_SCORE)
    results[[length(results) + 1L]] <- data.frame(
      TRAIT = trait, METHOD = "Shizhong_CT", THRESHOLD = threshold,
      METRIC = "pearson", VALUE = correlation, REQUIRED = 0.999,
      PASS = correlation >= 0.999, stringsAsFactors = FALSE
    )
  }
}

for (trait in c("BD", "MDD")) {
  archive_label <- if (trait == "BD") "BPD" else "MDD"
  old_path <- file.path(archive_root, "scores_bcftools",
                        paste0(archive_label, ".graphpred.tsv"))
  old <- read.delim(old_path, check.names = FALSE, stringsAsFactors = FALSE)
  new <- d[d$TRAIT == trait & d$GWAS_VERSION == "no23andMe" &
             d$METHOD == "GraphPred", c("SAMPLE", "RAW_SCORE")]
  merged <- merge(old[, 1:2], new, by.x = 1L, by.y = "SAMPLE")
  max_error <- max(abs(merged[[2L]] - merged$RAW_SCORE))
  results[[length(results) + 1L]] <- data.frame(
    TRAIT = trait, METHOD = "GraphPred", THRESHOLD = "ALL",
    METRIC = "max_abs_error", VALUE = max_error, REQUIRED = 1e-6,
    PASS = max_error <= 1e-6, stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, results)
write.table(results, output, sep = "\t", quote = FALSE, row.names = FALSE)
if (any(!results$PASS)) stop("PRS reproduction validation failed", call. = FALSE)
