#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: export_prsice_model_tables.R <model_manifest.tsv>")
}

manifest_path <- args[[1]]
manifest <- read.table(manifest_path, header = TRUE, sep = "\t",
                       quote = "", comment.char = "", check.names = FALSE)
required <- c("label", "base_file", "prsice_prefix", "output_dir")
missing <- setdiff(required, names(manifest))
if (length(missing)) stop("manifest missing columns: ", paste(missing, collapse = ","))

read_ws <- function(path) {
  read.table(path, header = TRUE, sep = "", quote = "", comment.char = "",
             check.names = FALSE)
}

write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

find_snp_file <- function(prefix) {
  hits <- Sys.glob(paste0(prefix, "*.snp*"))
  hits <- hits[!grepl("\\.(log|prsice|all_score|best|summary)$", hits)]
  if (!length(hits)) stop("no PRSice SNP file found for prefix: ", prefix)
  hits[which.min(nchar(hits))]
}

threshold_name <- function(x) paste0("in_", gsub("[^A-Za-z0-9.]", "", x))

manifest_rows <- list()
idx <- 0L

for (i in seq_len(nrow(manifest))) {
  label <- manifest$label[i]
  base_file <- manifest$base_file[i]
  prefix <- manifest$prsice_prefix[i]
  output_dir <- manifest$output_dir[i]
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  base <- read_ws(base_file)
  prsice <- read.table(paste0(prefix, ".prsice"), header = TRUE, sep = "\t",
                       quote = "", comment.char = "", check.names = FALSE)
  snp_file <- find_snp_file(prefix)
  snps <- read_ws(snp_file)
  snp_col <- if ("SNP" %in% names(snps)) "SNP" else names(snps)[1]
  keep <- data.frame(SNP = unique(as.character(snps[[snp_col]])))

  model <- merge(keep, base, by = "SNP", all.x = TRUE, sort = FALSE)
  needed <- c("CHR", "BP", "A1", "A2", "BETA", "P")
  if (any(!stats::complete.cases(model[needed]))) {
    stop(label, ": retained SNPs missing base annotations")
  }
  model <- model[c("CHR", "BP", "SNP", "A1", "A2", "BETA", "P")]
  model$label <- label
  model <- model[c("label", "CHR", "BP", "SNP", "A1", "A2", "BETA", "P")]

  for (j in seq_len(nrow(prsice))) {
    thr <- as.character(prsice$Threshold[j])
    col <- threshold_name(thr)
    model[[col]] <- as.numeric(model$P) <= as.numeric(thr)
  }

  post_path <- file.path(output_dir, paste0(label, ".postclump_model.tsv.gz"))
  write_tsv_gz(model, post_path)

  threshold_rows <- vector("list", nrow(prsice))
  for (j in seq_len(nrow(prsice))) {
    thr <- as.character(prsice$Threshold[j])
    col <- threshold_name(thr)
    y <- model[model[[col]], c("label", "CHR", "BP", "SNP", "A1", "A2", "BETA", "P")]
    y$threshold <- thr
    y <- y[c("label", "threshold", "CHR", "BP", "SNP", "A1", "A2", "BETA", "P")]
    threshold_rows[[j]] <- y
    idx <- idx + 1L
    manifest_rows[[idx]] <- data.frame(
      label = label,
      threshold = thr,
      prsice_num_snp = as.integer(prsice$Num_SNP[j]),
      model_num_snp = nrow(y),
      postclump_model = post_path,
      threshold_model = file.path(output_dir, paste0(label, ".threshold_model.tsv.gz")),
      prsice_snp_file = snp_file,
      stringsAsFactors = FALSE
    )
  }
  write_tsv_gz(do.call(rbind, threshold_rows),
               file.path(output_dir, paste0(label, ".threshold_model.tsv.gz")))
}

out_manifest <- do.call(rbind, manifest_rows)
write.table(out_manifest, file.path(unique(manifest$output_dir)[1], "model_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

count_check <- out_manifest[c("label", "threshold", "prsice_num_snp", "model_num_snp")]
count_check$match <- count_check$prsice_num_snp == count_check$model_num_snp
write.table(count_check, file.path(unique(manifest$output_dir)[1], "threshold_count_check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
if (any(!count_check$match)) stop("threshold counts do not match PRSice output")
