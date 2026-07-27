options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: export_prsice_model_tables.R <export_dir>")

export_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- file.path(export_dir, "out")
checks_dir <- file.path(export_dir, "checks")
dir.create(checks_dir, showWarnings = FALSE, recursive = TRUE)

labels <- c("BPD", "MDD")

read_hashes <- function(path) {
  if (!file.exists(path)) return(setNames(character(), character()))
  x <- read.table(path, header = FALSE, sep = "", col.names = c("hash", "path"),
                  stringsAsFactors = FALSE)
  stats::setNames(x$hash, x$path)
}

hash_for <- function(hashes, pattern) {
  hits <- hashes[grepl(pattern, names(hashes))]
  if (length(hits) == 0) return(NA_character_)
  unname(hits[[1]])
}

read_commands <- function(path) {
  if (!file.exists(path)) return(setNames(character(), character()))
  x <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
  stats::setNames(x$command, x$disorder)
}

read_version <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  x <- readLines(path, warn = FALSE)
  hit <- grep("^[0-9][.][0-9][.][0-9]", x, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

read_export_version <- function(export_dir) {
  primary <- read_version(file.path(export_dir, "logs", "tool_versions.txt"))
  if (!is.na(primary)) return(primary)
  read_version(file.path(export_dir, "logs", "prsice_version.txt"))
}

read_ws <- function(path, col_classes = NA) {
  read.table(path, header = TRUE, sep = "", check.names = FALSE, quote = "",
             comment.char = "", fill = TRUE, colClasses = col_classes)
}

zfile <- function(path, mode = "wt") {
  gzfile(path, open = mode)
}

find_snp_file <- function(prefix) {
  hits <- Sys.glob(paste0(prefix, "*.snp*"))
  hits <- hits[!grepl("\\.(log|prsice|all_score|best|summary)$", hits)]
  if (length(hits) == 0) stop("no PRSice SNP file found for prefix: ", prefix)
  hits[which.min(nchar(hits))]
}

threshold_name <- function(x) {
  paste0("in_", gsub("\\+", "", x, fixed = FALSE))
}

read_prsice_counts <- function(path) {
  x <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                  colClasses = "character")
  x$Num_SNP <- as.integer(x$Num_SNP)
  x
}

read_retained_snps <- function(path) {
  x <- read_ws(path)
  snp_col <- if ("SNP" %in% names(x)) "SNP" else names(x)[1]
  unique(as.character(x[[snp_col]]))
}

write_table_gz <- function(x, path) {
  con <- zfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

compare_scores <- function(label) {
  old_path <- file.path("prsice", "out", paste0(label, ".all_score"))
  new_path <- file.path(out_dir, paste0(label, ".all_score"))
  old <- read_ws(old_path)
  new <- read_ws(new_path)
  id_col <- if ("IID" %in% names(old)) "IID" else names(old)[1]
  score_cols <- intersect(grep("^Pt_", names(old), value = TRUE),
                          grep("^Pt_", names(new), value = TRUE))
  old <- old[order(old[[id_col]]), c(id_col, score_cols)]
  new <- new[order(new[[id_col]]), c(id_col, score_cols)]
  if (!identical(old[[id_col]], new[[id_col]])) {
    stop("sample IDs differ in all_score comparison for ", label)
  }
  diffs <- abs(as.matrix(old[score_cols]) - as.matrix(new[score_cols]))
  data.frame(disorder = label,
             compared_scores = length(score_cols),
             max_abs_diff = max(diffs, na.rm = TRUE))
}

manifest_rows <- list()
score_checks <- list()
input_hashes <- read_hashes(file.path(export_dir, "logs", "input_hashes.sha256"))
commands <- read_commands(file.path(export_dir, "logs", "prsice_commands.tsv"))
prsice_version <- read_export_version(export_dir)
prsice_binary_sha256 <- hash_for(input_hashes, "PRSice_linux$")

for (label in labels) {
  prefix <- file.path(out_dir, label)
  base_path <- file.path("prsice", "base", paste0(label, ".base.tsv.gz"))
  prsice_path <- paste0(prefix, ".prsice")
  snp_path <- find_snp_file(prefix)

  base <- read_ws(base_path)
  keep <- data.frame(SNP = read_retained_snps(snp_path))
  model <- merge(keep, base, by = "SNP", all.x = TRUE, sort = FALSE)
  missing_base <- sum(!stats::complete.cases(model[c("CHR", "BP", "A1", "A2", "BETA", "P")]))
  if (missing_base != 0) stop(label, ": retained SNPs missing from base file: ", missing_base)
  model <- model[c("CHR", "BP", "SNP", "A1", "A2", "BETA", "P")]
  model$disorder <- label
  model <- model[c("disorder", "CHR", "BP", "SNP", "A1", "A2", "BETA", "P")]

  counts <- read_prsice_counts(prsice_path)
  for (i in seq_len(nrow(counts))) {
    thr_text <- counts$Threshold[i]
    thr <- as.numeric(thr_text)
    col <- threshold_name(thr_text)
    model[[col]] <- as.numeric(model$P) <= thr
  }

  post_path <- file.path(out_dir, paste0(label, ".postclump_model.tsv.gz"))
  write_table_gz(model, post_path)

  threshold_rows <- vector("list", nrow(counts))
  for (i in seq_len(nrow(counts))) {
    thr_text <- counts$Threshold[i]
    col <- threshold_name(thr_text)
    y <- model[model[[col]], c("disorder", "CHR", "BP", "SNP", "A1", "A2", "BETA", "P")]
    y$threshold <- thr_text
    y <- y[c("disorder", "threshold", "CHR", "BP", "SNP", "A1", "A2", "BETA", "P")]
    threshold_rows[[i]] <- y
    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      disorder = label,
      threshold = thr_text,
      prsice_version = prsice_version,
      prsice_binary_sha256 = prsice_binary_sha256,
      base_sha256 = hash_for(input_hashes, paste0("prsice/base/", label, "[.]base[.]tsv[.]gz$")),
      prsice_command = unname(commands[[label]]),
      prsice_num_snp = counts$Num_SNP[i],
      model_num_snp = nrow(y),
      postclump_model = post_path,
      threshold_model = file.path(out_dir, paste0(label, ".threshold_model.tsv.gz")),
      prsice_snp_file = snp_path,
      stringsAsFactors = FALSE
    )
  }
  threshold_model <- do.call(rbind, threshold_rows)
  threshold_path <- file.path(out_dir, paste0(label, ".threshold_model.tsv.gz"))
  write_table_gz(threshold_model, threshold_path)

  score_checks[[length(score_checks) + 1L]] <- compare_scores(label)
}

manifest <- do.call(rbind, manifest_rows)
write.table(manifest, file.path(out_dir, "model_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

score_check <- do.call(rbind, score_checks)
write.table(score_check, file.path(checks_dir, "score_reproducibility.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

count_check <- manifest[c("disorder", "threshold", "prsice_num_snp", "model_num_snp")]
count_check$match <- count_check$prsice_num_snp == count_check$model_num_snp
write.table(count_check, file.path(checks_dir, "threshold_count_check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

if (any(!count_check$match)) stop("threshold counts do not match PRSice .prsice output")
if (any(score_check$max_abs_diff > 1e-8)) stop("rerun scores differ from previous PRSice outputs")
