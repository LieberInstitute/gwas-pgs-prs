#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, required = TRUE) {
  idx <- match(name, args)
  if (is.na(idx) || idx == length(args)) {
    if (required) stop("Missing argument: ", name, call. = FALSE)
    return(NULL)
  }
  args[[idx + 1L]]
}

ct_manifest <- get_arg("--ct-manifest")
prsice_manifest <- get_arg("--prsice-manifest")
graph_manifest <- get_arg("--graph-manifest")
expected_file <- get_arg("--expected-samples")
output <- get_arg("--output")

expected <- scan(expected_file, what = character(), quiet = TRUE)
if (length(expected) != 119L || anyDuplicated(expected)) {
  stop("Expected-sample file must contain 119 unique IDs", call. = FALSE)
}

read_manifest <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

normalize_ids <- function(ids) {
  direct <- ids %in% expected
  doubled <- paste(expected, expected, sep = "_")
  mapped <- expected[match(ids, doubled)]
  ids[!direct] <- mapped[!direct]
  if (anyNA(ids) || any(!ids %in% expected)) {
    stop("Could not normalize all sample IDs", call. = FALSE)
  }
  ids
}

new_rows <- function(sample, trait, version, method, setup, threshold,
                     score, count) {
  data.frame(
    SAMPLE = sample,
    TRAIT = trait,
    GWAS_VERSION = version,
    METHOD = method,
    SETUP = setup,
    THRESHOLD = threshold,
    RAW_SCORE = as.numeric(score),
    SCORED_VARIANTS = as.integer(round(count)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

registry <- list()
append_registry <- function(x) {
  registry[[length(registry) + 1L]] <<- x
}

## collect PLINK2 SCORE1_AVG values from every legacy threshold file
ct <- read_manifest(ct_manifest)
for (i in seq_len(nrow(ct))) {
  thresholds <- strsplit(ct$THRESHOLDS[[i]], ",", fixed = TRUE)[[1L]]
  for (threshold in thresholds) {
    path <- file.path(ct$OUT_DIR[[i]], "score", paste0("score.", threshold, ".sscore"))
    if (!file.exists(path)) stop("Missing C+T score: ", path, call. = FALSE)
    score <- read.table(path, header = TRUE, check.names = FALSE,
                        stringsAsFactors = FALSE, comment.char = "")
    id_col <- intersect(c("#IID", "IID"), names(score))
    if (length(id_col) != 1L || !all(c("ALLELE_CT", "SCORE1_AVG") %in% names(score))) {
      stop("Unexpected PLINK2 score columns: ", path, call. = FALSE)
    }
    append_registry(new_rows(
      normalize_ids(score[[id_col]]), ct$TRAIT[[i]], ct$GWAS_VERSION[[i]],
      "Shizhong_CT", ct$SETUP[[i]], threshold, score$SCORE1_AVG,
      score$ALLELE_CT / 2
    ))
  }
}

## collect fixed-threshold PRSice sums and Num_SNP values
prsice <- read_manifest(prsice_manifest)
for (i in seq_len(nrow(prsice))) {
  score_path <- paste0(prsice$OUT_PREFIX[[i]], ".all_score")
  count_path <- paste0(prsice$OUT_PREFIX[[i]], ".prsice")
  score <- read.table(score_path, header = TRUE, check.names = FALSE,
                      stringsAsFactors = FALSE)
  counts <- read.table(count_path, header = TRUE, check.names = FALSE,
                       stringsAsFactors = FALSE)
  thresholds <- strsplit(prsice$THRESHOLDS[[i]], ",", fixed = TRUE)[[1L]]
  score_thresholds <- suppressWarnings(as.numeric(sub("^Pt_", "", names(score))))
  count_thresholds <- suppressWarnings(as.numeric(counts$Threshold))

  for (threshold in thresholds) {
    numeric_threshold <- as.numeric(threshold)
    score_idx <- which(!is.na(score_thresholds) &
                         abs(score_thresholds - numeric_threshold) <=
                         max(1e-15, numeric_threshold * 1e-10))
    count_idx <- which(!is.na(count_thresholds) &
                         abs(count_thresholds - numeric_threshold) <=
                         max(1e-15, numeric_threshold * 1e-10))
    if (length(score_idx) != 1L || length(count_idx) != 1L) {
      stop("Could not map PRSice threshold ", threshold, " for ", score_path,
           call. = FALSE)
    }
    append_registry(new_rows(
      normalize_ids(score$IID), prsice$TRAIT[[i]], prsice$GWAS_VERSION[[i]],
      "PRSice", prsice$SETUP[[i]], threshold, score[[score_idx]],
      rep(counts$Num_SNP[[count_idx]], nrow(score))
    ))
  }
}

## collect GraphPred score and plugin-reported matched-variant count
graph <- read_manifest(graph_manifest)
for (i in seq_len(nrow(graph))) {
  score <- read.delim(graph$OUT_FILE[[i]], check.names = FALSE,
                      stringsAsFactors = FALSE)
  if (ncol(score) < 3L) {
    stop("GraphPred output lacks --counts column: ", graph$OUT_FILE[[i]],
         call. = FALSE)
  }
  append_registry(new_rows(
    normalize_ids(score[[1L]]), graph$TRAIT[[i]], graph$GWAS_VERSION[[i]],
    "GraphPred", graph$SETUP[[i]], "ALL", score[[2L]], score[[3L]]
  ))
}

registry <- do.call(rbind, registry)
series <- interaction(registry[, c("TRAIT", "GWAS_VERSION", "METHOD", "SETUP", "THRESHOLD")],
                      drop = TRUE, lex.order = TRUE)

for (idx in split(seq_len(nrow(registry)), series)) {
  ids <- registry$SAMPLE[idx]
  if (length(ids) != 119L || anyDuplicated(ids) || !setequal(ids, expected)) {
    label <- paste(registry[idx[1L], c("TRAIT", "GWAS_VERSION", "METHOD", "SETUP", "THRESHOLD")],
                   collapse = "/")
    stop("Score series does not contain exactly the expected 119 donors: ", label,
         call. = FALSE)
  }
  if (any(!is.finite(registry$RAW_SCORE[idx]))) {
    stop("Non-finite score in registry", call. = FALSE)
  }
}

registry <- registry[order(registry$TRAIT, registry$GWAS_VERSION, registry$METHOD,
                           registry$SETUP, suppressWarnings(as.numeric(registry$THRESHOLD)),
                           registry$SAMPLE), ]
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(registry, output, sep = "\t", quote = FALSE, row.names = FALSE)
