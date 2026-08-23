options(stringsAsFactors = FALSE)

## resolve repository-relative inputs so the script can run from any directory
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(script_arg)) dirname(normalizePath(script_arg)) else getwd()
repo_dir <- normalizePath(file.path(script_dir, ".."))
gwas_dir <- Sys.getenv("GWAS_DIR", file.path(repo_dir, "GWAS"))
out_dir <- Sys.getenv("OUT_DIR", file.path(repo_dir, "generated", "gwas_summary"))
bcftools <- Sys.getenv("BCFTOOLS", "bcftools")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## these are the exact European public panels used in the manuscript analysis
panels <- data.frame(
  disorder = c("BD", "MDD", "SCZD"),
  study = c("PGC_BIP2024", "PGC_MDD2025", "PGC_SCZ2022"),
  ancestry = "EUR",
  has_23andme = FALSE,
  n_case = c(59287L, 412305L, 53386L),
  n_control = c(781022L, 1588397L, 77258L),
  panel_file = c(
    "BPD/bip2024_eur_no23andMe.hg38.bcf",
    "MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf",
    "SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf"
  ),
  paper_doi = c(
    "10.1038/s41586-024-08468-9",
    "10.1016/j.cell.2024.12.002",
    "10.1038/s41586-022-04434-5"
  ),
  release_doi = c(
    "10.6084/m9.figshare.27216117",
    "10.6084/m9.figshare.27061255",
    "10.6084/m9.figshare.19426775"
  )
)
panels$n_total <- panels$n_case + panels$n_control
panels$path <- file.path(gwas_dir, panels$panel_file)

if (!all(file.exists(panels$path))) {
  stop("missing GWAS BCF: ", paste(panels$path[!file.exists(panels$path)], collapse = ", "))
}

record_count <- function(path) {
  ## count every BCF record without applying local FILTER or p-value criteria
  value <- system2(bcftools, c("index", "-n", path), stdout = TRUE, stderr = FALSE)
  if (length(value) != 1L || !grepl("^[0-9]+$", value)) {
    stop("could not count BCF records: ", path)
  }
  as.numeric(value)
}

read_exploratory <- function(path) {
  ## one query supplies both the 1e-5 and nested 5e-8 exact-variant sets
  tmp <- tempfile(fileext = ".tab")
  on.exit(unlink(tmp))
  fmt <- "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%LP]\\n"
  status <- system2(
    bcftools,
    c("query", "-i", shQuote("FORMAT/LP>=5"), "-f", shQuote(fmt), shQuote(path)),
    stdout = tmp,
    stderr = FALSE
  )
  if (!identical(status, 0L)) stop("bcftools query failed: ", path)

  x <- read.delim(
    tmp,
    header = FALSE,
    col.names = c("chrom", "pos", "ref", "alt", "lp"),
    colClasses = c("character", "integer", "character", "character", "character")
  )

  ## expand multiallelic records so ALT and Number=A LP values remain paired
  alts <- strsplit(x$alt, ",", fixed = TRUE)
  lps <- strsplit(x$lp, ",", fixed = TRUE)
  if (any(lengths(alts) != lengths(lps))) {
    stop("ALT and LP allele counts differ: ", path)
  }
  idx <- rep(seq_len(nrow(x)), lengths(alts))
  x <- data.frame(
    chrom = x$chrom[idx],
    pos = x$pos[idx],
    ref = x$ref[idx],
    alt = unlist(alts, use.names = FALSE),
    lp = suppressWarnings(as.numeric(unlist(lps, use.names = FALSE)))
  )
  x <- x[is.finite(x$lp), ]
  x$variant_id <- paste(x$chrom, x$pos, x$ref, x$alt, sep = ":")

  ## retain one value per canonical GRCh38 allele if a panel contains duplicates
  lp <- tapply(x$lp, x$variant_id, max)
  data.frame(variant_id = names(lp), lp = as.numeric(lp))
}

variants <- setNames(lapply(panels$path, read_exploratory), panels$disorder)
panels$n_variant_all <- vapply(panels$path, record_count, numeric(1))
panels$n_variant_1e5 <- vapply(variants, function(x) sum(x$lp >= 5), numeric(1))
panels$n_variant_5e8 <- vapply(
  variants,
  function(x) sum(x$lp >= -log10(5e-8)),
  numeric(1)
)

## build inclusive and mutually exclusive exact-variant overlap counts
overlap_row <- function(scheme, p_bd, p_mdd, p_sczd) {
  cutoffs <- -log10(c(BD = p_bd, MDD = p_mdd, SCZD = p_sczd))
  sets <- Map(function(x, cutoff) x$variant_id[x$lp >= cutoff], variants, cutoffs)
  ids <- Reduce(union, sets)
  member <- vapply(sets, function(x) ids %in% x, logical(length(ids)))
  bd <- member[, "BD"]
  mdd <- member[, "MDD"]
  sczd <- member[, "SCZD"]

  data.frame(
    scheme = scheme,
    p_bd = p_bd,
    p_mdd = p_mdd,
    p_sczd = p_sczd,
    n_bd = sum(bd),
    n_mdd = sum(mdd),
    n_sczd = sum(sczd),
    bd_only = sum(bd & !mdd & !sczd),
    mdd_only = sum(!bd & mdd & !sczd),
    sczd_only = sum(!bd & !mdd & sczd),
    bd_mdd_only = sum(bd & mdd & !sczd),
    bd_sczd_only = sum(bd & !mdd & sczd),
    mdd_sczd_only = sum(!bd & mdd & sczd),
    all_three = sum(bd & mdd & sczd),
    bd_mdd = sum(bd & mdd),
    bd_sczd = sum(bd & sczd),
    mdd_sczd = sum(mdd & sczd),
    union = length(ids)
  )
}

overlaps <- do.call(rbind, list(
  overlap_row("strict_all", 5e-8, 5e-8, 5e-8),
  overlap_row("mood_explore", 1e-5, 1e-5, 5e-8),
  overlap_row("explore_all", 1e-5, 1e-5, 1e-5)
))

## fail if exclusive regions do not reconstruct the set union and totals
exclusive_cols <- c(
  "bd_only", "mdd_only", "sczd_only", "bd_mdd_only",
  "bd_sczd_only", "mdd_sczd_only", "all_three"
)
stopifnot(rowSums(overlaps[exclusive_cols]) == overlaps$union)
stopifnot(overlaps$bd_mdd == overlaps$bd_mdd_only + overlaps$all_three)
stopifnot(overlaps$bd_sczd == overlaps$bd_sczd_only + overlaps$all_three)
stopifnot(overlaps$mdd_sczd == overlaps$mdd_sczd_only + overlaps$all_three)

metrics <- panels[c(
  "disorder", "study", "ancestry", "has_23andme", "n_total", "n_case",
  "n_control", "n_variant_all", "n_variant_1e5", "n_variant_5e8",
  "panel_file", "paper_doi", "release_doi"
)]

legend <- data.frame(
  table = c(
    rep("metrics", 10), rep("overlaps", 10), rep("scheme", 3)
  ),
  term = c(
    "n_total", "n_case", "n_control", "has_23andme", "n_variant_all",
    "n_variant_1e5", "n_variant_5e8", "panel_file", "paper_doi", "release_doi",
    "p_bd,p_mdd,p_sczd", "n_bd,n_mdd,n_sczd", "*_only", "*_mdd_only",
    "*_sczd_only", "all_three", "bd_mdd", "bd_sczd", "mdd_sczd", "union",
    "strict_all", "mood_explore", "explore_all"
  ),
  description = c(
    "European participants represented by cases plus controls",
    "European cases", "European controls", "whether the panel includes 23andMe",
    "all GRCh38 BCF records with no local p-value or FILTER criterion",
    "unique exact variants with p <= 1e-5", "unique exact variants with p <= 5e-8",
    "repository-relative GRCh38 BCF", "study publication DOI", "public release DOI",
    "per-disorder p-value cutoffs", "per-disorder variant-set sizes",
    "variants exclusive to one disorder", "two-disorder region excluding SCZD",
    "two-disorder region excluding the remaining mood disorder",
    "variants shared by BD, MDD, and SCZD", "inclusive BD and MDD intersection",
    "inclusive BD and SCZD intersection", "inclusive MDD and SCZD intersection",
    "unique variants present in one or more sets", "5e-8 for all disorders",
    "1e-5 for BD and MDD; 5e-8 for SCZD", "1e-5 for all disorders"
  )
)

write_tab <- function(x, name) {
  write.table(
    x,
    file.path(out_dir, name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

write_tab(metrics, "gwas_panel_metrics.tab")
write_tab(overlaps, "gwas_variant_overlaps.tab")
write_tab(legend, "gwas_table_legend.tab")

message("wrote GWAS summary tables to ", out_dir)
