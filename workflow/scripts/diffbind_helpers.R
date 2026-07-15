# Helpers for DESeq2 differential binding (ATACseq_Dx.ipynb).

#' Read a featureCounts table -> coords + integer count matrix (cols = sample names).
read_featurecounts_matrix <- function(path) {
  df <- read.delim(path, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  meta <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  count_cols <- setdiff(colnames(df), meta)
  samples <- sub("\\.nobl\\.bam$", "", basename(count_cols))
  counts <- as.matrix(df[, count_cols, drop = FALSE])
  # featureCounts with --countReadPairs yields integer counts; coerce to integer for DESeq2.
  storage.mode(counts) <- "integer"
  colnames(counts) <- samples
  rownames(counts) <- df$Geneid
  list(coords = df[, c("Geneid", "Chr", "Start", "End")], counts = counts, samples = samples)
}

#' Condition factor: samples whose name contains "Control" -> Ctrl (reference), else NICD3.
condition_from_samples <- function(samples) {
  cond <- ifelse(grepl("Control", samples), "Ctrl", "NICD3")
  factor(cond, levels = c("Ctrl", "NICD3"))
}

#' Pairing factor: the replicate index a Control/NICD3 pair shares -- the digits
#' before the trailing _S<lane> in the sample name (e.g. *_1_S11 / *_1_S12 -> "1").
pair_from_samples <- function(samples) {
  factor(sub(".*_([0-9]+)_S[0-9]+$", "\\1", samples))
}

#' TRUE where a peak overlaps any promoter interval from a BED file.
classify_promoter <- function(coords, promoter_bed_path) {
  pb <- read.delim(promoter_bed_path, header = FALSE, stringsAsFactors = FALSE)
  prom  <- GenomicRanges::GRanges(pb[[1]], IRanges::IRanges(pb[[2]] + 1L, pb[[3]]))  # BED 0-based -> 1-based
  peaks <- GenomicRanges::GRanges(coords$Chr, IRanges::IRanges(coords$Start, coords$End))
  IRanges::overlapsAny(peaks, prom)
}

#' DESeq2 NICD3-vs-Ctrl. With `sf = NULL` (default) DESeq2 uses its own
#' median-of-ratios size-factor estimate; pass explicit per-sample `sf` to apply
#' them verbatim instead.
#' If `pair` (a per-sample blocking factor) is given, the design becomes
#' `~pair + condition` (paired: control for pair-to-pair variation, test condition).
run_deseq2_group <- function(counts, condition, sf = NULL, coords, pair = NULL) {
  stopifnot(identical(rownames(counts), coords$Geneid))
  coldata <- data.frame(condition = condition, row.names = colnames(counts))
  if (is.null(pair)) {
    design <- ~condition
  } else {
    coldata$pair <- pair
    design <- ~pair + condition
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = design)
  if (is.null(sf)) {
    dds <- DESeq2::estimateSizeFactors(dds)                # default: median-of-ratios
  } else {
    DESeq2::sizeFactors(dds) <- sf[colnames(counts)]       # explicit factors by name; do NOT estimate
  }
  dds <- DESeq2::estimateDispersions(dds)
  dds <- DESeq2::nbinomWaldTest(dds)
  res    <- DESeq2::results(dds, contrast = c("condition", "NICD3", "Ctrl"))
  shrunk <- DESeq2::lfcShrink(dds, coef = "condition_NICD3_vs_Ctrl", type = "apeglm")
  data.frame(coords,
             baseMean       = res$baseMean,
             log2FoldChange = shrunk$log2FoldChange,
             lfcSE          = shrunk$lfcSE,
             pvalue         = res$pvalue,
             padj           = res$padj,
             row.names = NULL, check.names = FALSE)
}
