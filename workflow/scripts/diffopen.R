#!/usr/bin/env Rscript
# ============================================================================
# Differential chromatin openness with a SELECTABLE normalization mode
# ============================================================================
# Runs DESeq2 on the consensus fragment-count matrix, differing only in how the
# per-sample size factors are established:
#
#   --mode none     Standard DESeq2 median-of-ratios over ALL consensus peaks.
#                   Baseline. Assumes most regions do not change; a true
#                   genome-wide shift is absorbed and reported as zero.
#
#   --mode ctcf     Median-of-ratios restricted to constitutive CTCF anchors
#                   (peaks overlapping the CTCF cCRE BED). Spike-in free.
#                   Robust to composition changes, but like `none` it defines the
#                   global level as invariant, so a true global shift is
#                   undetectable by construction.
#
#   --mode rnastable  Median-of-ratios restricted to promoter-class consensus
#                   peaks that also sit over the TSS window of an RNA-seq-STABLE
#                   gene (neither up- nor down-regulated). Spike-in free, keyed to
#                   the sample's own transcriptome. Like `none`/`ctcf` it defines
#                   the anchor set as invariant, so a true global shift is
#                   undetectable by construction.
#
# NOTE: this is a spike-in-free workflow, so EVERY mode defines its anchor set as
# invariant and therefore cannot detect a genuine uniform genome-wide shift in
# accessibility — that needs an external spike-in. Read results as RELATIVE
# changes, and compare modes before trusting any single one.
#
# Inputs (produced by this workflow):
#   results/consensus/consensus_counts.txt     featureCounts matrix + coords
#   config/samples.csv                         design (sample_id, type, group)
#   ref/constitutive_ctcf_hg38.bed             mode=ctcf only
#
# Usage:
#   Rscript workflow/scripts/diffopen.R --mode ctcf \
#       --counts  results/consensus/consensus_counts.txt \
#       --samples config/samples.csv \
#       --ctcf    ref/GRCh38-cCREs.CTCF-only.bed \
#       --outdir  results/diffopen/ctcf [--ref-label Control]
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

# ---- tiny --key value arg parser -------------------------------------------
parse_args <- function(args) {
  out <- list(mode = "none", `ref-label` = "Control", `min-anchors` = "200",
              `trim-k` = "2.5", `trim-iter` = "2", `min-class-peaks` = "100",
              `tss-window` = "2000", `rna-gene-col` = "gene",
              `rna-lfc-col` = "log2FoldChange", `rna-padj-col` = "padj",
              `rna-basemean-col` = "baseMean", `rna-basemean-min` = "10",
              `rna-padj-min` = "0.5", `rna-lfc-max` = "0.5",
              `promoter-class-required` = "true")
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i]); out[[key]] <- args[i + 1]; i <- i + 2
  }
  out
}

# ---- IO ----------------------------------------------------------------------

#' featureCounts table -> list(coords, counts[G x n], samples). Cols named by sample.
read_featurecounts_matrix <- function(path) {
  df <- read.delim(path, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  meta <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  count_cols <- setdiff(colnames(df), meta)
  samples <- sub("\\.nobl\\.bam$", "", basename(count_cols))
  counts <- as.matrix(df[, count_cols, drop = FALSE]); storage.mode(counts) <- "integer"
  colnames(counts) <- samples; rownames(counts) <- df$Geneid
  list(coords = df[, c("Geneid", "Chr", "Start", "End")], counts = counts, samples = samples)
}

#' samples.csv -> data.frame(sample_id, condition[factor, ref first], pair[factor]).
read_design <- function(path, ref_label) {
  s <- read.csv(path, stringsAsFactors = FALSE)
  cond_raw <- s$type
  lv <- c(ref_label, setdiff(unique(cond_raw), ref_label))     # reference level first
  condition <- factor(cond_raw, levels = lv)
  pair <- factor(sub(".*_([0-9]+)_S[0-9]+$", "\\1", s$sample_id))
  data.frame(sample_id = s$sample_id, condition = condition, pair = pair,
             stringsAsFactors = FALSE)
}

#' Read a DESeq2/edgeR results table (TSV by default; CSV by extension).
read_de_table <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE))
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  else
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Transcript GRanges (with gene_name) from the cached models RDS, or by parsing
#' the GTF if no cache is supplied. Reuses the diffopen_gene_models cache so the
#' 1.3 GB GTF is not re-parsed.
load_gene_models_tx <- function(models, gtf) {
  if (!is.null(models) && nzchar(models) && file.exists(models)) {
    gm <- readRDS(models)
    return(gm$tx)
  }
  if (is.null(gtf)) stop("rnastable: need --models (cache) or --gtf")
  suppressPackageStartupMessages({ library(rtracklayer) })
  rtracklayer::import(gtf, feature.type = "transcript")
}

#' DESeq2/edgeR results -> character vector of transcriptionally STABLE genes.
#' Stable = expressed AND not significant AND small fold change:
#'   baseMean >= basemean_min  AND  (padj >= padj_min OR padj is NA)  AND
#'   abs(log2FoldChange) <= lfc_max
#' NA padj counts as non-significant; NA baseMean/log2FC never qualifies.
stable_genes_from_de <- function(de, gene_col = "gene", lfc_col = "log2FoldChange",
                                 padj_col = "padj", basemean_col = "baseMean",
                                 basemean_min = 10, padj_min = 0.5, lfc_max = 0.5) {
  for (col in c(gene_col, lfc_col, padj_col, basemean_col))
    if (!col %in% colnames(de)) stop("RNA-seq DE table missing column: ", col)
  bm   <- suppressWarnings(as.numeric(de[[basemean_col]]))
  padj <- suppressWarnings(as.numeric(de[[padj_col]]))
  lfc  <- suppressWarnings(as.numeric(de[[lfc_col]]))
  keep <- !is.na(bm)  & bm >= basemean_min &
          (is.na(padj) | padj >= padj_min) &
          !is.na(lfc)  & abs(lfc) <= lfc_max
  genes <- unique(as.character(de[[gene_col]][keep]))
  genes[!is.na(genes) & nzchar(genes)]
}

#' Logical over consensus peaks: TRUE where the peak overlaps any feature in a BED.
bed_overlap <- function(coords, bed_path) {
  suppressPackageStartupMessages({
    library(GenomicRanges); library(IRanges)
  })
  b <- read.delim(bed_path, header = FALSE, stringsAsFactors = FALSE)
  feat  <- GRanges(b[[1]], IRanges(b[[2]] + 1L, b[[3]]))           # BED 0-based -> 1-based
  peaks <- GRanges(coords$Chr, IRanges(coords$Start, coords$End))
  overlapsAny(peaks, feat)
}

#' Backwards-compatible alias used by the CTCF size-factor path.
ctcf_overlap <- function(coords, ctcf_path) bed_overlap(coords, ctcf_path)

#' Classify each consensus peak as promoter / enhancer / other.
#'
#' PROMOTER PRECEDENCE: a peak overlapping both a promoter and an enhancer
#' feature is called a promoter. Without this, peaks straddling an adjacent
#' promoter+enhancer pair would be double-counted across the two classes.
classify_peaks <- function(coords, promoter_bed, enhancer_bed) {
  prom <- bed_overlap(coords, promoter_bed)
  enh  <- bed_overlap(coords, enhancer_bed) & !prom
  factor(ifelse(prom, "promoter", ifelse(enh, "enhancer", "other")),
         levels = c("promoter", "enhancer", "other"))
}

#' Strand-aware TSS +/- window for transcripts whose gene_name is stable.
#' Reuses the same TSS definition as diffopen_annotate.R (resize to start).
stable_tss_windows <- function(tx, stable_genes, window) {
  suppressPackageStartupMessages({ library(GenomicRanges) })
  keep <- !is.na(tx$gene_name) & tx$gene_name %in% stable_genes
  tss  <- GenomicRanges::resize(tx[keep], width = 1, fix = "start")
  GenomicRanges::resize(tss, width = 2L * as.integer(window) + 1L, fix = "center")
}

#' Row indices of consensus peaks that anchor the rnastable normalization.
#' A peak anchors iff it overlaps a stable-gene TSS window AND (when required)
#' is promoter-class. `promoter_is` is a logical over the coords rows.
rnastable_anchor_idx <- function(coords, promoter_is, tss_windows,
                                 promoter_class_required = TRUE) {
  suppressPackageStartupMessages({ library(GenomicRanges); library(IRanges) })
  peaks <- GRanges(coords$Chr, IRanges(coords$Start, coords$End))
  over  <- IRanges::overlapsAny(peaks, tss_windows)
  if (promoter_class_required) which(promoter_is & over) else which(over)
}

# ---- size factors, one function per mode ------------------------------------

#' Median-of-ratios restricted to a subset of rows (the CTCF anchors).
#' sf_i = exp( median_g [ log(K_gi + 0.5) - mean_j log(K_gj + 0.5) ] ), geomean 1.
median_of_ratios <- function(counts, idx) {
  lc  <- log(counts[idx, , drop = FALSE] + 0.5)
  ref <- rowMeans(lc)                                  # per-anchor reference
  sf  <- exp(apply(lc - ref, 2, stats::median))
  sf / exp(mean(log(sf)))
}

#' CTCF size factors with an iterative *invariance* trim.
#'
#' The anchors are constitutive in ENCODE, but a given experiment can still move
#' a few of them (a CTCF site inside a responsive enhancer, a copy-number
#' difference, a blacklist-adjacent artifact) and a handful of large movers can
#' drag a median-of-ratios estimate. So: estimate -> measure how far each anchor
#' actually shifted between the two conditions -> drop anchors more than
#' `trim_k` MADs from the median shift -> re-estimate on the survivors.
#'
#' Caveat: this is self-referential, so it cannot detect a genuine *uniform*
#' genome-wide shift (neither can plain median-of-ratios). No mode in this
#' spike-in-free workflow can — that needs an external spike-in.
#'
#' @return list(sf, idx = surviving anchors, n_start)
size_factors_ctcf <- function(counts, anchor_idx, condition,
                              trim_k = 2.5, iter = 2, min_anchors = 200) {
  idx <- anchor_idx
  sf  <- median_of_ratios(counts, idx)
  lv  <- levels(droplevels(condition))
  if (length(lv) < 2) return(list(sf = sf, idx = idx, n_start = length(anchor_idx)))

  isA <- condition == lv[1]          # reference level (contrast denominator)
  isB <- condition == lv[2]          # contrast numerator
  if (!any(isA) || !any(isB))
    return(list(sf = sf, idx = idx, n_start = length(anchor_idx)))

  for (it in seq_len(max(0L, iter))) {
    # depth-corrected signal at the current anchors
    ln    <- sweep(log2(counts[idx, , drop = FALSE] + 0.5), 2, log2(sf), "-")
    delta <- rowMeans(ln[, isB, drop = FALSE]) - rowMeans(ln[, isA, drop = FALSE])
    med   <- stats::median(delta)
    s_mad <- stats::mad(delta)
    if (!is.finite(s_mad) || s_mad == 0) break
    keep <- abs(delta - med) <= trim_k * s_mad
    if (sum(keep) < min_anchors || all(keep)) break     # don't over-trim / nothing to do
    idx <- idx[keep]
    sf  <- median_of_ratios(counts, idx)
  }
  list(sf = sf, idx = idx, n_start = length(anchor_idx))
}

#' rnastable size factors: median-of-ratios over consensus peaks that are
#' promoter-class AND over a stable gene's TSS window, with the same iterative
#' invariance trim as the ctcf mode. Refuses below `min_anchors`.
size_factors_rnastable <- function(counts, coords, promoter_is, tx, de,
                                   gene_col, lfc_col, padj_col, basemean_col,
                                   basemean_min, padj_min, lfc_max,
                                   window, min_anchors, promoter_class_required,
                                   condition, trim_k, iter) {
  stable <- stable_genes_from_de(de, gene_col, lfc_col, padj_col, basemean_col,
                                 basemean_min, padj_min, lfc_max)
  present    <- intersect(stable, unique(stats::na.omit(as.character(tx$gene_name))))
  match_rate <- if (length(stable)) length(present) / length(stable) else NA_real_
  windows <- stable_tss_windows(tx, stable, window)
  idx     <- rnastable_anchor_idx(coords, promoter_is, windows, promoter_class_required)
  n_anchor <- length(idx)
  if (n_anchor < min_anchors)
    stop(sprintf(paste0("only %d rnastable anchors (< %d) [promoter-class peaks over ",
                        "stable-gene TSS windows]. Loosen thresholds (raise --rna-lfc-max ",
                        "or lower --rna-padj-min), widen --tss-window, or set ",
                        "--promoter-class-required false."), n_anchor, min_anchors))
  fit <- size_factors_ctcf(counts, idx, condition,
                           trim_k = trim_k, iter = iter, min_anchors = min_anchors)
  list(sf = fit$sf, n_stable = length(stable), match_rate = match_rate,
       n_anchor = n_anchor, n_kept = length(fit$idx))
}

# ---- DESeq2 -----------------------------------------------------------------

#' Write a results table plus pre-filtered nominal-significance subsets.
#'
#' With n=3 the per-peak FDR is very conservative, so the interpretable signal
#' lives in the nominal-p sets (read their up/down DIRECTION balance, not their
#' size -- the counts can fall below chance expectation). p<0.01 is the tighter
#' set typically carried into gene annotation / enrichment.
write_results <- function(tbl, outdir, stem) {
  wr <- function(x, f) utils::write.table(x, file.path(outdir, f), sep = "\t",
                                          quote = FALSE, row.names = FALSE)
  wr(tbl, paste0(stem, ".tsv"))
  for (thr in c(0.05, 0.01)) {
    sub <- tbl[!is.na(tbl$pvalue) & tbl$pvalue < thr, , drop = FALSE]
    sub <- sub[order(sub$pvalue), , drop = FALSE]
    wr(sub, sprintf("%s_nominal_p%02d.tsv", stem, round(thr * 100)))
  }
}

#' DESeq2 size factors as DESeq2 itself would estimate them (mode=none), computed
#' ONCE on the full matrix so every peak class shares them.
size_factors_deseq2 <- function(counts, condition) {
  cd  <- data.frame(condition = condition, row.names = colnames(counts))
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, cd, ~condition)
  DESeq2::sizeFactors(DESeq2::estimateSizeFactors(dds))
}

#' Fit one peak class with the GLOBAL size factors injected (never re-estimated:
#' size factors are a library-level property, so promoter/enhancer/all must share
#' them). Each class still gets its own dispersion trend and its own within-class
#' FDR, which is the point of splitting.
#'
#' @param idx row indices of the class (NULL = all peaks)
fit_class <- function(counts, coords, condition, pair, size_factors, idx = NULL,
                      label = "all") {
  if (!is.null(idx)) {
    counts <- counts[idx, , drop = FALSE]
    coords <- coords[idx, , drop = FALSE]
  }
  coldata <- data.frame(condition = condition, pair = pair, row.names = colnames(counts))
  # Use the paired design only when every pair is seen exactly once per condition.
  paired <- nlevels(droplevels(pair)) > 1 &&
            all(table(coldata$pair, coldata$condition) == 1)
  design <- if (paired) ~pair + condition else ~condition
  message(sprintf("[%s] %d peaks | design %s%s", label, nrow(counts),
                  deparse(design), if (paired) " (paired)" else " (unpaired)"))

  dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design)
  DESeq2::sizeFactors(dds) <- size_factors
  dds <- DESeq2::estimateDispersions(dds)
  dds <- DESeq2::nbinomWaldTest(dds)

  cf  <- grep("^condition_", DESeq2::resultsNames(dds), value = TRUE)[1]
  res <- DESeq2::results(dds, name = cf)
  shr <- tryCatch(DESeq2::lfcShrink(dds, coef = cf, type = "apeglm"),
                  error = function(e) { message("lfcShrink failed; using unshrunk LFC"); res })
  list(
    # Both the shrunk and the raw effect size are reported. apeglm's posterior
    # can collapse -- here it returned an lfcSE constant to 4 decimals across
    # every quantile of a 105k-peak matrix, squashing log2FoldChange to ~0 for
    # the bulk while leaving the top hits intact. The _MLE columns are the
    # unshrunk GLM estimate, which is what PyDESeq2 reports and what should be
    # used for any effect-size threshold or cross-tool comparison. The p-values
    # come from the unshrunk Wald test either way and are unaffected.
    table = data.frame(coords,
                       baseMean           = res$baseMean,
                       log2FoldChange     = shr$log2FoldChange,
                       lfcSE              = shr$lfcSE,
                       log2FoldChange_MLE = res$log2FoldChange,
                       lfcSE_MLE          = res$lfcSE,
                       stat               = res$stat,
                       pvalue             = res$pvalue,
                       padj               = res$padj,
                       row.names = NULL, check.names = FALSE),
    contrast = cf,
    n        = nrow(counts),
    n_sig    = sum(res$padj < 0.05, na.rm = TRUE),
    n_nom    = sum(res$pvalue < 0.05, na.rm = TRUE),
    n_nom01  = sum(res$pvalue < 0.01, na.rm = TRUE),
    med_lfc  = stats::median(shr$log2FoldChange, na.rm = TRUE),
    up_frac  = mean(shr$log2FoldChange[which(res$pvalue < 0.05)] > 0, na.rm = TRUE)
  )
}

# ---- main -------------------------------------------------------------------

main <- function() {
  a <- parse_args(commandArgs(trailingOnly = TRUE))
  mode <- a$mode
  if (!mode %in% c("none", "ctcf", "rnastable"))
    stop("--mode must be one of: none, ctcf, rnastable (got '", mode, "')")
  for (k in c("counts", "samples", "outdir"))
    if (is.null(a[[k]])) stop("missing required --", k)
  if (mode == "ctcf"    && is.null(a$ctcf))    stop("--mode ctcf requires --ctcf")
  if (mode == "rnastable") {
    if (is.null(a$`rna-table`)) stop("--mode rnastable requires --rna-table")
    if (is.null(a$models) && is.null(a$gtf))
      stop("--mode rnastable requires --models (cache) or --gtf")
    if (is.null(a$`promoter-bed`) || is.null(a$`enhancer-bed`))
      stop("--mode rnastable requires --promoter-bed and --enhancer-bed")
  }
  dir.create(a$outdir, showWarnings = FALSE, recursive = TRUE)

  fc  <- read_featurecounts_matrix(a$counts)
  des <- read_design(a$samples, a$`ref-label`)

  samp <- des$sample_id
  stopifnot(all(samp %in% colnames(fc$counts)))
  counts <- fc$counts[, samp, drop = FALSE]

  if (nlevels(droplevels(des$condition)) < 2)
    stop("need >= 2 conditions in samples.csv 'type' column for a differential test")

  # ---- size factors for the requested mode ----
  # if/else rather than switch(): `<<-` inside a switch branch would assign to the
  # global env, not this frame, silently leaving the counters at NA.
  n_anchors  <- NA_integer_
  n_kept     <- NA_integer_
  rna_stable <- NA_integer_
  rna_match  <- NA_real_
  if (mode == "none") {
    # Estimated once on the FULL matrix, then shared by every peak class.
    sf <- size_factors_deseq2(counts, des$condition)
  } else if (mode == "ctcf") {
    idx <- which(ctcf_overlap(fc$coords, a$ctcf))
    n_anchors <- length(idx)
    message(sprintf("CTCF-overlapping consensus peaks: %d / %d", n_anchors, nrow(counts)))
    if (n_anchors < as.integer(a$`min-anchors`))
      stop(sprintf("only %d CTCF anchors (< %s); refusing to normalize on so few",
                   n_anchors, a$`min-anchors`))
    fit <- size_factors_ctcf(counts, idx, des$condition,
                             trim_k      = as.numeric(a$`trim-k`),
                             iter        = as.integer(a$`trim-iter`),
                             min_anchors = as.integer(a$`min-anchors`))
    sf     <- fit$sf
    n_kept <- length(fit$idx)
    message(sprintf("anchors kept after invariance trim: %d / %d", n_kept, n_anchors))
  } else if (mode == "rnastable") {
    tx  <- load_gene_models_tx(a$models, a$gtf)
    de  <- read_de_table(a$`rna-table`)
    prom_is <- classify_peaks(fc$coords, a$`promoter-bed`, a$`enhancer-bed`) == "promoter"
    rr <- size_factors_rnastable(
      counts, fc$coords, prom_is, tx, de,
      gene_col = a$`rna-gene-col`, lfc_col = a$`rna-lfc-col`,
      padj_col = a$`rna-padj-col`, basemean_col = a$`rna-basemean-col`,
      basemean_min = as.numeric(a$`rna-basemean-min`),
      padj_min = as.numeric(a$`rna-padj-min`),
      lfc_max = as.numeric(a$`rna-lfc-max`),
      window = as.integer(a$`tss-window`),
      min_anchors = as.integer(a$`min-anchors`),
      promoter_class_required = identical(tolower(a$`promoter-class-required`), "true"),
      condition = des$condition,
      trim_k = as.numeric(a$`trim-k`), iter = as.integer(a$`trim-iter`))
    sf         <- rr$sf
    n_anchors  <- rr$n_anchor
    n_kept     <- rr$n_kept
    rna_stable <- rr$n_stable
    rna_match  <- rr$match_rate
    message(sprintf("rnastable: %d stable genes (%.1f%% matched GTF) -> %d anchor peaks -> %d kept",
                    rna_stable, 100 * rna_match, n_anchors, n_kept))
  }

  # ---- fit: pooled first, then promoter / enhancer separately ----
  # Each class is fit on its own so it gets its own dispersion trend and its own
  # within-class FDR; the size factors above are shared, never re-estimated.
  fits <- list(all = fit_class(counts, fc$coords, des$condition, des$pair, sf,
                               NULL, "all"))
  res  <- fits$all$table

  cls <- NULL
  if (!is.null(a$`promoter-bed`) && !is.null(a$`enhancer-bed`)) {
    cls <- classify_peaks(fc$coords, a$`promoter-bed`, a$`enhancer-bed`)
    message("peak classes: ",
            paste(sprintf("%s=%d", levels(cls), as.integer(table(cls))), collapse = "  "))
    for (k in c("promoter", "enhancer")) {
      idx <- which(cls == k)
      if (length(idx) < as.integer(a$`min-class-peaks`)) {
        message(sprintf("skipping %s: only %d peaks", k, length(idx))); next
      }
      fits[[k]] <- fit_class(counts, fc$coords, des$condition, des$pair, sf, idx, k)
    }
  }

  # ---- outputs ----
  write_results(res, a$outdir, "differential_openness")
  for (k in c("promoter", "enhancer")) {
    if (!is.null(fits[[k]]))
      write_results(fits[[k]]$table, a$outdir, sprintf("diffopen_%s", k))
  }
  sft <- data.frame(sample = names(sf), size_factor = as.numeric(sf),
                    row.names = NULL, check.names = FALSE)
  utils::write.table(sft, file.path(a$outdir, "size_factors.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  # MA plot
  grDevices::png(file.path(a$outdir, "MA_plot.png"), 1100, 850, res = 130)
  sig <- !is.na(res$padj) & res$padj < 0.05
  plot(log2(res$baseMean + 1), res$log2FoldChange, pch = 16, cex = 0.3, col = "#9aa8a5",
       xlab = "log2 mean normalized count", ylab = "log2 fold change",
       main = sprintf("MA — normalization: %s", mode))
  graphics::points(log2(res$baseMean[sig] + 1), res$log2FoldChange[sig],
                   pch = 16, cex = 0.4, col = "#b26a1b")
  graphics::abline(h = 0, col = "#c0c0c0", lty = 3)
  grDevices::dev.off()

  # ---- summary (incl. the size-factor spread diagnostic) ----
  spread <- max(sf) / min(sf)
  class_lines <- character(0)
  if (!is.null(cls)) {
    tb <- table(cls)
    class_lines <- c(
      "",
      sprintf("peak classes (promoter precedence): promoter=%d  enhancer=%d  other=%d",
              tb[["promoter"]], tb[["enhancer"]], tb[["other"]]),
      "  class             n    padj<0.05   p<0.05   p<0.01   median log2FC   %% up (p<0.05)")
    for (k in c("all", "promoter", "enhancer")) {
      f <- fits[[k]]
      if (is.null(f)) next
      class_lines <- c(class_lines,
        sprintf("  %-10s %7d %8d %8d %8d %14.4f %12.1f%%",
                k, f$n, f$n_sig, f$n_nom, f$n_nom01, f$med_lfc, 100 * f$up_frac))
    }
    class_lines <- c(class_lines,
      "",
      "  ---- how to read this table -------------------------------------------",
      "  n            peaks in the class. Promoter/enhancer are assigned by overlap",
      "               with the Ensembl Regulatory Build BEDs, PROMOTER PRECEDENCE",
      "               (a peak hitting both is counted as promoter, never twice).",
      "               Each class is fit separately, so it gets its own dispersion",
      "               trend and its own within-class FDR. Size factors are NOT",
      "               re-estimated per class -- they are a library-level property,",
      "               computed once on all peaks and shared.",
      "",
      "  padj<0.05    Benjamini-Hochberg significant WITHIN the class. With n=3 per",
      "               condition the per-peak power is very low, so expect few or",
      "               none; this is the only column you may quote as 'significant'.",
      "",
      "  p<0.05 /     Nominal, NOT multiple-testing corrected. Do not report these",
      "  p<0.01       as significant sites, and do not read their COUNT as evidence:",
      "               DESeq2's p-values are conservative at n=3, so the count can sit",
      "               BELOW the ~5%% expected by chance. They exist to expose direction.",
      "",
      "  %% up         Of the p<0.05 peaks, the fraction with log2FoldChange > 0, i.e.",
      "               MORE OPEN IN TREATMENT. This is a DIRECTION BALANCE, not a",
      "               significance measure. Under the null it sits near 50%%.",
      "               Interpretation:",
      "                 ~50%%        no coherent directional program",
      "                 70-95%%      a real, coordinated shift -- and it should get",
      "                             STRONGER as you tighten p<0.05 -> p<0.01 (noise",
      "                             regresses to 50%%, true signal sharpens)",
      "                 exactly 100%% SUSPECT. A perfect one-directional split across",
      "                             independent classes is the fingerprint of a global",
      "                             scaling artifact (a bad size-factor set), not",
      "                             biology. Cross-check the size-factor spread above.",
      "",
      "  Compare modes before trusting any single one: `none` and `ctcf` should",
      "  broadly agree; where they disagree, that axis is normalization-dependent",
      "  and should not be reported as a finding.",
      "  ------------------------------------------------------------------------")
  }
  summ <- c(
    sprintf("normalization mode        : %s", mode),
    sprintf("contrast                  : %s", fits$all$contrast),
    sprintf("consensus peaks           : %d", nrow(counts)),
    if (mode == "ctcf")
      sprintf("CTCF anchors              : %d overlapping -> %d kept after invariance trim (%.1f%% dropped)",
              n_anchors, n_kept, 100 * (1 - n_kept / n_anchors)) else NULL,
    if (mode == "rnastable") c(
      sprintf("stable genes (RNA-seq)    : %d (%.1f%% matched GTF gene_name)", rna_stable, 100 * rna_match),
      sprintf("rnastable anchors         : %d promoter-class over stable TSS -> %d kept after trim",
              n_anchors, n_kept),
      "note                      : like none/ctcf, assumes anchors invariant; cannot see a uniform global shift"
    ) else NULL,
    sprintf("size-factor spread (max/min): %.2fx", spread),
    sprintf("median |log2FC|           : %.4f", stats::median(abs(res$log2FoldChange), na.rm = TRUE)),
    sprintf("median log2FC (global tilt): %.4f", stats::median(res$log2FoldChange, na.rm = TRUE)),
    sprintf("differential (padj<0.05)  : %d", sum(res$padj < 0.05, na.rm = TRUE)),
    sprintf("nominal (pvalue<0.05)     : %d", sum(res$pvalue < 0.05, na.rm = TRUE)),
    class_lines,
    "",
    "Interpretation: a large size-factor spread that tracks condition means the",
    "normalization is confounded -- compare the modes before trusting any one of them.")
  writeLines(summ, file.path(a$outdir, "run_summary.txt"))
  cat(paste(summ, collapse = "\n"), "\n")
  message("done -> ", a$outdir)
}

if (sys.nframe() == 0) main()
