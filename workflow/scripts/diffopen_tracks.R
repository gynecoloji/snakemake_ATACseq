#!/usr/bin/env Rscript
# ============================================================================
# Gviz genome-browser tracks for the top differential-openness regions
# ============================================================================
# One figure per gene (PNG + PDF) for the most UP- and most DOWN-regulated
# regions in the treatment group, showing per-sample coverage, the GENCODE gene
# models, and the differential region highlighted.
#
# Coverage comes from the depth-normalized RPGC bigWigs (results/bigwig/).
# NOTE: these are NOT re-scaled by the differential size factors, so read the
# tracks as "where the signal is", not as the quantitative effect size -- the
# log2FC in the table is the quantitative statement.
#
# Gated: only regions from a tier with more than --min-genes hits are plotted,
# matching the enrichment gate.
#
# Usage:
#   Rscript diffopen_tracks.R --genedir results/diffopen/ctcf/genes \
#       --bigwigdir results/bigwig --models results/diffopen/gene_models.rds \
#       --outdir results/diffopen/ctcf/tracks [--top 5] [--tier p01] [--pad 5000]
# ============================================================================
suppressPackageStartupMessages({
  library(Gviz); library(GenomicRanges); library(rtracklayer)
})
options(ucscChromosomeNames = FALSE)   # GENCODE/UCSC chr names, no UCSC lookup

pa <- function(args) {
  o <- list(top = "5", tier = "p01", pad = "5000", `min-genes` = "10")
  i <- 1; while (i <= length(args)) { o[[sub("^--", "", args[i])]] <- args[i + 1]; i <- i + 2 }
  o
}

#' Filename-safe gene label for a peak row (falls back to the peak id).
peak_label <- function(r) {
  g <- ifelse(is.na(r$coding_gene) | !nzchar(r$coding_gene), r$Geneid, r$coding_gene)
  gsub("[^A-Za-z0-9._-]", "_", g)
}

#' @param dense TRUE squashes the gene track into a single row. Gviz sizes are
#'   RELATIVE proportions, so a locus with more overlapping genes than the track
#'   can stack fails with "Too many stacks to draw" no matter how tall the device
#'   is; dense stacking is the only reliable escape. Used as a retry, not a
#'   default, because it loses the per-gene rows.
plot_one <- function(r, bws, ex, outdir, pad, tag, dense = FALSE, gene = NULL) {
  chr <- as.character(r$Chr); ws <- r$Start - pad; we <- r$End + pad
  if (is.null(gene)) gene <- peak_label(r)

  # common y-scale so samples are visually comparable
  rngs <- lapply(bws, function(f) {
    v <- try(rtracklayer::import(f, which = GRanges(chr, IRanges(ws, we))), silent = TRUE)
    if (inherits(v, "try-error") || !length(v)) 0 else max(v$score, na.rm = TRUE)
  })
  ymax <- max(unlist(rngs), 1)

  # sizes are accumulated alongside the tracks: an enhancer in a gene desert has
  # no exons in the window, so the GENCODE track is absent and a fixed-length
  # sizes vector would not match the trackList.
  trks <- list(Gviz::GenomeAxisTrack(cex = .8)); szs <- 1; genesz <- 0
  cols <- grDevices::hcl.colors(length(bws), "Dark 3")
  for (i in seq_along(bws)) {
    trks[[length(trks) + 1]] <- Gviz::DataTrack(
      # genome must be a character: Gviz's @genome slot rejects NA (logical).
      # ucscChromosomeNames=FALSE above means this label is never used for lookup.
      range = bws[[i]], genome = "hg38", chromosome = chr, from = ws, to = we,
      name = names(bws)[i], type = "histogram", col.histogram = cols[i],
      fill.histogram = cols[i], ylim = c(0, ymax), cex.axis = .55, cex.title = .6)
    szs <- c(szs, 2)
  }
  exw <- ex[as.character(GenomeInfoDb::seqnames(ex)) == chr &
            GenomicRanges::start(ex) <= we & GenomicRanges::end(ex) >= ws]
  if (length(exw)) {
    # collapseTranscripts="meta" merges every transcript of a gene into one
    # meta-transcript. Without it, a gene-dense window (GENCODE has 40+ isoforms
    # at some loci) needs more stacked rows than the device has height for and
    # Gviz aborts with "Too many stacks to draw".
    trks[[length(trks) + 1]] <- Gviz::GeneRegionTrack(
      exw, chromosome = chr, name = "GENCODE", transcriptAnnotation = "symbol",
      collapseTranscripts = "meta", stacking = if (dense) "dense" else "squish",
      fill = "#4a6fa5", col = NA, cex.title = .6, cex.group = .6)
    # A gene-dense window still needs one stacked row per overlapping gene even
    # after collapsing isoforms (the ACTB and CDKN2A neighbourhoods carry a dozen
    # lncRNAs). Give the track proportional height, or Gviz aborts with
    # "Too many stacks to draw".
    ngene  <- length(unique(exw$gene_name[!is.na(exw$gene_name)]))
    genesz <- max(2, min(9, ngene / 2.5))   # local: `if` does not open a scope in R
    szs <- c(szs, genesz)
  }
  trks[[length(trks) + 1]] <- Gviz::AnnotationTrack(
    GRanges(chr, IRanges(r$Start, r$End)), name = "peak",
    fill = "#b26a1b", col = NA, cex.title = .6)
  szs <- c(szs, .7)

  ttl <- sprintf("%s  %s:%s-%s  log2FC %+.2f  p=%.2g", gene, chr,
                 format(r$Start, big.mark = ","), format(r$End, big.mark = ","),
                 r$log2FoldChange, r$pvalue)
  # grow the canvas with the gene track so the extra rows actually have room
  hpx <- 1000 + round(110 * max(0, genesz - 2))
  for (dev in c("png", "pdf")) {
    f <- file.path(outdir, sprintf("gviz_%s_%s.%s", tag, gene, dev))
    if (dev == "png") grDevices::png(f, 1400, hpx, res = 150)
    else grDevices::pdf(f, width = 9.3, height = hpx / 150)
    Gviz::plotTracks(trks, from = ws, to = we, chromosome = chr, main = ttl,
                     cex.main = .85, sizes = szs)
    grDevices::dev.off()
  }
  gene
}

main <- function() {
  a <- pa(commandArgs(trailingOnly = TRUE))
  for (k in c("genedir", "bigwigdir", "models", "outdir")) if (is.null(a[[k]])) stop("missing --", k)
  dir.create(a$outdir, showWarnings = FALSE, recursive = TRUE)
  top <- as.integer(a$top); pad <- as.integer(a$pad); ming <- as.integer(a$`min-genes`)

  bwf <- sort(list.files(a$bigwigdir, pattern = "\\.bw$", full.names = TRUE))
  if (!length(bwf)) stop("no bigWigs in ", a$bigwigdir)
  names(bwf) <- sub("\\.bw$", "", basename(bwf))
  message(sprintf("%d bigWigs", length(bwf)))
  ex <- readRDS(a$models)$ex

  made <- character(0)
  for (cls in c("promoter", "enhancer")) {
    f <- file.path(a$genedir, sprintf("%s_%s_gene_annotation.tsv", cls, a$tier))
    if (!file.exists(f)) { message("skip ", cls, ": no ", basename(f)); next }
    d <- utils::read.delim(f, stringsAsFactors = FALSE)
    if (nrow(d) <= ming) {
      message(sprintf("SKIP %s/%s: %d regions (<= %d)", cls, a$tier, nrow(d), ming)); next
    }
    up <- utils::head(d[order(d$pvalue), ][d[order(d$pvalue), ]$log2FoldChange > 0, ], top)
    dn <- utils::head(d[order(d$pvalue), ][d[order(d$pvalue), ]$log2FoldChange < 0, ], top)
    # one pathological locus must not cost the whole figure set: record and move on
    # Two peaks can share a nearest gene; without make.unique the second figure
    # would silently overwrite the first.
    draw <- function(df, dir_) {
     if (!nrow(df)) return(invisible(NULL))
     labs <- make.unique(peak_label(df), sep = "_")
     for (i in seq_len(nrow(df))) {
      tg <- sprintf("%s_%s", cls, dir_)
      g  <- try(plot_one(df[i, ], bwf, ex, a$outdir, pad, tg, gene = labs[i]),
                silent = TRUE)
      if (inherits(g, "try-error")) {
        grDevices::graphics.off()          # drop the half-open device
        # gene-dense locus: retry with the gene track squashed to one row
        g <- try(plot_one(df[i, ], bwf, ex, a$outdir, pad, tg, dense = TRUE,
                          gene = labs[i]), silent = TRUE)
        if (inherits(g, "try-error")) {
          grDevices::graphics.off()
          message(sprintf("  WARN %s/%s #%d not drawn: %s", cls, dir_, i,
                          sub("\n.*", "", conditionMessage(attr(g, "condition")))))
        } else message(sprintf("  note %s/%s #%d drawn with dense gene track", cls, dir_, i))
      }
      if (!inherits(g, "try-error")) made <<- c(made, g)
     }
    }
    draw(up, "up"); draw(dn, "down")
    message(sprintf("  %s/%s: %d up + %d down attempted", cls, a$tier, nrow(up), nrow(dn)))
  }
  message("done -> ", a$outdir, "  (", length(made), " figures x png+pdf)")
}

if (sys.nframe() == 0) main()
