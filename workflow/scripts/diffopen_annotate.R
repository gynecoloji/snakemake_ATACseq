#!/usr/bin/env Rscript
# ============================================================================
# Assign genes to differential-openness regions (nearest transcript TSS)
# ============================================================================
# Uses TRANSCRIPT-level TSS, not gene-level 5' ends: a gene's leftmost
# coordinate can sit far from its canonical promoter (GENCODE PGK1 spans
# 77.91-78.13 Mb, so a peak at PGK1's promoter is 193 kb from the gene START
# and mis-assigns by gene-level distance). Transcript TSS also captures
# alternative promoters.
#
# Reports the nearest TSS of ANY biotype and, separately, the nearest
# PROTEIN-CODING TSS -- enhancer nearest-genes are often lncRNAs/pseudogenes,
# so the coding column is the actionable one for enrichment.
#
# Emits, per class x significance tier, an annotation table plus up/down gene
# lists ready for enrichment, and caches the parsed gene models as RDS so the
# Gviz step does not re-parse a 1.3 GB GTF.
#
# Usage:
#   Rscript diffopen_annotate.R --indir results/diffopen/ctcf \
#       --gtf ref/gencode.v36.annotation.gtf --outdir results/diffopen/ctcf/genes \
#       [--models results/diffopen/gene_models.rds]
# ============================================================================
suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer)
})

pa <- function(args) {
  o <- list(models = "")
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    nxt <- if (i < length(args)) args[i + 1] else NA_character_
    # valueless flag (e.g. --models-only): the next token is another flag, or
    # there is none. Advancing by 2 here would swallow it and desynchronize
    # every argument that follows.
    if (is.na(nxt) || grepl("^--", nxt)) { o[[key]] <- TRUE; i <- i + 1 }
    else                                 { o[[key]] <- nxt;  i <- i + 2 }
  }
  o
}

# significance tiers, in the order they are reported
TIERS <- list(padj05 = function(d) !is.na(d$padj)   & d$padj   < 0.05,
              p01    = function(d) !is.na(d$pvalue) & d$pvalue < 0.01,
              p05    = function(d) !is.na(d$pvalue) & d$pvalue < 0.05)

main <- function() {
  a <- pa(commandArgs(trailingOnly = TRUE))
  # --models-only: just parse the GTF and cache the models, then exit. Lets the
  # workflow build the shared cache in ONE wildcard-free rule instead of racing
  # three per-mode jobs to write the same file.
  models_only <- !is.null(a$`models-only`)
  if (models_only) {
    if (is.null(a$gtf) || is.null(a$models)) stop("--models-only needs --gtf and --models")
    message("parsing GTF (transcripts + exons) ...")
    g <- rtracklayer::import(a$gtf, feature.type = c("transcript", "exon"))
    dir.create(dirname(a$models), showWarnings = FALSE, recursive = TRUE)
    saveRDS(list(tx = g[g$type == "transcript"], ex = g[g$type == "exon"]), a$models)
    message("cached -> ", a$models); return(invisible(NULL))
  }
  for (k in c("indir", "gtf", "outdir")) if (is.null(a[[k]])) stop("missing --", k)
  dir.create(a$outdir, showWarnings = FALSE, recursive = TRUE)

  # ---- gene models (parse the GTF once) ----
  cached <- nzchar(a$models) && file.exists(a$models)
  if (cached) {
    message("loading cached gene models: ", a$models)
    gm <- readRDS(a$models)
  } else {
    message("parsing GTF (transcripts + exons) ...")
    g  <- rtracklayer::import(a$gtf, feature.type = c("transcript", "exon"))
    gm <- list(tx = g[g$type == "transcript"], ex = g[g$type == "exon"])
    if (nzchar(a$models)) {
      dir.create(dirname(a$models), showWarnings = FALSE, recursive = TRUE)
      saveRDS(gm, a$models); message("cached -> ", a$models)
    }
  }
  tx  <- gm$tx
  tss <- GenomicRanges::resize(tx, width = 1, fix = "start")   # strand-aware TSS
  is_pc <- !is.na(tx$gene_type) & tx$gene_type == "protein_coding"
  tss_pc <- tss[is_pc]
  message(sprintf("  %d transcript TSS (%d protein-coding)", length(tss), length(tss_pc)))

  #' nearest-gene columns for a peak GRanges
  near <- function(pk, ref, tag) {
    h <- GenomicRanges::distanceToNearest(pk, ref, ignore.strand = TRUE)
    out <- data.frame(gene = NA_character_, dist = NA_integer_)[rep(1, length(pk)), ]
    q <- S4Vectors::queryHits(h); s <- S4Vectors::subjectHits(h)
    out$gene[q] <- ref$gene_name[s]
    out$dist[q] <- S4Vectors::mcols(h)$distance
    stats::setNames(out, paste0(tag, c("_gene", "_dist")))
  }

  summary_rows <- list()
  for (cls in c("promoter", "enhancer")) {
    f <- file.path(a$indir, sprintf("diffopen_%s.tsv", cls))
    if (!file.exists(f)) { message("skip ", cls, ": no ", f); next }
    d <- utils::read.delim(f, stringsAsFactors = FALSE)

    for (tier in names(TIERS)) {
      sel <- TIERS[[tier]](d)
      sub <- d[sel, , drop = FALSE]
      n <- nrow(sub)
      if (n == 0) {
        message(sprintf("  %s / %s: 0 regions", cls, tier))
        summary_rows[[length(summary_rows) + 1]] <-
          data.frame(class = cls, tier = tier, n_regions = 0, n_up = 0, n_down = 0,
                     n_genes_up = 0, n_genes_down = 0)
        next
      }
      pk <- GenomicRanges::GRanges(sub$Chr, IRanges::IRanges(sub$Start, sub$End))
      ann <- cbind(sub, near(pk, tss, "nearest"), near(pk, tss_pc, "coding"))
      ann <- ann[order(ann$pvalue), ]
      utils::write.table(ann, file.path(a$outdir, sprintf("%s_%s_gene_annotation.tsv", cls, tier)),
                         sep = "\t", quote = FALSE, row.names = FALSE)

      up <- unique(stats::na.omit(ann$coding_gene[ann$log2FoldChange > 0]))
      dn <- unique(stats::na.omit(ann$coding_gene[ann$log2FoldChange < 0]))
      writeLines(up, file.path(a$outdir, sprintf("genes_%s_%s_up.txt", cls, tier)))
      writeLines(dn, file.path(a$outdir, sprintf("genes_%s_%s_down.txt", cls, tier)))
      message(sprintf("  %s / %s: %d regions -> %d up / %d down coding genes",
                      cls, tier, n, length(up), length(dn)))
      summary_rows[[length(summary_rows) + 1]] <-
        data.frame(class = cls, tier = tier, n_regions = n,
                   n_up = sum(ann$log2FoldChange > 0), n_down = sum(ann$log2FoldChange < 0),
                   n_genes_up = length(up), n_genes_down = length(dn))
    }
  }

  # universe for enrichment: every coding gene reachable from any tested peak
  allpk <- do.call(rbind, lapply(c("promoter", "enhancer"), function(cls) {
    f <- file.path(a$indir, sprintf("diffopen_%s.tsv", cls))
    if (file.exists(f)) utils::read.delim(f, stringsAsFactors = FALSE)[, c("Chr", "Start", "End")]
  }))
  if (!is.null(allpk)) {
    u <- near(GenomicRanges::GRanges(allpk$Chr, IRanges::IRanges(allpk$Start, allpk$End)),
              tss_pc, "coding")
    writeLines(unique(stats::na.omit(u$coding_gene)), file.path(a$outdir, "universe_genes.txt"))
  }

  s <- do.call(rbind, summary_rows)
  utils::write.table(s, file.path(a$outdir, "annotation_summary.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  print(s)
  message("done -> ", a$outdir)
}

if (sys.nframe() == 0) main()
