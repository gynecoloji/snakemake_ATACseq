#!/usr/bin/env Rscript
# ============================================================================
# GO enrichment for differential-openness gene sets (offline, clusterProfiler)
# ============================================================================
# Runs enrichGO separately for promoter and enhancer, for each significance
# tier (padj<0.05, p<0.01, p<0.05), and separately for UP and DOWN sets.
#
# Deliberately clusterProfiler + org.Hs.eg.db rather than Enrichr/gseapy: no
# network call at runtime, so the workflow stays reproducible offline.
#
# GATE: a set is only tested when it has more than --min-genes genes (default
# 10). Small sets give unstable, uninterpretable enrichment -- with n=3 the
# padj<0.05 tier will usually be far below the gate and is skipped by design.
#
# The universe is the coding genes reachable from ANY tested peak in the same
# analysis (written by diffopen_annotate.R), not all of org.Hs.eg.db -- using
# the whole genome as background inflates significance for ATAC-derived sets.
#
# Usage:
#   Rscript diffopen_enrich.R --genedir results/diffopen/ctcf/genes \
#       --outdir results/diffopen/ctcf/enrichment [--min-genes 10] [--ont BP]
# ============================================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
})

pa <- function(args) {
  o <- list(`min-genes` = "10", ont = "BP", `p-cut` = "0.05")
  i <- 1; while (i <= length(args)) { o[[sub("^--", "", args[i])]] <- args[i + 1]; i <- i + 2 }
  o
}

rd <- function(p) if (file.exists(p)) unique(stats::na.omit(readLines(p))) else character(0)

main <- function() {
  a <- pa(commandArgs(trailingOnly = TRUE))
  for (k in c("genedir", "outdir")) if (is.null(a[[k]])) stop("missing --", k)
  dir.create(a$outdir, showWarnings = FALSE, recursive = TRUE)
  ming <- as.integer(a$`min-genes`)

  universe <- rd(file.path(a$genedir, "universe_genes.txt"))
  message(sprintf("universe: %d coding genes", length(universe)))
  if (!length(universe)) stop("empty universe -- run diffopen_annotate.R first")

  rows <- list()
  for (cls in c("promoter", "enhancer")) {
    for (tier in c("padj05", "p01", "p05")) {
      for (dir_ in c("up", "down")) {
        genes <- rd(file.path(a$genedir, sprintf("genes_%s_%s_%s.txt", cls, tier, dir_)))
        tag <- sprintf("%s_%s_%s", cls, tier, dir_)
        if (length(genes) <= ming) {
          message(sprintf("  SKIP %-28s %d genes (<= %d)", tag, length(genes), ming))
          rows[[length(rows) + 1]] <- data.frame(set = tag, n_genes = length(genes),
                                                 tested = FALSE, n_terms = 0, top_term = NA_character_)
          next
        }
        e <- try(clusterProfiler::enrichGO(
          gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
          keyType = "SYMBOL", ont = a$ont, pAdjustMethod = "BH",
          pvalueCutoff = as.numeric(a$`p-cut`), qvalueCutoff = 0.2,
          readable = FALSE), silent = TRUE)
        if (inherits(e, "try-error") || is.null(e) || !nrow(as.data.frame(e))) {
          message(sprintf("  none %-28s %d genes, 0 enriched terms", tag, length(genes)))
          rows[[length(rows) + 1]] <- data.frame(set = tag, n_genes = length(genes),
                                                 tested = TRUE, n_terms = 0, top_term = NA_character_)
          next
        }
        df <- as.data.frame(e)
        utils::write.table(df, file.path(a$outdir, sprintf("GO_%s.tsv", tag)),
                           sep = "\t", quote = FALSE, row.names = FALSE)
        # dot plot of the top terms
        grDevices::png(file.path(a$outdir, sprintf("GO_%s.png", tag)), 1200, 900, res = 130)
        print(enrichplot::dotplot(e, showCategory = min(15, nrow(df))) +
                ggplot2::ggtitle(sprintf("GO:%s  %s  (%d genes)", a$ont, tag, length(genes))))
        grDevices::dev.off()
        message(sprintf("  OK   %-28s %d genes -> %d terms | top: %s",
                        tag, length(genes), nrow(df), df$Description[1]))
        rows[[length(rows) + 1]] <- data.frame(set = tag, n_genes = length(genes),
                                               tested = TRUE, n_terms = nrow(df),
                                               top_term = df$Description[1])
      }
    }
  }
  s <- do.call(rbind, rows)
  utils::write.table(s, file.path(a$outdir, "enrichment_summary.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  print(s)
  message("done -> ", a$outdir)
}

if (sys.nframe() == 0) main()
