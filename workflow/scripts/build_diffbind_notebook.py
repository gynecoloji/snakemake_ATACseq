import nbformat as nbf
nb = nbf.v4.new_notebook()
md = nbf.v4.new_markdown_cell
code = nbf.v4.new_code_cell

cells = []
cells.append(md(
"# Differential binding: NICD3 vs Ctrl (DESeq2)\n\n"
"DESeq2 differential binding on the workflow's **consensus peaks**, run separately on "
"**promoter** and **distal** peaks, using DESeq2's default **median-of-ratios** normalization.\n\n"
"- Counts: `results/consensus/consensus_counts.txt`\n"
"- Promoter set: `ref/promoter_chr1-22X.bed` (Ensembl reg-build); distal = non-overlapping\n"
"- **Paired design** `~pair + condition`: Ctrl/NICD3 share a replicate index (pair 1/2/3), "
"so each replicate acts as its own block\n"
"- Contrast NICD3 vs Ctrl: **positive log2FC = more open in NICD3**\n"
"- Sample **PCA** on the variance-stabilized (VST) count matrix\n"
"- Runs on `atacseq-diffbind.sif` (R kernel)."))

cells.append(code(
"suppressMessages({library(DESeq2); library(GenomicRanges); library(ggplot2)})\n"
"source('workflow/scripts/diffbind_helpers.R')\n"
"outdir <- 'results/diff_region'; dir.create(outdir, recursive=TRUE, showWarnings=FALSE)\n"
"THRESH_PADJ <- 0.05; THRESH_LFC <- 1"))

cells.append(code(
"# 1. counts + coords\n"
"fc <- read_featurecounts_matrix('results/consensus/consensus_counts.txt')\n"
"counts <- fc$counts; coords <- fc$coords\n"
"cat('peaks:', nrow(counts), ' samples:', paste(colnames(counts), collapse=', '), '\\n')"))

cells.append(code(
"# 2. condition + pairing + promoter/distal split\n"
"condition <- condition_from_samples(colnames(counts))\n"
"pair      <- pair_from_samples(colnames(counts))   # replicate index shared by each Ctrl/NICD3 pair\n"
"is_prom <- classify_promoter(coords, 'ref/promoter_chr1-22X.bed')\n"
"print(data.frame(sample=colnames(counts), condition=as.character(condition), pair=as.character(pair)))\n"
"cat('promoter peaks:', sum(is_prom), ' distal peaks:', sum(!is_prom), '\\n')"))

cells.append(code(
"# 2b. Sample PCA on the DESeq2-normalized (VST) count matrix\n"
"dds_all <- DESeqDataSetFromMatrix(counts, data.frame(condition=condition, row.names=colnames(counts)), design=~condition)\n"
"dds_all <- estimateSizeFactors(dds_all)   # DESeq2 median-of-ratios size factors\n"
"vsd <- vst(dds_all, blind=TRUE)            # variance-stabilize on the default size factors\n"
"pca <- plotPCA(vsd, intgroup='condition', returnData=TRUE)\n"
"pv <- round(100 * attr(pca, 'percentVar'))\n"
"pca_p <- ggplot(pca, aes(PC1, PC2, color=condition, label=sub('^GSF[0-9]+-', '', name))) +\n"
"  geom_point(size=3) + geom_text(vjust=-0.8, size=3, show.legend=FALSE) +\n"
"  xlab(paste0('PC1: ', pv[1], '% var')) + ylab(paste0('PC2: ', pv[2], '% var')) +\n"
"  ggtitle('Sample PCA (VST)') + theme_bw()\n"
"ggsave(file.path(outdir, 'PCA_samples.png'), pca_p, width=6, height=5, dpi=120)\n"
"print(pca_p)"))

cells.append(code(
"# 3. DESeq2 per group (median-of-ratios; paired design ~pair + condition)\n"
"res_prom  <- run_deseq2_group(counts[is_prom, , drop=FALSE],  condition, sf=NULL, coords[is_prom, ],  pair=pair)\n"
"res_dist  <- run_deseq2_group(counts[!is_prom, , drop=FALSE], condition, sf=NULL, coords[!is_prom, ], pair=pair)\n"
"write.table(res_prom, file.path(outdir,'promoter_NICD3_vs_Ctrl_deseq2.tsv'), sep='\\t', quote=FALSE, row.names=FALSE)\n"
"write.table(res_dist, file.path(outdir,'distal_NICD3_vs_Ctrl_deseq2.tsv'),   sep='\\t', quote=FALSE, row.names=FALSE)"))

cells.append(code(
"# 4. DB calls + summary\n"
"summarize_db <- function(res, label){\n"
"  sig <- subset(res, !is.na(padj) & padj < THRESH_PADJ & abs(log2FoldChange) > THRESH_LFC)\n"
"  write.table(sig, file.path(outdir, paste0(label,'_NICD3_vs_Ctrl_sig.tsv')), sep='\\t', quote=FALSE, row.names=FALSE)\n"
"  cat(sprintf('%-9s tested=%d  DB=%d  up(NICD3)=%d  down=%d\\n', label, sum(!is.na(res$padj)),\n"
"      nrow(sig), sum(sig$log2FoldChange>0), sum(sig$log2FoldChange<0)))\n"
"  sig\n"
"}\n"
"sig_prom <- summarize_db(res_prom,'promoter'); sig_dist <- summarize_db(res_dist,'distal')"))

cells.append(code(
"# 5. MA + volcano per group\n"
"ma_plot <- function(res,label){ggplot(res, aes(baseMean, log2FoldChange)) + geom_point(aes(color=!is.na(padj)&padj<THRESH_PADJ), size=.4)+\n"
"  scale_x_log10()+scale_color_manual(values=c('grey70','red'),guide='none')+labs(title=paste('MA',label))+theme_bw()}\n"
"vol_plot <- function(res,label){ggplot(res, aes(log2FoldChange, -log10(padj))) + geom_point(aes(color=!is.na(padj)&padj<THRESH_PADJ&abs(log2FoldChange)>THRESH_LFC), size=.4)+\n"
"  scale_color_manual(values=c('grey70','red'),guide='none')+geom_vline(xintercept=c(-THRESH_LFC,THRESH_LFC),lty=2)+labs(title=paste('Volcano',label))+theme_bw()}\n"
"for(g in list(list(res_prom,'promoter'), list(res_dist,'distal'))){\n"
"  ggsave(file.path(outdir,paste0('MA_',g[[2]],'.png')), ma_plot(g[[1]],g[[2]]), width=5,height=4,dpi=120)\n"
"  ggsave(file.path(outdir,paste0('volcano_',g[[2]],'.png')), vol_plot(g[[1]],g[[2]]), width=5,height=4,dpi=120)}\n"
"print(ma_plot(res_prom,'promoter')); print(vol_plot(res_dist,'distal'))"))

cells.append(code(
"# 6. nearest-gene annotation of DB regions (ChIPseeker + hg38 TxDb)\n"
"suppressMessages({library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db)})\n"
"txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene\n"
"annotate_sig <- function(sig,label){\n"
"  if(nrow(sig)==0){cat('no DB regions for',label,'\\n'); return(invisible())}\n"
"  gr <- GRanges(sig$Chr, IRanges(sig$Start, sig$End))\n"
"  an <- as.data.frame(annotatePeak(gr, TxDb=txdb, annoDb='org.Hs.eg.db', verbose=FALSE))\n"
"  out <- cbind(sig, annotation=an$annotation, SYMBOL=an$SYMBOL, distanceToTSS=an$distanceToTSS)\n"
"  write.table(out, file.path(outdir, paste0(label,'_NICD3_vs_Ctrl_sig_annotated.tsv')), sep='\\t', quote=FALSE, row.names=FALSE)\n"
"  cat(label,'annotated:',nrow(out),'regions\\n')}\n"
"annotate_sig(sig_prom,'promoter'); annotate_sig(sig_dist,'distal')"))

cells.append(md(
"## Target loci -- signal tracks\n\n"
"Gviz genome-browser tracks of ATAC signal across all 6 samples at a curated target panel, "
"with gene models from `ref/gencode.v36.annotation.gtf`. Two groups:\n\n"
"- **Notch/NICD targets & regulators** (expected to respond to NICD3): HES4, HES5, NOTCH3, HES1, "
"HEY1, HEY2, HEYL, JAG1 (=JAGGED1), ZEB1, PTEN, NRARP, MYC, DLGAP5 (HURP).\n"
"- **Housekeeping & lineage markers** (controls / context): GAPDH, TUBB, IL6, HPRT1, ACTB (beta-actin), "
"PECAM1, CDH5, VWF (=vWF), TAGLN, ACTA2.\n\n"
"Aliases were resolved to gencode symbols (JAGGED1->JAG1, vWF->VWF). **18S rRNA is omitted** -- the "
"18S/45S rDNA repeats are not in the gencode v36 primary annotation (they sit on acrocentric/unplaced "
"contigs, not chr1-22/X), so there is no gene model and no signal on the analysis BAMs. "
"**One figure per gene** from the **RPGC** depth-normalized bigWigs "
"(`results/bigwig/*.bw`, saved `<gene>_rpgc`) -- each with the 6 sample tracks (shared y-axis) "
"plus the gene model. Tracks are colored by condition "
"(Ctrl blue / NICD3 red). Read them by group: for the **Notch targets**, look for NICD3 opening "
"the promoter vs Ctrl -- if it does not open even here, the perturbation was weak (consistent with "
"DB near 0); the **housekeeping** genes (GAPDH, ACTB, ...) should be strongly and equally open in "
"all samples (a sanity check that signal/normalization is comparable), while the **lineage markers** "
"(PECAM1/CDH5/VWF endothelial, ACTA2/TAGLN smooth-muscle) give cell-identity context. "
"Gene models show **all annotated transcripts**, each on its own row and colored by "
"**transcript type** (PC = protein_coding, NMD = nonsense-mediated decay, RI = retained_intron, "
"PT = processed_transcript); exons are drawn as **arrows pointing in the direction of "
"transcription** (strand); each transcript label also carries a `(+)`/`(-)` strand tag so the "
"orientation is explicit even for transcripts too small to draw an arrow. Figures are written "
"to `results/browser_tracks/` as PNG and PDF."))

cells.append(code(
"# 7. Target loci: gene models from the GTF + per-sample bigWigs\n"
"suppressMessages({library(Gviz); library(rtracklayer)})\n"
"target_genes <- c(\n"
"  # Notch/NICD targets & regulators (aliases resolved: JAGGED1->JAG1)\n"
"  'HES4','HES5','NOTCH3','HES1','HEY1','HEY2','HEYL','JAG1','ZEB1','PTEN','NRARP','MYC','DLGAP5',\n"
"  # Housekeeping & lineage markers -- controls / context (aliases: ACTB=beta-actin, vWF->VWF)\n"
"  'GAPDH','TUBB','IL6','HPRT1','ACTB','PECAM1','CDH5','VWF','TAGLN','ACTA2')\n"
"# NOTE: '18S rRNA' is not in the gencode v36 primary annotation (rDNA repeats live on\n"
"# acrocentric/unplaced contigs, not chr1-22/X) so it has no gene model and is omitted.\n"
"gene_group <- setNames(c(rep('Notch target/regulator', 13), rep('control / lineage marker', 10)), target_genes)\n"
"gtf_path <- 'ref/gencode.v36.annotation.gtf'\n"
"trkdir <- 'results/browser_tracks'; dir.create(trkdir, recursive=TRUE, showWarnings=FALSE)\n"
"q <- '\"'\n"
"pat <- paste0('gene_name ', q, '(', paste(target_genes, collapse='|'), ')', q)\n"
"small_gtf <- tempfile(fileext='.gtf')\n"
"system(paste('grep -E', shQuote(pat), gtf_path, '>', small_gtf))   # pull just these genes (fast)\n"
"gtf_gr <- import(small_gtf); genes_gr <- gtf_gr[gtf_gr$type=='gene']\n"
"bw_rpgc <- setNames(sprintf('results/bigwig/%s.bw', colnames(counts)), colnames(counts))   # RPGC depth-normalized\n"
"cond_col <- setNames(ifelse(condition=='Ctrl','#4477aa','#ee6677'), colnames(counts))  # Ctrl blue / NICD3 red\n"
"BT_COL <- c(protein_coding='#4477aa', nonsense_mediated_decay='#ee6677',\n"
"            retained_intron='#ccbb44', processed_transcript='#228833', lncRNA='#aa3377')  # by transcript type\n"
"BT_SHORT <- c(protein_coding='PC', nonsense_mediated_decay='NMD',\n"
"              retained_intron='RI', processed_transcript='PT', lncRNA='lnc')\n"
"print(data.frame(transcript_type=names(BT_COL), short=unname(BT_SHORT[names(BT_COL)]), color=unname(BT_COL)))\n"
"print(data.frame(symbol=genes_gr$gene_name, chr=as.character(seqnames(genes_gr)),\n"
"                 start=start(genes_gr), end=end(genes_gr), strand=as.character(strand(genes_gr))))"))

cells.append(code(
"# 8. One figure per gene from the RPGC bigWigs (6 tracks + gene model) -> PNG + PDF\n"
"FLANK <- 5000\n"
"sn <- function(x) sub('^Control','Ctrl', sub('-V5','', sub('_S[0-9]+$','', sub('^GSF[0-9]+-','',x))))  # compact label\n"
"gene_model <- function(sym, chr){   # all annotated transcripts, coloured by type, arrows = strand\n"
"  ex <- gtf_gr[gtf_gr$type=='exon' & gtf_gr$gene_name==sym]\n"
"  bt <- ex$transcript_type; short <- ifelse(bt %in% names(BT_SHORT), BT_SHORT[bt], bt)\n"
"  st <- as.character(strand(ex)); arr <- ifelse(st=='+', '(+)', ifelse(st=='-', '(-)', ''))\n"
"  lbl <- paste0(ex$transcript_name, ' [', short, '] ', arr)\n"
"  grt <- GeneRegionTrack(ex, chromosome=chr, genome='hg38', name=sym, feature=bt,\n"
"           gene=ex$gene_id, exon=ex$exon_id, transcript=lbl, symbol=ex$gene_name,\n"
"           transcriptAnnotation='transcript', stacking='full', col='grey40', shape='arrow')\n"
"  list(grt=grt, n_tx=max(length(unique(ex$transcript_id)), 1))\n"
"}\n"
"plot_norm <- function(sym, bws, suffix, normlabel){\n"
"  gg <- genes_gr[genes_gr$gene_name==sym]; chr <- as.character(seqnames(gg))\n"
"  s <- start(gg)-FLANK; e <- end(gg)+FLANK; win <- GRanges(chr, IRanges(s, e))\n"
"  ymax <- max(vapply(bws, function(f){v<-import(f, which=win)$score; if(length(v)) max(v) else 0}, numeric(1)))\n"
"  ymax <- max(ceiling(ymax*1.05), 1)\n"
"  gm <- gene_model(sym, chr); grt <- gm$grt; n_tx <- gm$n_tx\n"
"  dts <- lapply(colnames(counts), function(sm) DataTrack(range=bws[[sm]], genome='hg38',\n"
"           chromosome=chr, type='histogram', name=sn(sm),\n"
"           col.histogram=cond_col[[sm]], fill.histogram=cond_col[[sm]], ylim=c(0,ymax)))\n"
"  trks <- c(list(GenomeAxisTrack()), dts, list(grt))\n"
"  sizes <- c(0.7, rep(2.4, length(dts)), max(2.0, 0.7*n_tx))   # gene model grows with transcript count\n"
"  args <- c(list(from=s, to=e, chromosome=chr, main=paste0(sym, ' [', gene_group[sym], '] -- ', normlabel), cex.main=1, sizes=sizes,\n"
"                 background.title='white', col.title='black', col.axis='black'), as.list(BT_COL))\n"
"  H <- 1050 + 40*n_tx\n"
"  png(file.path(trkdir, paste0(sym,'_',suffix,'.png')), width=1050, height=H, res=120)\n"
"  do.call(plotTracks, c(list(trks), args)); invisible(dev.off())\n"
"  pdf(file.path(trkdir, paste0(sym,'_',suffix,'.pdf')), width=1050/120, height=H/120)\n"
"  do.call(plotTracks, c(list(trks), args)); invisible(dev.off())\n"
"  do.call(plotTracks, c(list(trks), args))            # also render inline\n"
"}\n"
"present <- target_genes[target_genes %in% genes_gr$gene_name]\n"
"missing <- setdiff(target_genes, present)\n"
"if(length(missing)) cat('WARNING: no gene model found for:', paste(missing, collapse=', '), '\\n')\n"
"for(sym in present){\n"
"  plot_norm(sym, bw_rpgc, 'rpgc', 'RPGC depth-normalized')\n"
"}\n"
"cat('wrote', length(present), 'locus figures (RPGC, PNG + PDF) to', trkdir, '\\n')"))

nb.cells = cells
nb.metadata.kernelspec = {"name":"ir","display_name":"R","language":"R"}
nb.metadata.language_info = {"name":"R"}
nbf.write(nb, "ATACseq_Dx.ipynb")
print("wrote ATACseq_Dx.ipynb with", len(cells), "cells (ir kernel)")
