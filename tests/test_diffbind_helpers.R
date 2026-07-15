library(testthat)
source("../workflow/scripts/diffbind_helpers.R")

test_that("read_featurecounts_matrix parses coords + integer counts + sample names", {
  tf <- tempfile(fileext = ".txt")
  writeLines(c(
    "# Program:featureCounts v2.1.1",
    "Geneid\tChr\tStart\tEnd\tStrand\tLength\tresults/blacklist_filtered/S1.nobl.bam\tresults/blacklist_filtered/S2.nobl.bam",
    "peak_1\tchr1\t100\t600\t.\t500\t10\t20",
    "peak_2\tchr1\t800\t1300\t.\t500\t30\t40"
  ), tf)
  fc <- read_featurecounts_matrix(tf)
  expect_equal(fc$samples, c("S1", "S2"))
  expect_equal(dim(fc$counts), c(2L, 2L))
  expect_true(is.integer(fc$counts))
  expect_equal(fc$counts["peak_1", "S2"], 20L)
  expect_equal(fc$coords$Chr, c("chr1", "chr1"))
})

test_that("condition_from_samples maps Control->Ctrl (ref) and else->NICD3", {
  cond <- condition_from_samples(c("GSF4007-Control_1_S11", "GSF4007-NICD3-V5_1_S12"))
  expect_equal(as.character(cond), c("Ctrl", "NICD3"))
  expect_equal(levels(cond), c("Ctrl", "NICD3"))   # Ctrl is the reference
})

test_that("classify_promoter flags peaks overlapping the promoter bed", {
  pf <- tempfile(fileext = ".bed")
  writeLines(c("chr1\t1000\t2000\tP1", "chr2\t5000\t6000\tP2"), pf)   # 0-based BED
  coords <- data.frame(
    Chr   = c("chr1", "chr1", "chr2"),
    Start = c(1500,   9000,   5500),      # peak1 in P1, peak2 far, peak3 in P2
    End   = c(1600,   9100,   5600),
    stringsAsFactors = FALSE
  )
  expect_equal(classify_promoter(coords, pf), c(TRUE, FALSE, TRUE))
})

test_that("run_deseq2_group runs with manual size factors and returns result columns", {
  set.seed(1)
  samples <- c("Control_A","Control_B","Control_C","NICD3_A","NICD3_B","NICD3_C")
  np <- 200
  base <- matrix(rpois(np*6, 100), nrow = np, dimnames = list(paste0("peak_",1:np), samples))
  base[1:20, 4:6] <- base[1:20, 4:6] + 400L      # 20 peaks up in NICD3
  cond <- condition_from_samples(samples)
  sf   <- setNames(c(4.0, 2.0, 3.33, 1.67, 1.0, 5.0), samples)   # explicit per-sample size factors
  coords <- data.frame(Geneid=rownames(base), Chr="chr1",
                       Start=seq(1,by=1000,length.out=np), End=seq(500,by=1000,length.out=np))
  res <- run_deseq2_group(base, cond, sf, coords)
  expect_equal(nrow(res), np)
  expect_true(all(c("log2FoldChange","padj","baseMean","Chr","Start","End") %in% colnames(res)))
  expect_true(median(res$log2FoldChange[1:20], na.rm=TRUE) > 0)   # spiked-up peaks trend positive
})

test_that("run_deseq2_group applies the supplied size factors (not re-estimated)", {
  set.seed(2)
  samples <- c("Control_A","Control_B","Control_C","NICD3_A","NICD3_B","NICD3_C")
  base <- matrix(rpois(100*6, 100), nrow = 100, dimnames = list(paste0("p_",1:100), samples))
  cond <- condition_from_samples(samples)
  coords <- data.frame(Geneid = rownames(base), Chr = "chr1",
                       Start = seq_len(100), End = seq_len(100) + 100)
  sf1 <- setNames(c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5), samples)   # uniform
  sf2 <- setNames(c(0.2, 0.2, 0.2, 2.0, 2.0, 2.0), samples)   # condition-confounded
  r1 <- run_deseq2_group(base, cond, sf1, coords)
  r2 <- run_deseq2_group(base, cond, sf2, coords)
  # baseMean = mean of size-factor-normalized counts -> MUST differ between sf1 and sf2.
  # If run_deseq2_group re-estimated size factors (e.g. via DESeq()), it would ignore
  # the supplied sf and both runs would be identical -> this assertion would fail.
  expect_false(isTRUE(all.equal(r1$baseMean, r2$baseMean)))
})

test_that("pair_from_samples extracts the shared replicate index as the pairing factor", {
  s <- c("GSF4007-Control_1_S11","GSF4007-Control_2_S13","GSF4007-Control_3_S15",
         "GSF4007-NICD3-V5_1_S12","GSF4007-NICD3-V5_2_S14","GSF4007-NICD3-V5_3_S16")
  p <- pair_from_samples(s)
  expect_true(is.factor(p))
  expect_equal(as.character(p), c("1","2","3","1","2","3"))
})

test_that("run_deseq2_group with sf=NULL uses DESeq2 default size factors (median-of-ratios)", {
  set.seed(4)
  samples <- c("Control_A","Control_B","Control_C","NICD3_A","NICD3_B","NICD3_C")
  base <- matrix(rpois(200*6, 100), nrow = 200, dimnames = list(paste0("p_",1:200), samples))
  base[, 4:6] <- base[, 4:6] * 3L                 # a library-size difference default norm should absorb
  cond <- condition_from_samples(samples)
  coords <- data.frame(Geneid = rownames(base), Chr = "chr1",
                       Start = seq_len(200), End = seq_len(200) + 100)
  # reference: DESeq2's own median-of-ratios normalized baseMean
  dds <- DESeq2::DESeqDataSetFromMatrix(base, data.frame(condition = cond, row.names = samples), ~condition)
  dds <- DESeq2::estimateSizeFactors(dds)
  ref_bm <- rowMeans(DESeq2::counts(dds, normalized = TRUE))
  r <- run_deseq2_group(base, cond, sf = NULL, coords)
  expect_equal(nrow(r), 200)
  expect_true(all(c("log2FoldChange","padj","baseMean") %in% colnames(r)))
  expect_equal(r$baseMean, unname(ref_bm), tolerance = 1e-6)   # sf=NULL -> DESeq2 default normalization
})

test_that("run_deseq2_group accepts a pairing factor (paired design ~pair+condition)", {
  set.seed(3)
  samples <- c("Control_1_S11","Control_2_S13","Control_3_S15",
               "NICD3_1_S12","NICD3_2_S14","NICD3_3_S16")
  base <- matrix(rpois(100*6, 100), nrow = 100, dimnames = list(paste0("p_",1:100), samples))
  cond <- condition_from_samples(samples)
  pair <- pair_from_samples(samples)
  sf   <- setNames(rep(1.0, 6), samples)
  coords <- data.frame(Geneid = rownames(base), Chr = "chr1",
                       Start = seq_len(100), End = seq_len(100) + 100)
  r <- run_deseq2_group(base, cond, sf, coords, pair = pair)
  expect_equal(nrow(r), 100)
  expect_true(all(c("log2FoldChange","padj","baseMean") %in% colnames(r)))
})
