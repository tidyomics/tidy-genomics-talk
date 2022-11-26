library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
g <- genes(txdb)

suppressPackageStartupMessages(library(plyranges))
g <- g %>%
  mutate(symbol = mapIds(org.Hs.eg.db, gene_id,
                         "SYMBOL", "ENTREZID"))

# for visualizing, restrict to a range of chr4
chrom <- "chr4"
rng <- c(98.8e6, 99.8e6) # where we will zoom into, 1 Mb
rng_big <- c(95e6, 100e6) # where features live, 5 Mb

# filtering the genes:
r <- data.frame(seqnames=chrom,start=rng[1]+1,end=rng[2]) %>%
  as_granges()
g %>%
  filter_by_overlaps(r) %>%
  sort() %>%
  arrange(strand)

source("boot_and_match_script.R")

library(plotgardener)
plotSomeGenes()

gr <- makeClusterRanges(chrom, rng_big, 150, 5)
seqlengths(gr) <- seqlengths(g)["chr4"]
crp <- colorRampPalette(c("dodgerblue2", "firebrick2"))
rplt <- plotRanges(gr, params=p, y=2, height=1, fill=colorby("score", palette=crp), order="random", baseline=TRUE)
textp <- pgParams(x=.1, rot=90, just="left")

plotText("original", params=textp, y=3)
shuf <- shuffle(gr, rng_big)
splt <- plotRanges(shuf, params=p, y=1, height=1, fill=colorby("score", palette=crp), order="random", baseline=TRUE)
plotText("shuffled", params=textp, y=2)

library(nullranges)
seg <- makeSegmentation(chrom, rng, g)
boot <- bootRanges(gr, blockLength=1e5, R=1, seg=seg, proportionLength=FALSE)
bplt <- plotRanges(boot, params=p, y=0, height=1, fill=colorby("score", palette=crp), order="random", baseline=TRUE)
plotText("boot", params=textp, y=1)

pageGuideHide()
