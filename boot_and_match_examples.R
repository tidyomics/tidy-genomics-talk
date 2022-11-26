library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
g <- genes(txdb)

# add symbols to genes
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

suppressPackageStartupMessages(library(plotgardener))
plotSomeGenes(chrom, rng, showGuides=FALSE)

# make n features in clumps of ~lambda
gr <- makeClusterRanges(chrom, rng_big, n=150, lambda=5, g)

# define some plotting parameters 
pal <- colorRampPalette(c("dodgerblue2", "firebrick2"))
p <- pgParams(
  chrom=chrom, chromstart=rng[1], chromend=rng[2], width=5.5,
  x=.25, fill=colorby("score", palette=pal),
  order="random", baseline=TRUE, height=1
)
textp <- pgParams(x=.1, rot=90, just="left")

# plot the original GRanges
plotRanges(gr, params=p, y=2)
plotText("original", params=textp, y=3)

# uniform shuffling
shuf <- shuffle(gr, rng_big)

# plot shuffled ranges
plotRanges(shuf, params=p, y=1)
plotText("shuffled", params=textp, y=2)

# segmented block bootstrapping
library(nullranges)
seg <- makeSegmentation(chrom, rng, g)
boot <- bootRanges(gr, blockLength=1e5, R=1,
                   seg=seg, proportionLength=FALSE)

# plot bootstrapped ranges
plotRanges(boot, params=p, y=0)
plotText("boot", params=textp, y=1)
