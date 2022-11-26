library(plotgardener)
pageCreate(width=6, height=4)
chrom <- "chr4"; rng <- c(98.8e6, 99.8e6); rng_big <- c(95e6,100e6)
p <- pgParams(chrom=chrom, chromstart=rng[1], chromend=rng[2], width=5.5, x=.25)
gplt <- plotGenes(params=p, y=3, height=.75)
annoGenomeLabel(plot=gplt, x=.25, y=3.75, scale="Mb")
gr <- makeClusterRanges(chrom, rng_big, 150, 5)
crp <- colorRampPalette(c("dodgerblue2", "firebrick2"))
rplt <- plotRanges(gr, params=p, y=2, height=1, fill=colorby("score", palette=crp), order="random", baseline=TRUE)
plotText("original", x=.25, y=3, rot=90, just="left")
splt <- plotRanges(shuffle(gr, rng_big), params=p, y=1, height=1, fill=colorby("score", palette=crp), order="random", baseline=TRUE)
plotText("shuffled", x=.25, y=2, rot=90, just="left")
pageGuideHide()

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
g <- genes(txdb)
g <- g %>% mutate(symbol = mapIds(org.Hs.eg.db, gene_id, "SYMBOL", "ENTREZID"))
suppressPackageStartupMessages(library(plyranges))
r <- data.frame(seqnames=chrom,start=rng[1]+1,end=rng[2]) %>% as_granges()
g %>%
  filter_by_overlaps(r) %>%
  sort() %>%
  arrange(strand)

makeClusterRanges <- function(chrom, rng, n, lambda) {
  niter <- n/lambda
  out <- lapply(seq_len(niter), function(i) {
    nranges <- max(rpois(1, lambda), 1)
    pos <- round(runif(1, rng[1], rng[2]))
    mu <- rnorm(1, 0, 2)
    start <- pos + round(runif(nranges, -2e4, 2e4))
    score <- rnorm(nranges, mu, .5)
    data.frame(seqnames=chrom, start, width=1e4, score)
  })
  gr <- do.call(rbind, out) %>%
    as_granges() %>%
    sort()
  gr
}

shuffle <- function(gr, rng) {
  new_pos <- round(runif(length(gr),rng[1],rng[2]))
  data.frame(seqnames=seqnames(gr), start=new_pos, end=new_pos + 1e4, score=gr$score) %>%
    as_granges()
}
