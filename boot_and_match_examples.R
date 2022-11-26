library(plotgardener)
pageCreate(width=6, height=3)
chrom <- "chr4"; rng <- c(98800000, 99700000)
p <- pgParams(chrom=chrom, chromstart=rng[1], chromend=rng[2], width=5.5)
genesPlot <- plotGenes(params=p, x=.25, y=2, height=.75)
annoGenomeLabel(plot=genesPlot, x=.25, y=2.75, scale="Mb")
gr <- makeClusterRanges(chrom, rng)
crp <- colorRampPalette(c("dodgerblue2", "firebrick2"))
plotRanges(gr, params=p, x=.25, y=1.5, height=.5, fill=colorby("score", palette=crp))


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

makeClusterRanges <- function(chrom, rng, n=30, lambda=5) {
  niter <- n/lambda
  out <- lapply(seq_len(niter), function(i) {
    nranges <- max(rpois(1, lambda), 1)
    pos <- round(runif(1, rng[1]+5e4, rng[2]-5e4))
    mu <- rnorm(1, 0, 2)
    start <- pos + round(runif(nranges, -1e4, 1e4))
    score <- rnorm(nranges, mu, .2)
    data.frame(seqnames=chrom, start, width=1e4, score)
  })
  gr <- do.call(rbind, out) %>%
    as_granges() %>%
    sort()
  gr
}
