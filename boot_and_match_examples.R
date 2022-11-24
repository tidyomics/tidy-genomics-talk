library(plotgardener)
pageCreate(width=6, height=3)
p <- pgParams(chrom="chr1", chromstart=10e6, chromend=11e6, width=5.5)
genesPlot <- plotGenes(
    params = p, x = 0.25, y = 0.25, height = 0.75,  just = c("left", "top")
)
annoGenomeLabel(
    plot = genesPlot, x = 0.25, y = 1.0,
    scale = "Mb", just = c("left", "top")
)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
g <- genes(txdb)
g <- g %>% mutate(symbol = mapIds(org.Hs.eg.db, gene_id, "SYMBOL", "ENTREZID"))
suppressPackageStartupMessages(library(plyranges))
r <- data.frame(seqnames="chr1",start=10e6+1,end=11e6) %>% as_granges()
g %>%
  filter_by_overlaps(r) %>%
  sort() %>%
  arrange(strand)
