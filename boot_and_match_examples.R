##########################################
## visualizing genes and other features ##
##########################################

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
rng_big <- c(90e6, 110e6) # where features live, 20 Mb

# filtering the genes:
r <- data.frame(seqnames=chrom, start=rng[1]+1, end=rng[2]) %>%
  as_granges()
g <- g %>%
  filter_by_overlaps(r) %>%
  sort() %>%
  arrange(strand)

source("boot_and_match_script.R")

suppressPackageStartupMessages(library(plotgardener))
plotSomeGenes(chrom, rng, showGuides=FALSE)

# make n features in clumps of ~lambda
gr <- makeClusterRanges(chrom, rng_big, n=300, lambda=5, g)

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
seg <- makeSegmentation(chrom, rng_big, g)
boot <- bootRanges(gr, blockLength=1e5, R=1,
                   seg=seg, proportionLength=FALSE)

# plot bootstrapped ranges
plotRanges(boot, params=p, y=0)
plotText("boot", params=textp, y=1)

##############################
## bootstrapping statistics ##
##############################

# first just counts

g %>%
  mutate(n_overlaps = count_overlaps(., gr))

g %>% join_overlap_left(gr) %>%
  group_by(symbol) %>%
  summarize(n_overlaps = sum(!is.na(id)))

# working with metadata

g %>% join_overlap_left(gr) %>%
  group_by(symbol) %>%
  summarize(ave_score = mean(score))

# simple violin plot

library(tibble)
library(ggplot2)
g %>% join_overlap_left(gr) %>%
  filter(!is.na(id)) %>%
  mutate(type = "original") %>%
  group_by(symbol, type) %>%
  summarize(ave_score = mean(score)) %>%
  as_tibble() %>%
  ggplot(aes(type, ave_score)) +
  geom_violin() +
  geom_jitter()

# adding more draws from the distribution

niter <- 50
sim_list <- replicate(niter, makeClusterRanges(chrom, rng_big, n=300, lambda=5, g))
sim_long <- bind_ranges(sim_list, .id="iter")

g %>% join_overlap_left(sim_long) %>%
  filter(!is.na(id)) %>%
  mutate(type = "original") %>%
  group_by(symbol, iter, type) %>%
  summarize(ave_score = mean(score)) %>%
  as_tibble() %>%
  ggplot(aes(type, ave_score)) +
  geom_violin() +
  geom_jitter()

# shuffling and bootstrapping multiple times

shuf_list <- replicate(niter, shuffle(gr, rng))
shuf_long <- bind_ranges(shuf_list, .id="iter")

boot_long <- bootRanges(gr, blockLength=1e5, R=niter,
                   seg=seg, proportionLength=FALSE)

# bind together

lvls <- c("sim","shuffle","boot")
all <- bind_ranges(sim=sim_long, shuffle=shuf_long,
                   boot=boot_long, .id="type") %>%
  mutate(type = factor(type, levels=lvls))

head(table(all$iter, all$type))

# final plot of distributions:
# multiple draws, shuffling one instance, bootstrapping one instances

g %>% join_overlap_left(all) %>%
  filter(!is.na(id)) %>%
  group_by(symbol, iter, type) %>%
  summarize(ave_score = mean(score)) %>%
  as_tibble() %>%
  ggplot(aes(type, ave_score)) +
  geom_violin() +
  geom_jitter(width=.25, alpha=.15)
