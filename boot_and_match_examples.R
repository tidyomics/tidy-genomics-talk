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

# filtering the genes to this range
r <- data.frame(seqnames=chrom, start=rng[1]+1, end=rng[2]) %>%
  as_granges()

# just look at the genes in this small range
g <- g %>%
  filter_by_overlaps(r) %>%
  sort() %>%
  arrange(strand)

source("boot_and_match_script.R")

suppressPackageStartupMessages(library(plotgardener))

plotSomeGenes(chrom, rng, showGuides=FALSE)

# make n features in clumps of ~lambda
seqlens <- seqlengths(g)[chrom]
gr <- makeClusterRanges(chrom, rng_big, n=300, lambda=5, seqlens)

# define some plotting parameters for plotgardener,
# e.g. a palette for feature 'score':
pal <- colorRampPalette(c("dodgerblue2", "firebrick2"))

# shared genomic location, width & height, x position, fill, etc.
p <- pgParams(
  chrom=chrom, chromstart=rng[1], chromend=rng[2],
  width=5.5, height=1, x=.25,
  fill=colorby("score", palette=pal),
  order="random", baseline=TRUE, 
)

# shared parameters for text labels
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
# blocks 100kb, not proportion to segment length
library(nullranges)
seg <- makeSegmentation(chrom, rng_big, seqlens)
boot <- bootRanges(gr, blockLength=1e5, R=1,
                   seg=seg, proportionLength=FALSE)

# plot bootstrapped ranges
plotRanges(boot, params=p, y=0)
plotText("boot", params=textp, y=1)

# for genome-wide analysis, consider excluding gaps, repeats, etc.
# see https://dozmorovlab.github.io/excluderanges for details

#library(AnnotationHub)
#ah <- AnnotationHub()
#query(ah, "excluderanges")

###########################
## bootstrapping example ##
###########################

# first just counts as statistic

g %>%
  mutate(n_overlaps = count_overlaps(., gr))

g %>% join_overlap_left(gr) %>%
  group_by(symbol) %>%
  summarize(n_overlaps = sum(!is.na(id)))

# working with metadata

g %>% join_overlap_left(gr) %>%
  group_by(symbol) %>%
  summarize(sum_score = sum(score))

# simple violin plot

library(tibble)
library(ggplot2)
# inner instead of left: leaves out no-overlap genes
g %>% join_overlap_inner(gr) %>%
  mutate(type = "original") %>%
  group_by(symbol, type) %>%
  summarize(sum_score = sum(score)) %>%
  as_tibble() %>%
  ggplot(aes(type, sum_score)) +
  geom_violin() +
  geom_point()

# adding more draws from the distribution for simulated features

niter <- 50
sim_list <- replicate(niter, {
  makeClusterRanges(chrom, rng_big, n=300, lambda=5, seqlens)
})
sim_long <- bind_ranges(sim_list, .id="iter")

g %>% join_overlap_inner(sim_long) %>%
  mutate(type = "original") %>%
  group_by(symbol, iter, type) %>%
  summarize(sum_score = sum(score)) %>%
  as_tibble() %>%
  ggplot(aes(type, sum_score)) +
  geom_violin() +
  geom_jitter()

# shuffling and bootstrapping multiple times

shuf_list <- replicate(niter, shuffle(gr, rng_big))
shuf_long <- bind_ranges(shuf_list, .id="iter")

boot_long <- bootRanges(gr, blockLength=1e5, R=niter,
                   seg=seg, proportionLength=FALSE)

# bind together

lvls <- c("sim","shuffle","boot")
all <- bind_ranges(sim=sim_long, shuffle=shuf_long,
                   boot=boot_long, .id="type") %>%
  mutate(type = factor(type, levels=lvls))

# show table of features per iteration
head(table(all$iter, all$type))

# final plot of distributions:
# multiple draws, shuffling one instance, bootstrapping one instances

g %>% join_overlap_inner(all) %>%
  group_by(symbol, iter, type) %>%
  summarize(sum_score = sum(score)) %>%
  as_tibble() %>%
  ggplot(aes(type, sum_score)) +
  geom_violin() +
  geom_jitter(width=.25, alpha=.15)

######################
## matching example ##
######################

# start with gene plot again
plotSomeGenes(chrom, rng, showGuides=FALSE)

# make some features with particular distribution
# 1) near gene TSS, 2) tend to have large 'score' values
focal <- makeFocalFeatures(g, chrom, rng)

# 5 color palette for 'score'
pal <- colorRampPalette(c("blue","green","yellow","red"))

# new plot parameters
p <- pgParams(
  chrom=chrom, chromstart=rng[1], chromend=rng[2],
  width=5.5, height=1, x=.25,
  fill=colorby("score", palette=pal, range=c(1,5)),
  order="random", baseline=TRUE, 
)

# plot the original 'focal' GRanges
plotRanges(focal, params=p, y=2)
plotText("focal", params=textp, y=3)

# make a 'pool' of features to select from
pool <- makePool(5000, chrom, rng, seqlens)

# plot the pool (subset)
plotRanges(pool[1:200], params=p, y=1)
plotText("pool", params=textp, y=2)

# add another feature: distance to nearest TSS
tss <- g %>% anchor_5p() %>% mutate(width=1)

both <- bind_ranges(focal = focal, pool = pool, .id="type") %>%
  add_nearest_distance(tss) %>%
  mutate(log10dist = log10(distance + 1000))

hist(both$log10dist)

# three different methods for matching
m <- matchRanges(both[both$type == "focal"],
                 both[both$type == "pool"],
                 covar=~score + log10dist,
                 method="nearest", replace=TRUE)

plotCovariate(m, covar="score")
plotCovariate(m, covar="log10dist")

# plot the matched set (need to replot the others)
plotSomeGenes(chrom, rng, showGuides=FALSE)
plotRanges(focal, params=p, y=2)
plotText("focal", params=textp, y=3)
plotRanges(pool[1:200], params=p, y=1)
plotText("pool", params=textp, y=2)
plotRanges(matched(m), params=p, y=0)
plotText("matched", params=textp, y=1)
