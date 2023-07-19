# code examples in tidy format (plyranges) and base Bioconductor
# Michael Love
# July 12 2023

###############
## example 1 ##
###############

# first example is from the plyranges paper Figure 3
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1597-8
# "an overlap and aggregate operation that returns the same result"

library(plyranges)

# create data in R rather than reading in BED files

gwas <- data.frame(seqnames=1,
                   start=round(runif(100,0,100)),
                   width=1, rsID=paste0("rs",1:100)) %>%
  as_granges()

exons <- data.frame(seqnames=1,
                    start=round(runif(100,0,100)),
                    width=5, exonID=paste0("e",1:100)) %>%
  as_granges()

# tidy

res1 <- exons %>%
  join_overlap_inner(gwas) %>%
  group_by(rsID) %>%
  summarise(n = n_distinct(exonID))

# base bioc

hits <- findOverlaps(exons, gwas, ignore.strand = FALSE)
olap <- splitAsList(exons$exonID[queryHits(hits)], gwas$rsID[subjectHits(hits)])
n <- lengths(unique(olap))
res2 <- DataFrame(rsID = names(n), n = as.integer(n))


all.equal(res1, res2)

###############
## example 2 ##
###############

# distance from one set of features (5p) to nearest other set (center)
# group by the type of features and plot histogram

x <- data.frame(seqnames=1,
                start=round(runif(100,0,1e4)),
                width=round(runif(100,5,15))) %>%
  as_granges() %>%
  sort()
x <- x %>%
  mutate(xID = paste0("x",1:100),
         group = paste0("g",rep(1:2,each=50)))

y <- data.frame(seqnames=1,
                start=round(runif(100,0,1e4)),
                width=round(runif(100,5,15))) %>%
  as_granges() %>%
  sort()
y <- y %>%
  mutate(yID=paste0("y",1:100))

library(tibble)
library(ggplot2)

# tidy

x %>%
  anchor_5p() %>%
  mutate(width=1) %>%
  add_nearest_distance(y %>% anchor_center %>% mutate(width=1)) %>%
  as_tibble() %>%
  ggplot(aes(distance, group=group, fill=group)) +
  geom_histogram(position="dodge")

# base bioc

x_5p <- resize(x, width=1)
y_mid <- y - ifelse(width(y) %% 2 == 0, width(y)/2-.5, floor(width(y)/2))
hits <- distanceToNearest(x_5p, y_mid)
x$distance[queryHits(hits)] <- mcols(hits)$distance
df <- as.data.frame(mcols(x)[,c("group","distance")])
ggplot(df, aes(distance, group=group, fill=group)) +
  geom_histogram(position="dodge")


###############
## example 3 ##
###############

# find disjoint regions within groups of features, filter to the overlapping pieces

# tidy

x %>%
  join_overlap_inner(range(x) %>%
                     tile_ranges(width=1000) %>%
                     mutate(tile=seq_along(.))) %>%
  group_by(tile) %>%
  disjoin_ranges(total = n()) %>%
  filter(total > 1)

# base bioc

tiles <- tile(range(x), width=1000)[[1]]
tiles$tile <- seq_along(tiles)
hits <- findOverlaps(x, tiles)
res <- lapply(1:length(tiles), function(t) {
  x_sub <- x[queryHits(hits)[subjectHits(hits) == t]]
  d <- disjoin(x_sub)
  cov <- as(coverage(x_sub), "GRanges")
  d[d %over% cov[cov$score > 1]]
})
do.call(c, res)
