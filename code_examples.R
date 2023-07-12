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
                width=round(runif(100,5,15)),
                xID=paste0("x",1:100),
                group=paste0("g",rep(1:2,each=50))) %>%
  as_granges()

y <- data.frame(seqnames=1,
                start=round(runif(100,0,1e4)),
                width=round(runif(100,5,15)),
                yID=paste0("y",1:100)) %>%
  as_granges()

library(tibble)
library(ggplot2)

x %>%
  anchor_5p %>%
  mutate(width=1) %>%
  add_nearest_distance(y %>% anchor_center %>% mutate(width=1)) %>%
  as_tibble() %>%
  ggplot(aes(distance, group=group, fill=group)) +
  geom_histogram(position="dodge")

x_5p <- resize(x, width=1)
y_mid <- y - ifelse(width(y) %% 2 == 0, width(y)/2-.5, floor(width(y)/2))
hits <- distanceToNearest(x_5p, y_mid)
x$distance[queryHits(hits)] <- mcols(hits)$distance
mcols(x)[,c("group","distance")] %>%
  as_tibble() %>%
  ggplot(aes(distance, group=group, fill=group)) +
  geom_histogram(position="dodge")


###############
## example 3 ##
###############

# find disjoin regions within groups of features

