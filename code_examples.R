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

