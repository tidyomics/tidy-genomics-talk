plotSomeGenes <- function(chrom, rng) {
  pageCreate(width=6, height=4)
  p <- pgParams(chrom=chrom, chromstart=rng[1], chromend=rng[2], width=5.5, x=.25)
  cols <- c("dodgerblue","navy")
  gplt <- plotGenes(params=p, y=3, height=.75, fill=cols, fontcolor=cols)
  annoGenomeLabel(plot=gplt, x=.25, y=3.75, scale="Mb")
}

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
  new_pos <- round(runif(length(gr), rng[1], rng[2]))
  data.frame(seqnames=seqnames(gr), start=new_pos, end=new_pos + 1e4, score=gr$score) %>%
    as_granges()
}


makeSegmentation <- function(chrom, rng, g) {
  seg <- data.frame(seqnames=chrom, start=c(1,rng[1]+1,rng[2]+1),
                    end=c(rng[1],rng[2],seqlengths(g)[[chrom]]),
                    state=c(1,2,1)) %>%
    as_granges()
}

