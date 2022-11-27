plotSomeGenes <- function(chrom, rng, showGuides) {
  pageCreate(width=6, height=4, showGuides=showGuides)
  p <- pgParams(chrom=chrom, chromstart=rng[1], chromend=rng[2], width=5.5)
  cols <- c("dodgerblue","navy")
  gplt <- plotGenes(params=p, x=.25, y=3, height=.75, fill=cols, fontcolor=cols)
  annoGenomeLabel(plot=gplt, x=.25, y=3.75, scale="Mb")
}

makeClusterRanges <- function(chrom, rng, n, lambda, seqlens) {
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
    sort() %>%
    mutate(id = seq_along(.))
  seqlengths(gr) <- seqlens
  gr
}

shuffle <- function(gr, rng, width=1e4) {
  new_pos <- round(runif(length(gr), rng[1], rng[2]))
  data.frame(seqnames=seqnames(gr), start=new_pos, end=new_pos + width,
             score=gr$score, id=gr$id) %>%
    as_granges()
}


makeSegmentation <- function(chrom, rng, seqlens) {
  seg <- data.frame(seqnames=chrom, start=c(1,rng[1]+1,rng[2]+1),
                    end=c(rng[1],rng[2],seqlens),
                    state=c(1,2,1)) %>%
    as_granges()
}

makeFocalFeatures <- function(g, chrom, rng) {
  tss <- g %>%
    anchor_5p() %>%
    mutate(width = 1e4) %>%
    select(-c(gene_id, symbol))
  bind_ranges(replicate(3, tss)) %>%
    shift(round(runif(length(.), -1e4, 1e4))) %>%
    mutate(
      score = factor(sample(3:5, length(.), replace=TRUE, prob=1:3), levels=1:5)
    ) %>%
    unname()
}

makePool <- function(n, chrom, rng, seqlens) {
  gr <- data.frame(
    seqnames=chrom,
    start=round(runif(n, rng[1], rng[2])), width=1e4,
    score = factor(sample(1:5, n, replace=TRUE),levels=1:5)) %>%
    as_granges()
  seqlengths(gr) <- seqlens
  gr
}
