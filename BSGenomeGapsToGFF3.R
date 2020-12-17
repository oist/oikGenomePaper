library("BSgenome.Odioica.local.OKI2018.I69")

getGaps <- function(chr, BSgenome) {
  m <- matchPattern("N", BSgenome[[chr]])
  m <- reduce(m)
  if (length(m) == 0)
    return (GRanges())
  GRanges(chr, IRanges(start(m), end(m)), seqinfo=seqinfo(BSgenome))        
}

chrs <- seqlevels(OKI2018_I69)

ranges <- sapply(chrs, getGaps, OKI2018_I69)

ranges <- unlist(as(ranges, "GRangesList"))

rtracklayer::export.gff3(ranges, "~/OKI2018_I69.gaps.gff3")
