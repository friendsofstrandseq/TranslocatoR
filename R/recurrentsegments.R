#' Get recurrent segments
#' 
#' @param counts.seg segmented count matrix (segfinder.R output)
#' 
#' @return a list of recurrent segments that are not whole chromosomes and at least 1Mb long
#' @author Alex van Vliet
#' 
#' @export

rec.seg <- function(counts.seg) {
  
  # identify recurrent (occur more than once) breakpoints that cover at least 1Mb and aren't whole chromosomes
  bps.recseg <- counts.seg[,.N, by=.(chrom, start, end)][mixedorder(chrom)]
  bps.recseg[, whole.chrom := start == min(start) & end == max(end), by=chrom]
  bps.recseg[, len := end - start]
  bps.recseg <- bps.recseg[whole.chrom==F & N > 1 & len >= 1000000]
  
  return(bps.recseg)
}
