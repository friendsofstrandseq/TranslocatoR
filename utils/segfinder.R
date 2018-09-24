#' Segment finder for TranslocatoR
#' 
#' @param m count matrix or segment matrix with sample, cell, chrom, start, end and class
#'
#' @return a list of breakpoints where strand state changes take place
#' @author Alex van Vliet

seg.finder <- function(m) {
  
  # data in correct order for identifying consecutive bins (i.e. adjacent bins that have the same state)
  m <- m[order(sample, cell, chrom, start, end)]
  # next line is from Venla's SCE detector
  m$cnsc <- cumsum(m[, .(consecutive = c(1,abs(diff(as.numeric(factor(class)))))), by = .(sample, cell, chrom)]$consecutive)
  
  # gather data in segments: take start, end and state of segments
  m.seg <- m[, .SD[c(1,.N)], by=cnsc]
  m.seg[, c("start","end"):=.(min(start), max(end)), by=cnsc]
  m.seg <- m.seg[,.SD[1], by=.(sample, cell, chrom, start, end)]
  
  return(m.seg)
}
