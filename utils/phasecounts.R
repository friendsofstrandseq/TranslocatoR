#' Phase MosaiCatcher count data using StrandPhaseR output
#' 
#' @param counts clean count matrix
#' @param phased.segments StrandPhaseR output
#' 
#' @return count matrix with WC -> CW as computed by StrandPhaseR
#' @author Alex van Vliet
#' 

phase.counts <- function(counts, phased.segments) {
  
  # overlap the phased data with the count data and switch all WC to CW where necessary
  setkey(counts, sample, cell, chrom, start, end)
  setkey(phased, sample, cell, chrom, start, end)
  switch.counts <- foverlaps(counts, phased, mult="last")
  switch.counts[i.class == 'WC' & class == "CW", i.class := 'CW']
  switch.counts <- switch.counts[,.(sample, cell, chrom, i.start, i.end, c, w, i.class)]
  setnames(switch.counts, c("i.start", "i.end", "i.class"), c("start", "end", "class"))
  
  cat("The strand inheritance before phasing:\n")
  print(counts[,prop.table(table(class))])
  cat("The strand inheritance after phasing:\n")
  print(switch.counts[,prop.table(table(class))])
  
  return(switch.counts)
}
