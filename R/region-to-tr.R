#' From a list of regions get the corresponding strand state for that region in each cell
#' 
#' @param regions all regions for the current sample
#' @param phased the phased MosaiCatcher data
#' 
#' @return a matrix with as many columns as regions. Each row is a cell, each column is the strand state for each cell in that region.
#' @author Alex van Vliet
#' 
#' @export

region.to.tr <- function(regions, phased) {
  
  # overlap regions with phased data to get the state in each region for each chromosome
  setkey(regions, chrom, start, end)
  setkey(phased, chrom, start, end)
  overlaps <- foverlaps(phased, regions)
  
  # the 'sample' column will be NA for all rows without overlap
  overlaps <- overlaps[!is.na(sample)]
  # label for identification of region later
  overlaps[, label := paste0(chrom, ":", start, "-", end)]
  # if multiple phased segments overlap a region, choose the one with the most overlap
  overlaps[, startdiff := abs(start - i.start)]
  overlaps[, enddiff := abs(end - i.end)]
  overlaps[, totaldiff := startdiff + enddiff]
  overlaps <- overlaps[, .SD[which.min(totaldiff)], by=.(cell, chrom, start, end)]
  translos <- dcast(overlaps, cell~label, value.var = "class", fill=NA)
  
  return(translos)
}
