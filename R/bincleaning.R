#' Bin cleaning: remove blacklisted bins and None bins
#' 
#' @param count.data count matrix for the current sample
#' @param blacklist bin blacklist for centromeres and adjacent hard-to-sequence regions
#' 
#' @return the count matrix with centromeric bins and bins with no reads (None bins) removed
#' @author Alex van Vliet
#' 
#' @export

clean.bins <- function(counts, blacklist) {
  
  # start by blacklisting regions and then removing leftover None bins
  setkey(blacklist, chrom, start, end)
  counts.blacklist <- foverlaps(counts, blacklist)
  
  # all bins that aren't in the blacklist will have a NA value in the (blacklist) start column so just select for those
  counts <- counts.blacklist[is.na(start)]
  counts <- counts[class!='None']
  
  cat(paste("Removed", nrow(counts.blacklist[!is.na(start)]), "blacklisted bins in", total.cells, "cells\n"))
  cat(paste("Removed", nrow(counts.blacklist[class=='None']), "None bins in", total.cells, "cells\n"))
  
  # remove superfluous columns, rename
  counts[, c("start", "end") := NULL]
  setnames(counts, c("i.start", "i.end"), c("start", "end"))
  
  return(counts)
}
