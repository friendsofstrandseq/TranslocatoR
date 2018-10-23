#' Get all cells that do not conform to a translocation pattern to check for SVs/SCEs
#' 
#' @param phased the StrandphaseR file with all phased segments
#' @param hapmatrix the matrix that contains the split haplotypes for each chromosome/segment 
#' @param translocations the file with all potential translocations
#' @param x the rownumber of the translocations file
#' 
#' @author Alex van Vliet
#' 
#' @export

get.outliers <- function(phased, hapmatrix, translocations, x) {
  
  partnerA <- translocations[["segA"]][x]
  partnerB <- translocations[["segB"]][x]
  posneg <- translocations[["cor"]][x]
  
  # if the correlation is positive, outliers are cells where states are not equal and vice versa
  if (posneg > 0) {
    outlier.cells <- hapmatrix[,.(cell, A = get(partnerA), B = get(partnerB))][A != B, cell]
    chroms <- c(str_extract(partnerA, "chr[0-9X-Y]+"), str_extract(partnerB, "chr[0-9X-Y]+"))
  }
  else {
    outlier.cells <- hapmatrix[,.(cell, A = get(partnerA), B = get(partnerB))][A == B, cell]
    chroms <- c(str_extract(partnerA, "chr[0-9X-Y]+"), str_extract(partnerB, "chr[0-9X-Y]+"))
  }
  
  outliers <- phased[cell %in% outlier.cells & chrom %in% chroms][order(cell)]
  
  return(outliers)
}
