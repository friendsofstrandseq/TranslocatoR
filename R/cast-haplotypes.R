#' Cast the cleaned, phased, segmented data into a new shape, split haplotypes
#' 
#' @param segs matrix of segments with their class
#' 
#' @return a matrix where the rows are cells and the columns are recurrent segments for each haplotype
#' @author Alex van Vliet
#' 
#' @export

cast.haplotypes <- function(segs) {
  
  # the label is the unique segment
  segs[, label:=paste0(chrom, ":", start, "-", end)]
  # split and factorize haplotypes
  segs[,c("H1","H2"):=.(substring(class, 0,1), substring(class,2))]
  segs[,c("H1_factor","H2_factor"):=.(as.numeric(factor(H1)),as.numeric(factor(H2)))]
  
  # make two matrices for each individual haplotype, then bind them together after renaming
  segs.castH1 <- dcast(segs, cell~label,value.var = 'H1_factor', drop=F)
  setcolorder(segs.castH1, mixedorder(names(segs.castH1)))
  segs.castH2 <- dcast(segs, cell~label,value.var = 'H2_factor', drop=F)
  setcolorder(segs.castH2, mixedorder(names(segs.castH2)))
  
  setnames(segs.castH1, names(segs.castH1)[-1], paste0(names(segs.castH1)[-1], ".H1"))
  setnames(segs.castH2, names(segs.castH2)[-1], paste0(names(segs.castH2)[-1], ".H2"))
  segs.H1H2 <- cbind(segs.castH1, segs.castH2[,-1])
  
  return(segs.H1H2)
}

## This is a similar function for pq analysis
cast.haplo.small <- function(final, haplo){
  tmp <- dcast(final, cell~chrom,value.var = paste0('H',haplo,'_factor'), drop=F)
  setcolorder(tmp, mixedsort(names(tmp)))
  return(tmp)
}
