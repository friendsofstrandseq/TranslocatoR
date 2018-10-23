#' Cast the phased majority states into a new shape, split haplotypes
#' 
#' @param segs matrix of sample + cell + chromosomes with their majority class
#' 
#' @return a matrix where the rows are cells and the columns are haplotypes
#' @author Alex van Vliet
#' 
#' @export

cast.majority <- function(segs) {
  
  # split and factorize haplotypes
  segs[,c("H1","H2"):=.(substring(class, 0,1), substring(class,2))]
  segs[,c("H1_factor","H2_factor"):=.(as.numeric(factor(H1)),as.numeric(factor(H2)))]
  
  # make two matrices for each individual haplotype, then bind them together after renaming
  segs.castH1 <- dcast(segs, cell~chrom,value.var = 'H1_factor', drop=F)
  setcolorder(segs.castH1, mixedorder(names(segs.castH1)))
  segs.castH2 <- dcast(segs, cell~chrom,value.var = 'H2_factor', drop=F)
  setcolorder(segs.castH2, mixedorder(names(segs.castH2)))
  
  setnames(segs.castH1, names(segs.castH1)[-1], paste0(names(segs.castH1)[-1], ".H1"))
  setnames(segs.castH2, names(segs.castH2)[-1], paste0(names(segs.castH2)[-1], ".H2"))
  segs.H1H2 <- cbind(segs.castH1, segs.castH2[,-1])
  
  return(segs.H1H2)
}