#' mathematical TranslocatoR helpers
#' 
#' @author Alex van Vliet
#' 
#' @export

# deprecated function
p.helper <- function(x) {
  tmp.table <- table(na.omit(x))
  if (nrow(tmp.table)< 2){
    return(c(NA,NA))
  }
  else if (ncol(tmp.table) < 2){
    return(c(NA,NA))
  }
  else {
    ft <- fisher.test(tmp.table)
    or.p <- c(ft$estimate, ft$p.value)
    return(or.p)
  }
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# get all combinations of segments that will not produce NA and get the direction of the correlation for plotting later
getcor <- function(casted.haplo) {
  cor_ <- cor(casted.haplo[, !"cell"], use="pairwise.complete.obs")
  meltcor <- na.omit(melt(get_lower_tri(cor_)))
  cor.dt <- sapply(meltcor[,c(1,2)], function(x) as.character(x))
  cor.dt <- data.table(cor.dt, meltcor$value, stringsAsFactors = F)
  setnames(cor.dt, c("Var1", "Var2", "V2"), c("segA", "segB", "cor"))
  return(cor.dt)
}

## discrete FDR correction:
## get cumulative distribution function (CDF) of binomial distribution
## for each segment comparison
## x = number of succesful trials (equal states)
## n = number of observations (informative cells)

p.helper.cdf <- function(x) {
  tmp.table <- table(na.omit(x))
  if (nrow(tmp.table)< 2){
    return(c(NA,NA))
  }
  else if (ncol(tmp.table) < 2){
    return(c(NA,NA))
  }
  else {
    x <- sum(tmp.table[c(1, 4)])
    n <- sum(tmp.table[1:4])
    ft <- fisher.test(tmp.table)
    p <- c(ft$p.value, x, n)
    return(p)
  }
}

p.helper.apply <- function(cormat, casthaplos, x) {
  total <- nrow(cormat)
  if (x%%round(total/10) == 0) {
    cat(paste(x, "combinations out of", total,"calculated\n"))}
  p.helper.cdf(casthaplos[, c(cormat[["segA"]][x], cormat[["segB"]][x]), with=F])
}

get.pvalue.dt <- function(factordt) {
  # TODO: suppress only stdev warnings
  cormat <- suppressWarnings(getcor(factordt))
  pvals <- sapply(1:nrow(cormat), p.helper.apply, cormat = cormat, casthaplos = factordt)
  pvals.dt <- as.data.table(matrix(pvals, nrow = nrow(cormat), ncol = 3, byrow = T))
  setnames(pvals.dt, c("V1", "V2", "V3"), c("p", "x", "n"))
  
  # apply FDR correction (Benjamini-Hochberg) for discrete p-values following a binomial distribution
  cat("FDR-adjusting...\n")
  pCDFlist <- lapply(1:nrow(pvals.dt), function(i){pbinom(0:pvals.dt[i, x], pvals.dt[i, n], 0.5)})
  pBH <- p.discrete.adjust(pvals.dt[,p], pCDFlist, method = "BH")
  pvals <- data.table(cormat, pvals.dt, pBH)
  
  return(pvals)
}
