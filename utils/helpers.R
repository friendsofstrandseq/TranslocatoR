### mathematical TranslocatoR helpers

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
getcor <- function(casted.haplo){
  spearman.cor <- cor(casted.haplo[, !"cell"], use="pairwise.complete.obs", method = "spearman")
  melt.spearman <- na.omit(melt(get_lower_tri(spearman.cor)))
  spearman.dt <- sapply(melt.spearman[,c(1,2)], function(x) as.character(x))
  spearman.dt <- as.data.table(cbind(spearman.dt, melt.spearman$value))
  setnames(spearman.dt, "V3", "cor")
  return(spearman.dt)
}

# TODO?
p.helper.apply <- function(melt.cormat, castH1H2, x) {
  total <- nrow(melt.cormat)
  
  p.helper(final.segs.H1H2[,.(get(melt.cormat[["Var1"]][x]), get(melt.cormat[["Var2"]][x]))])
}