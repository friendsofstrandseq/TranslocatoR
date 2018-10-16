#' Split haplotypes per chromosome arm
#' 
#' @param phased MosaiCatcher phased segments file
#' 
#' @return data.table with all haplotype/arm combinations
#' 
#' @author Alex van Vliet
#' 
#' @export

gethaplopq <- function(phased) {
  
  # uses MosaiCatcher phased by default
  # if you want to make your own "phased", run TranslocatoR with options = "segments" first
  phased[, c('H1','H2') := .(substring(class, 0, 1), substring(class,2))]
  phased[, c('H1_factor','H2_factor') := .(as.numeric(factor(H1)), as.numeric(factor(H2)))]
  cursamp.p <- phased[, .SD[1], by=.(sample, cell, chrom)]
  cursamp.q <- phased[, .SD[.N], by=.(sample, cell, chrom)]
  
  # get all of the haplotypes per arm
  cursamp.p.H1 <- cast.haplo.small(cursamp.p, 1)
  cursamp.q.H1 <- cast.haplo.small(cursamp.q, 1)
  cursamp.q.H2 <- cast.haplo.small(cursamp.q, 2)
  cursamp.p.H2 <- cast.haplo.small(cursamp.p, 2)
  
  # rename everything so that haplotypes and arms are distinguishable
  setnames(cursamp.p.H1, grep("chr", names(cursamp.p.H1), value=T), paste0(grep("chr", names(cursamp.p.H1), value=T), ".p.H1"))
  setnames(cursamp.p.H2, grep("chr", names(cursamp.p.H2), value=T), paste0(grep("chr", names(cursamp.p.H2), value=T), ".p.H2"))
  setnames(cursamp.q.H1, grep("chr", names(cursamp.q.H1), value=T), paste0(grep("chr", names(cursamp.q.H1), value=T), ".q.H1"))
  setnames(cursamp.q.H2, grep("chr", names(cursamp.q.H2), value=T), paste0(grep("chr", names(cursamp.q.H2), value=T), ".q.H2"))
  
  # make one big data.table out of the four haplotype/arm combinations
  cursamp.pqH1H2 <- merge(cursamp.p.H1, cursamp.p.H2, all = T)
  cursamp.pqH1H2 <- merge(cursamp.pqH1H2, cursamp.q.H1, all = T)
  cursamp.pqH1H2 <- merge(cursamp.pqH1H2, cursamp.q.H2, all = T)
  
  return(cursamp.pqH1H2)
}
