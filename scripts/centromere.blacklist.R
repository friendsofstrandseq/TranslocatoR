#' Generate centromere blacklist from UCSC Table Browser gaps + centromeres text files
#' 
#' @param centro the list of centromeres from the UCSC Table browser (http://genome.ucsc.edu/cgi-bin/hgTables > group : All Tables, table : centromeres)
#' @param gap the list of gaps from the UCSC Table browser (http://genome.ucsc.edu/cgi-bin/hgTables > group : All Tables, table : gaps)
#' 
#' @author Alex van Vliet
#'
#' TODO: autodownload with curl
#' uses devtools::install_github("r-lib/rprojroot")

library(rprojroot)
library(data.table)
library(gtools)

centro.blacklist <- function() {
  
  root.dir <- rprojroot::find_root("translocator.R")
  
  centromeres <- fread(file.path(root.dir, "data/centromeres.txt"))
  gaps <- fread(file.path(root.dir, "data/gaps.txt"))
  
  # remove unnecessary columns
  gaps[,c("#bin", "ix", "n") := NULL]
  centromeres[,c("#bin", "name") := NULL]
  
  # add type column for centromere
  centromeres[, type := "centromere"]
  
  allgaps <- rbind(centromeres, gaps, fill=T)
  allgaps <- allgaps[order(chrom, chromStart, chromEnd)] # not necessary
  # start and end centromere column
  allgaps <- merge(allgaps, allgaps[type=='centromere',.(centro.start=min(chromStart), centro.end=max(chromEnd)),by=chrom])
  
  # bool column for selection: extend blacklist to regions surrounding centromeres by 2Mb (contigs, scaffolds, heterochromatin, short arm in the gaps file)
  allgaps[,include := ((chromStart >= (centro.start - 2000000) & chromStart < centro.end)| 
                         (chromEnd <= (centro.end + 2000000) & chromEnd > centro.end) | 
                         type=='heterochromatin' |
                         type == "short_arm")]

  blacklist <- allgaps[include==T,.(start = min(chromStart), end = max(chromEnd)),by=chrom]
  
  return(blacklist[mixedorder(chrom)])
}
