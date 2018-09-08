#' Find translocations
#' 
#' Find translocated chromosomal segments by correlating haplotype strand states. Fully vectorized approach heavily reliant on data.table.
#' 
#' @param data.folder MosaiCatcher output folder
#' @param samples on which samples in the MosaiCatcher output folder would you like to run TranslocatoR
#' @param trfile list of textfiles with manually determined translocation states if desired, must contain cell ID and approximate start/end
#' @param output.folder absolute path to output folder for TranslocatoR data
#' @param blacklist whether to use the blacklist for centromeres and short arms for acrocentric chromosomes. defaults to True. use is strongly recommended.
#' @return matrix
#' 
#' @author Alex van Vliet
#' TODO: namespace http://r-pkgs.had.co.nz/namespace.html
#' TODO: what happens when segment is exactly half (not majority/minority)

### Libraries --> move to namespace
library(data.table)
library(gtools)
library(GenomicRanges)
library(rprojroot)
library(dplyr)

translocatoR <- function(data.folder, samples, trfile, output.folder, binsize, blacklist = T) {
  
  ###############################
  ### File checking, initializing
  ###############################
  
  root.dir <- rprojroot::find_root("translocator.R")
  
  if (blacklist == F) {
    cat("Proceeding without blacklist. This will affect performance negatively.")
  }
  else {
    blacklist <- fread(file.path(root.dir, "data/blacklist.txt"))
  }
  
  # default binsize 100kb
  if(missing(binsize)) {
    binsize <- as.integer(100000)
  }
  else {
    binsize <- as.integer(binsize)
  }
  
  data.folder <- tools::file_path_as_absolute(data.folder)
  if(!dir.exists(data.folder)){
    stop("Please check that you provided the correct path for the data folder. Stopping execution.")
  }
  else{
    cat(paste("Data folder is", data.folder))
  }
  
  for (sample in samples) {
    
    output.folder <- file.path(output.folder, sample)
    tryCatch({
      if(!dir.exists(file.path(output.folder))) {
        dir.create(file.path(output.folder)) }
      else {
        stop("The folder for this sample already exists. Stopping execution to prevent overwriting.") }
          }, warning = function(w) {print("Could not create output folder.")
      })
    
    counts.folder <- file.path(data.folder, "counts", sample)
    if(!dir.exists(counts.folder)) {
      stop(paste("The folder with read counts", counts.folder, "for sample", sample, "is not in the data folder. Stopping execution."))
    }
    
    phased.folder <- file.path(data.folder, "strand_states", sample)
    if(!dir.exists(phased.folder)) {
      stop(paste("The folder with phased haplotypes", phased.folder, "for sample", sample, "is not in the data folder. Stopping execution."))
    }
    
    segments.folder <- file.path(data.folder, "segmentation", sample)
    if(!dir.exists(segments.folder)) {
      stop(paste("The folder with segments", segments.folder, "for sample", sample, "is not in the data folder. Stopping execution."))
    }
    
    # read in files. uses normalized counts and assumes the exact output structure and names of MosaiCatcher
    tryCatch({
      counts <- fread(paste("zcat", file.path(counts.folder, paste0(binsize, "_fixed_norm.txt.gz"))))
    }, warning = function(w) {
      stop(paste("Can't unzip and/or read the counts file", paste0(binsize, "_fixed_norm.txt.gz")))
    })
    
    tryCatch({
      phased <- fread(file.path(phased.folder, "final.txt"))
    }, warning = function(w) {
      stop(paste("Can't read the phased counts file final.txt in", phased.folder))
    })
    
    tryCatch({
      segments <- fread(file.path(segments.folder, paste0(binsize, "_fixed_norm.txt")))
    }, warning = function(w) {
      stop(paste("Can't read the segments file", paste0(binsize, "_fixed_norm.txt"), "in", segments.folder))
    })
    
    #################
    ## Analysis start
    #################
    
    total.cells <- counts[, length(unique(cell))]
    cat(paste("Analyzing", total.cells, "cells in sample", sample))
    
    ### Combining phased data and counts data to get phased counts (flipping WC to CW where necessary)
    
    # start by blacklisting regions and then removing leftover None bins
    setkey(blacklist, chrom, start, end)
    counts.blacklist <- foverlaps(counts, blacklist)
    
    # all bins that aren't in the blacklist will have a NA value in the (blacklist) start column so just select for those
    cat(paste("Removing", nrow(counts.blacklist[!is.na(start)]), "blacklisted bins in", total.cells, "cells"))
    counts <- counts.blacklist[is.na(start)]
    cat(paste("Removing", nrow(counts.blacklist[class=='None']), "None bins in", total.cells, "cells"))
    counts <- counts[class!='None']
    
    # remove superfluous columns, rename
    counts[, c("start", "end") := NULL]
    setnames(counts, c("i.start", "i.end"), c("start", "end"))
    
    # identify stretchs of consecutive states (taken from Venla's SCE detector) -- remove here?
    counts <- counts[order(sample, cell, chrom, start, end)]
    counts$cnsc <- cumsum(counts[, .(consecutive = c(1,abs(diff(as.numeric(factor(class)))))), by = .(sample, cell, chrom)]$consecutive)
    
    # FOR DEBUG see all state changes
    # counts[,.SD[c(1,.N)],by=cnsc]
    
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
    switch.counts[,prop.table(table(class))]
    
    counts <- switch.counts
    
    # get segments with lowest sse
    segments.many <- segments[,.SD[sse==min(sse)],by=chrom][,.(chrom, start, end)]
    
    ## DEBUG recurring segments / discrepancy with segmented region --> chr13
    counts <- counts[order(sample, cell, chrom, start, end)]
    counts$cnsc <- cumsum(counts[, .(consecutive = c(1,abs(diff(as.numeric(factor(class)))))), by = .(sample, cell, chrom)]$consecutive)
    counts.recseg <- counts[, .SD[c(1,.N)], by=cnsc]
    counts.recseg[, c("start","end"):=.(min(start), max(end)), by=cnsc]
    counts.recseg <- counts.recseg[,.SD[1], by=.(sample, cell, chrom, start, end)]
    # recurrent breakpoints that cover at least 1Mb
    bps.recseg <- counts.recseg[,.N, by=.(chrom, start, end)][mixedorder(chrom)]
    bps.recseg[, whole.chrom := start == min(start) & end == max(end), by=chrom]
    bps.recseg[, len := end - start]
    bps.recseg <- bps.recseg[whole.chrom==F & N > 1 & len >= 1000000]
    setkey(bps.recseg, chrom, start, end)
    # by how many segments is each of these regions represented
    whichseg <- foverlaps(segments.many, bps.recseg)
    whichseg[, seg.id := paste0(i.start, "-", i.end)]
    howmanyseg <- na.omit(whichseg[,.SD[,uniqueN(seg.id)],by=.(chrom, start, end)])
    setnames(howmanyseg, "V1", "numseg")
    # if the region falls in just one segment and is not the majority it can never be represented
    # so a final check is to see if that one segment covers the majority of the region
    # if not, make a new segment
    newseg <- howmanyseg[numseg == 1]
    newseg <- merge(newseg, whichseg, by=c("chrom", "start", "end"))
    newseg[, seg.len := i.end - i.start]
    newseg[, majority := len/seg.len >= 0.5]
    newseg <- newseg[majority==F, .(start = min(start), end = max(end)),by=.(chrom, i.start, i.end)]
    newseg <- newseg[, cbind(newseg[, .(start = c(i.start, start, end))], newseg[, .(end = c(start, end, i.end))]), by= chrom]
    # now overlap: remove old segment(s) and replace new at the same time
    setkey(newseg, chrom, start, end)
    segments.many <- foverlaps(segments.many, newseg)
    segments.many[start >= i.start & end <= i.end, c("i.start", "i.end"):=.(start, end)]
    segments.many <- segments.many[,.(chrom, i.start, i.end)]
    setnames(segments.many, c("i.start", "i.end"), c("start", "end"))
    
    setkey(segments.many, chrom, start, end)
    setkey(counts, chrom, start, end)
    segment.counts <-foverlaps(counts, segments.many)
    segment.counts <- segment.counts[order(cell,chrom, start, end)]
    
    maj.segments <- segment.counts %>% group_by(sample, cell, chrom, start, end, class) %>% summarise(counter = n()) %>% filter(row_number()==1)
  ## sample loop 
  }
  
## end of function  
}  