#' Find translocations
#' 
#' Find translocated chromosomal segments by correlating haplotype strand states in Mosaicatcher output. 
#' Fully vectorized approach heavily reliant on data.table.
#' 
#' @param data.folder path to MosaiCatcher data folder
#' @param output.folder absolute path to output folder for TranslocatoR data
#' @param samples on which samples in the MosaiCatcher output folder would you like to run TranslocatoR
#' @param options can take values "segments", "pq", "majority". Defaults to "segments".
#' @param binsize which binsize to use, defaults to 100kb
#' @param trfile path to txt file with a suspected translocation. 
#' Can contain multiple translocations for one sample in separate columns, 
#' can also be a list of files for different samples. 
#' The first two columns must contain the cell IDs and the sample name.
#' @param regions list of regions in the format "chr#:<start>-<end>" for potential translocations
#' @param blacklist whether to use the blacklist for centromeres and short arms for acrocentric chromosomes. defaults to True. use is strongly recommended.
#' 
#' @example translocatoR(data.folder = "/g/korbel2/StrandSeq/20180727_MosiacatcherResults/all-2018-08-02", 
#' output.folder = "/g/korbel/vliet/translocations", samples = c("RPE-BM510", "C7-data"))
#' @return matrix
#' 
#' @author Alex van Vliet
#' TODO: namespace http://r-pkgs.had.co.nz/namespace.html
#' TODO: what happens when segment is exactly half (not majority/minority)
#' TODO: option to look at one or more segments (instead of seeing which segments are characterized by strand changes, get a segment and look which state it has in each chrom)
#' TODO: return list of interesting bps/potential inversions? -> selecting from bps
#' TODO: samples="all" to do all in the MC dir?

### Libraries --> move to namespace
library(data.table)
library(gtools)
library(GenomicRanges) # necessary?
library(rprojroot)
library(dplyr) # necessary?
library(ggplot2)
library(stringr)
library(discreteMTP)
library(R.utils) # for sourcing all helpers

translocatoR <- function(data.folder, output.folder, samples, options = "segments",
                         binsize = 100000L, regions = NULL, trfile = NULL, blacklist = T) {
  
  ###############################
  ### File checking, initializing
  ###############################
  
  # locate root directory
  root.dir <- rprojroot::find_root("translocator.R")
  # load helper functions
  sourceDirectory(file.path(root.dir, "utils"))
  
  if (blacklist == F) {
    cat("Proceeding without blacklist. This will affect performance negatively.")
  }
  else {
    blacklist <- fread(file.path(root.dir, "data/blacklist.txt"))
  }
  
  # default binsize 100kb, integer to prevent scientific notation
  if(class(binsize) != "integer") {
    binsize <- as.integer(binsize)
  }
  
  if (!is.null(trfile)) {
    translocations <- TRUE
    for (f in trfile) {
      tryCatch({
        tr.tmp <- fread(file.path(f))
      }, warning = function(w) {
        stop(paste("Can't read the translocation file", f, "you provided."))
      })
      # TODO: assertthat columns
      if (tr.tmp[, length(unique(sample))] != 1) {
        stop("Your translocation file contains more than one sample name in the sample column.")
      }
      sample.name <- tr.tmp[, unique(sample)]
      # name the translocation dt after the sample
      assign(paste0(sample.name, ".tr"), tr.tmp)
    }
  }
  ### solution for reading trfile?
  getroworder <- function(states, tr){
    roworder <- match(states[,cell], get(tr)[,cell])
    return(roworder)
  }
  cast.haplos <- function(states, trstates, haplo){
    tmp <- dcast(states, cell~chrom,value.var = paste0('H',haplo,'_factor'), drop=F)
    for (x in trstates) {
      tmp <- cbind(tmp, get(x)[getroworder(tmp, x), tr_factor])
      setnames(tmp, "V2", x)
    }
    setcolorder(tmp, mixedsort(colnames(tmp)))
    return(tmp)
  }
  
  cast.haplos <- function(states, trstates, haplo){
    tmp <- dcast(states, cell~chrom,value.var = paste0('H',haplo,'_factor'), drop=F)
    tmp <- tmp[mixedorder(cell)]
    setcolorder(tmp, mixedsort(colnames(tmp)))
    # first column is always cell
    whichtr <- names(trstates)
    trstates <- trstates[mixedorder(trstates)]
    for (x in 2:ncol(trstates)) {
      tmp <- merge(tmp, trstates[, .(cell, get(whichtr[x]))], by="cell", all.x=T)
      setnames(tmp, "V2", whichtr[x])
    }
    return(tmp)
  }
  ### end solution ###
  
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
        cat("Creating sample output folder")
        dir.create(file.path(output.folder), recursive=T) }
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
    cat("Reading files")
    ## TODO: if ends in .gz, zcat otherwise read txt
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
    
    #####################
    ## Start analysis
    #####################
    
    if (options == "pq") {
      
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
      
      if (translocations == T) {
        cursamp.pqH1H2 <- merge(cursamp.pqH1H2, get(paste0(sample.name),".tr")[, -"sample"], by = "cell", all.x=T)
      }
      
      cormat <- suppressWarnings(getcor(cursamp.pqH1H2))
      pvals <- sapply(1:nrow(cormat), p.helper.apply, cormat = cormat, casthaplos = TALL03pqH1H2)
      pvals.dt <- as.data.table(matrix(pvals, nrow = nrow(cormat), ncol = 3, byrow = T))
      setnames(pvals.dt, c("V1", "V2", "V3"), c("p", "x", "n"))
      
      # apply FDR correction (Benjamini-Hochberg) for discrete p-values following a binomial distribution
      cat("FDR-adjusting...")
      pCDFlist <- lapply(1:nrow(pvals.dt), function(i){pbinom(0:pvals.dt[i, x], pvals.dt[i, n], 0.5)})
      pBH <- p.discrete.adjust(pvals.dt[,p], pCDFlist, method = "BH")
      pvals <- data.table(cormat, pvals.dt, pBH)
      
      ###########
      ## Output
      ###########
    }
    if (options == "segments") {
      #################################
      ## Analysis: cleaning and phasing
      #################################
      
      total.cells <- counts[, length(unique(cell))]
      cat(paste("Analyzing", total.cells, "cells in sample", sample))
      
      # remove all None bins and centromeric bins
      counts <- clean.bins(counts, blacklist)
      
      # combine phased data and counts data to get phased counts (flipping WC to CW where necessary)
      counts <- phase.counts(counts, phased)
      
      ###########################################################
      ### Analysis: use MosaiCatcher segments to segment the data 
      ### and get the haplotype-aware state for each segment
      ###########################################################
      
      # find recurrent segments in the count data without looking at the segmented data
      seg.count <- seg.finder(counts)
      recurrent.segs <- rec.seg(seg.count)
      
      # get segments with lowest sse
      segments.many <- segments[,.SD[sse == min(sse)], by = chrom][,.(chrom, start, end)]
      
      # the following line is a segmentation security check: if there are any recurrent segments that are not represented by the 
      # current maximally segmented data, add this segment to the list
      # this step can likely be removed after segmentation in MosaiCatcher improves
      segments.many <- segment.check(segments.many, recurrent.segs)
      
      # combine the count data and the segments to get the phased state per segment
      setkey(segments.many, chrom, start, end)
      setkey(counts, chrom, start, end)
      segment.counts <- foverlaps(counts, segments.many)
      segment.counts <- segment.counts[order(cell,chrom, start, end)]
      
      #######################################################
      ## Analysis: compare all segments to each other to find 
      ## co-segregating segments
      #######################################################
      
      ## print how many cells had 2 or more choices? 
      # segment.counts[, .(nbins = max(.N)), by = .(sample, cell, chrom, start, end, class)][, .N, by=.(sample, cell, chrom, start, end)][N>1]
      
      # majority class per segment from count data
      segment.class <- segment.counts[, .(nbins = max(.N)), by = .(sample, cell, chrom, start, end, class)][, .(class = class[which.max(nbins)]), by=.(sample, cell, chrom, start, end)]
      
      # summarize the majority class per segment by retrieving start and end for each segment
      bp.segs <- seg.finder(segment.class)
      
      # assign each cell + chromosome combination a label of breakpoints
      unique.bps <- bp.segs[,unique(c(start, end)) ,by=.(cell, chrom)][, .(unique.bp = sapply(.SD, toString)), by=.(cell, chrom)]
      bp.segs <- merge(bp.segs, unique.bps)
      
      # by using the breakpoint label unique.bp we can tell how many cells have this exact breakpoint pattern
      bp.segs[, cells.unique.bp := length(unique(cell)), by=.(chrom, unique.bp)]
      
      # get bp.count: how often this exact breakpoint (start-end) occurs in the sample
      bp.segs <- bp.segs[, bp.count := .N, by=.(chrom, start, end)][mixedorder(chrom)]
      bp.segs <- bp.segs[order(cell, chrom, start, end)]
      bp.segs[, len := end - start]
      
      # for breakpoints that occur <=2 times: ignore (melt into surrounding segments) this breakpoint as it's probably an SCE and will not be useful for analysis later
      # THIS IS BUGGY. DOES NOT WORK IF THERE IS MORE THAN ONE RARE BP!!!
      bp.segs[bp.count <= 2, c("start","end", "majority"):=.(min(start), max(end), .SD[which.max(len), class]), by=.(cell, chrom)]
      final.segs <- unique(bp.segs[is.na(majority) | class==majority, .(sample, cell, chrom, start, end, class)])
      
      final.segs.H1H2 <- cast.haplotypes(final.segs)
      
      cormat <- suppressWarnings(getcor(final.segs.H1H2)) # TODO: silence only standard deviation warnings
      
      #################################
      ## Analysis: p-value calculations
      #################################
      
      # for every valid combination of segments (i.e. those that have shared informative cells) get the p-value
      # by applying Fisher's test to a Watson/Crick contingency table
      cat("Starting p-value calculations...")
      pvals <- sapply(1:nrow(cormat), p.helper.apply, cormat = cormat, casthaplos = final.segs.H1H2)
      pvals.dt <- as.data.table(matrix(pvals, nrow = nrow(cormat), ncol = 3, byrow = T))
      setnames(pvals.dt, c("V1", "V2", "V3"), c("p", "x", "n"))
      
      # apply FDR correction (Benjamini-Hochberg) for discrete p-values following a binomial distribution
      cat("FDR-adjusting...")
      pCDFlist <- lapply(1:nrow(pvals.dt), function(i){pbinom(0:pvals.dt[i, x], pvals.dt[i, n], 0.5)})
      pBH <- p.discrete.adjust(pvals.dt[,p], pCDFlist, method = "BH")
      pvals <- data.table(cormat, pvals.dt, pBH)
      
      pvals[, log10.p := -log10(p)]
      pvals[, B.p := p.adjust(p, "bonferroni")]
      
      #############
      ## Output
      #############
    
      write.table()
    
    # option: segment if statement  
    }
    
    #############
    ### Plotting
    #############
    
    ## plot by haplotype (H1-H1, H1-H2, H2-H1, H2-H2)
    
    ## Interchromosomal translocation candidates
    tr.can <- pvals[pBH < 0.01][str_extract(segA, 'chr[0-9X-Y]+') != str_extract(segB, 'chr[0-9X-Y]+')]
    
    ## Segments in consistent order, with lowest-number chrom first
    ## watch for chrX/Y
    tr.can[, c("chrom1", "chrom2") := .(as.numeric(str_extract(segA, '\\d+')), as.numeric(str_extract(segB, '\\d+')))]
    tr.can[chrom1 > chrom2, c("segB", "segA") := .(segA, segB)]
    tr.can <- tr.can[mixedorder(segA)]
    tr.can[, c("chrom1", "chrom2") := .(as.numeric(str_extract(segA, '\\d+')), as.numeric(str_extract(segB, '\\d+')))]
    
    ## get segments back as numeric
    tr.can[, c("start.chrom1", "end.chrom1", "start.chrom2", "end.chrom2") := 
             .(as.numeric(str_extract(segA, '(?<=:)\\d+')), as.numeric(str_extract(segA, '(?<=-)\\d+')), 
               as.numeric(str_extract(segB, '(?<=:)\\d+')), as.numeric(str_extract(segB, '(?<=-)\\d+')))]
    tr.can[, plottogether := .GRP, by=.(chrom1, chrom2)]
    tr.can[, haplo :=paste0(str_extract(segA, "H[1-2]"), "-",str_extract(segB, "H[1-2]"))]
    tr.can[cor > 0, posneg := "positive"]
    tr.can[cor < 0, posneg := "negative"]
    ## ggplot takes only data frames
    plotdf <- setDF(tr.can)
    
    ## updating defaults is the only way to get the text to change colour
    update_geom_defaults("text", list(colour = "white"))
    for (i in unique(plotdf$plottogether)) {
      curpartners <- paste0(unique(plotdf[plotdf$plottogether == i,"chrom1"]),"-", unique(plotdf[plotdf$plottogether == i,"chrom2"]))
      ggplot(plotdf[plotdf$plottogether == i,]) + geom_tile(aes(segA, segB, fill = posneg))+ 
        geom_text(aes(segA, segB, label=formatC(as.numeric(pBH), format = "e", digits = 2)))+
        scale_fill_manual(values = c("positive" = "darkgreen", "negative" = "darkred"))+
        facet_wrap(~haplo)+
        theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 8, hjust = 1))
      ggsave(filename = paste0(sample, "_co-segregation_", curpartners,".pdf"), plot = last_plot(), path = output.folder)
    }
    ## sample loop 
  }
  
## end of function  
}  
