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
library(assertthat)

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
  
  if (!is.null(regions)) {
    use.regions <- TRUE
    assert_that("chrom" %in% colnames(regions))
    assert_that("sample" %in% colnames(regions))
    assert_that("start" %in% colnames(regions))
    assert_that("end" %in% colnames(regions))
  }
  
  data.folder <- tools::file_path_as_absolute(data.folder)
  if(!dir.exists(data.folder)){
    stop("Please check that you provided the correct path for the data folder. Stopping execution.")
  }
  else{
    cat(paste("Data folder is", data.folder))
  }
  
  for (cursample in samples) {
    
    output.folder <- file.path(output.folder, cursample)
    tryCatch({
      if(!dir.exists(file.path(output.folder))) {
        cat("Creating sample output folder")
        dir.create(file.path(output.folder), recursive=T) }
      else {
        stop("The folder for this sample already exists. Stopping execution to prevent overwriting.") }
          }, warning = function(w) {print("Could not create output folder.")
      })
    
    pval.output <- file.path(output.folder, "pvalues")
    tryCatch({
      dir.create(file.path(output.folder, "pvalues"))
      }, warning = function(w) {
        stop("Could not create output folder for p-values.")
    })
    
    counts.folder <- file.path(data.folder, "counts", cursample)
    if(!dir.exists(counts.folder)) {
      stop(paste("The folder with read counts", counts.folder, "for sample", cursample, "is not in the data folder. Stopping execution."))
    }
    
    phased.folder <- file.path(data.folder, "strand_states", cursample)
    if(!dir.exists(phased.folder)) {
      stop(paste("The folder with phased haplotypes", phased.folder, "for sample", cursample, "is not in the data folder. Stopping execution."))
    }
    
    segments.folder <- file.path(data.folder, "segmentation", cursample)
    if(!dir.exists(segments.folder)) {
      stop(paste("The folder with segments", segments.folder, "for sample", cursample, "is not in the data folder. Stopping execution."))
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
    ## Analysis: pq + suspected translocation file
    #####################
    
    if (options == "pq") {
      
      phased.pq <- gethaplopq(phased)
      
      if (translocations == T) {
        tryCatch({
          phased.pq <- merge(phased.pq, get(paste0(sample.name),".tr")[, -"sample"], by = "cell", all.x=T)
        }, warning = function(w){
          stop("Check that the 'sample' column of your translocations file matches the 'sample' column used by MosaiCatcher.")
        })
      }
      
      # if a 'regions' file is provided, use those to make extra potential translocations
      if (regions == T) {
        regions.cursample <- regions[sample == cursample]
        setkey(regions.cursample, chrom, start, end)
        setkey(phased, chrom, start, end)
        overlaps <- foverlaps(phased, regions.cursample)
        
        overlaps[!is.na(sample), grp := .GRP, by=.(chrom, start, end)]
        overlaps <- overlaps[!is.na(grp)]
        overlaps[, label := paste0(chrom, ":", start, "-", end)]
        overlaps[, startdiff := abs(start - i.start)]
        overlaps[, enddiff := abs(end - i.end)]
        overlaps[, totaldiff := startdiff+enddiff]
        overlaps <- overlaps[, .SD[which.min(totaldiff)], by=.(cell, chrom, start, end)]
        translos <- dcast(overlaps, cell~label, value.var = "class", fill=NA)
        
        merge(phased, translos, by = "cell")
      }
      
      p.values <- get.pvalue.dt(phased.pq)
      
      ###########
      ## Output
      ###########
      
      write.table(p.values, file.path(pval.output, "significance-table.txt"), quote = F, row.names = F)
      write.table(phased.pq, file.path(output.folder, "haplotypes-per-arm.txt"), quote = F, row.names = F)
    
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
      
      #################################
      ## Analysis: p-value calculations
      #################################
      
      p.values <- get.pvalue.dt(final.segs.H1H2)
      
      #############
      ## Output
      #############
    
      write.table(p.values, file.path(pval.output, "significance-table.txt"), quote = F, row.names = F)
      write.table(final.segs.H1H2, file.path(output.folder, "haplotypes-per-segment.txt"), quote = F, row.names = F)
      write.table(recurrent.segs, file.path(output.folder, "recurrent-segments.txt"), quote = F, row.names = F)
    
    # option: segment if statement  
    }
    
    #############
    ### Plotting
    #############
    
    ## sample loop 
  }
  
## end of function  
}  
