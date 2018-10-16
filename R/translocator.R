#' Find translocations
#' 
#' Find translocated chromosomal segments by correlating haplotype strand states in Mosaicatcher output. 
#' Fully vectorized approach heavily reliant on data.table.
#' 
#' @import data.table
#' @import gtools
#' @import ggplot2
#' @import stringr
#' @import discreteMTP
#' @import assertthat
#'
#' @param data.folder path to MosaiCatcher data folder
#' @param output.folder absolute path to output folder for TranslocatoR data
#' @param samples on which samples in the MosaiCatcher output folder would you like to run TranslocatoR
#' @param options can take values "segments", "pq", "majority". Defaults to "segments".
#' @param binsize which binsize to use, defaults to 100kb
#' @param cutoff cutoff for significant FDR-corrected p-values, defaults to 0.01
#' @param trfile path to txt file with a suspected translocation. 
#' Can contain multiple translocations for one sample in separate columns, 
#' can also be a list of files for different samples. 
#' The first two columns must contain the cell IDs and the sample name.
#' @param regions list of regions in the format "chr#:<start>-<end>" for potential translocations
#' @param blacklist whether to use the blacklist for centromeres and short arms for acrocentric chromosomes. defaults to True. use is strongly recommended.
#' 
#' @examples
#' translocatoR(data.folder = "/g/korbel2/StrandSeq/20180727_MosiacatcherResults/all-2018-08-02", output.folder = "/g/korbel/vliet/translocations", samples = c("RPE-BM510", "C7-data"))
#' @return matrix
#' 
#' @author Alex van Vliet
#' 
#' @export

# TODO: namespace http://r-pkgs.had.co.nz/namespace.html
# TODO: what happens when segment is exactly half (not majority/minority)
# TODO: option to look at one or more segments (instead of seeing which segments are characterized by strand changes, get a segment and look which state it has in each chrom)
# TODO: return list of interesting bps/potential inversions? -> selecting from bps
# TODO: samples="all" to do all in the MC dir?
# 


translocatoR <- function(data.folder, output.folder, samples, options = "segments",
                         binsize = 100000L, cutoff = 0.01, regions = NULL, trfile = NULL, blacklist = T) {
  
  ###############################
  ### File checking, initializing
  ###############################
  
  # locate data directory
  translocator.data <- system.file(package = "TranslocatoR", "data", mustWork = T)
  
  if (blacklist == F) {
    cat("Proceeding without blacklist. This will affect performance negatively.\n")
  }
  else {
    blacklist <- fread(file.path(translocator.data, "blacklist.txt"))
  }
  
  # default binsize 100kb, integer to prevent scientific notation
  if(class(binsize) != "integer") {
    binsize <- as.integer(binsize)
  }
  
  translocations <- FALSE
  
  if (!is.null(trfile)) {
    translocations <- TRUE
    for (f in trfile) {
      tryCatch({
        tr.tmp <- fread(file.path(f))
      }, warning = function(w) {
        stop(paste("Can't read the translocation file", f, "you provided.\n"))
      })
      # TODO: assertthat columns
      if (tr.tmp[, length(unique(sample))] != 1) {
        stop("Your translocation file contains more than one sample name in the sample column.\n")
      }
      sample.name <- tr.tmp[, unique(sample)]
      # name the translocation dt after the sample
      assign(paste0(sample.name, ".tr"), tr.tmp)
    }
  }
  
  use.regions <- FALSE
  
  if (!is.null(regions)) {
    use.regions <- TRUE
    tryCatch({
      regions.f <- fread(file.path(regions))
    }, warning = function(w) {
      stop(paste("Can't read the regions file", f, "you provided.\n"))
    })
    assert_that("chrom" %in% colnames(regions.f))
    assert_that("sample" %in% colnames(regions.f))
    assert_that("start" %in% colnames(regions.f))
    assert_that("end" %in% colnames(regions.f))
  }
  
  data.folder <- tools::file_path_as_absolute(data.folder)
  if(!dir.exists(data.folder)){
    stop("Please check that you provided the correct path for the data folder. Stopping execution.\n")
  }
  else{
    cat(paste("Data folder is", data.folder, "\n"))
  }
  
  for (cursample in samples) {
    
    cat(paste("Current sample is", cursample, "\n"))
    output.folder.cursamp <- file.path(output.folder, cursample)
    tryCatch({
      if(!dir.exists(file.path(output.folder.cursamp))) {
        cat("Creating sample output folder\n")
        dir.create(file.path(output.folder.cursamp), recursive=T) }
    }, warning = function(w) {print("Could not create output folder.")
    })
    
    pval.output <- file.path(output.folder.cursamp, "pvalues")
    tryCatch({
      if (!dir.exists(file.path(output.folder.cursamp, "pvalues"))) {
        dir.create(file.path(output.folder.cursamp, "pvalues")) }
      }, warning = function(w) {
        stop("Could not create output folder for p-values.\n")
      })
    
    phased.folder <- file.path(data.folder, "strand_states", cursample)
    if(!dir.exists(phased.folder)) {
      stop(paste("The folder with phased haplotypes", phased.folder, "for sample", cursample, "is not in the data folder. Stopping execution.\n"))
    }
  
    # read in files. uses normalized counts and assumes the exact output structure and directory names of MosaiCatcher
    cat("Reading files\n")
    
    tryCatch({
      phased <- fread(file.path(phased.folder, "final.txt"))
    }, warning = function(w) {
      stop(paste("Can't read the phased counts file final.txt in", phased.folder, "\n"))
    })
    
    #####################
    ## Analysis: pq + suspected translocation file
    #####################
    
    if (options == "pq") {
      
      cat("Starting analysis\n")
      phased.pq <- gethaplopq(phased)
      
      if (translocations == T) {
        
        tryCatch({
          phased.pq <- merge(phased.pq, get(paste0(cursample,".tr"))[, -"sample"], by = "cell", all.x=T)
        }, warning = function(w){
          stop("Check that the 'sample' column of your translocations file matches the 'sample' column used by MosaiCatcher.\n")
        })
      }
      
      # if a 'regions' file is provided, use those to make extra potential translocations
      if (use.regions == T) {
        
        regions.cursample <- regions.f[sample == cursample]
        translos <- region.to.tr(regions.cursample, phased)
        
        # convert translocation file to factors
        names <- colnames(translos[, -"cell"])
        translos[, paste0(names, ".H1") := lapply(.SD, function(x) {as.numeric(factor(substring(x, 0, 1)))}), .SDcols = names]
        translos[, paste0(names, ".H2") := lapply(.SD, function(x) {as.numeric(factor(substring(x, 2)))}), .SDcols = names]
        
        phased.pq <- merge(phased.pq, translos[, .SD, .SDcols = names(translos) %like% ".H[1-2]" | names(translos) %like% "cell"], by = "cell")
      }
      
      p.values <- get.pvalue.dt(phased.pq)
      potential.tr <- p.values[str_extract(segA, "chr[0-9X-Y]") != str_extract(segB, "chr[0-9X-Y]")][pBH < cutoff]
      
      ###########
      ## Output
      ###########
      
      cat("Writing output files\n")
      write.table(p.values, file.path(pval.output, "significance-table.txt"), quote = F, row.names = F)
      write.table(potential.tr, file.path(pval.output, "co-segregations.txt"), quote = F, row.names = F)
      write.table(phased.pq, file.path(output.folder.cursamp, "haplotypes-per-arm.txt"), quote = F, row.names = F)
    
    } # option = pq if statement
    
    if (options == "segments") {
      
      ###############
      ## File checking and reading
      ###############
      
      counts.folder <- file.path(data.folder, "counts", cursample)
      if(!dir.exists(counts.folder)) {
        stop(paste("The folder with read counts", counts.folder, "for sample", cursample, "is not in the data folder. Stopping execution.\n"))
      }
      
      segments.folder <- file.path(data.folder, "segmentation", cursample)
      if(!dir.exists(segments.folder)) {
        stop(paste("The folder with segments", segments.folder, "for sample", cursample, "is not in the data folder. Stopping execution.\n"))
      }
      
      tryCatch({
        counts <- fread(paste("zcat", file.path(counts.folder, paste0(binsize, "_fixed_norm.txt.gz"))))
      }, warning = function(w) {
        stop(paste("Can't unzip and/or read the counts file", paste0(binsize, "_fixed_norm.txt.gz\n")))
      })
      
      tryCatch({
        segments <- fread(file.path(segments.folder, paste0(binsize, "_fixed_norm.txt")))
      }, warning = function(w) {
        stop(paste("Can't read the segments file", paste0(binsize, "_fixed_norm.txt"), "in", segments.folder, "\n"))
      })
      
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
