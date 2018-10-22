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
#' @param options can take values "segments", "pq", "majority". Defaults to "pq".
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

# TODO: what happens when segment is exactly half (not majority/minority)
# TODO: samples="all" to do all in the MC dir?

translocatoR <- function(data.folder, output.folder, samples, options = "pq",
                         binsize = 100000L, cutoff = 0.01, regions = NULL, trfile = NULL, blacklist = T) {
  
  ###############################
  ### File checking, initializing
  ###############################
  
  # default binsize 100kb, integer to prevent scientific notation
  if(class(binsize) != "integer") {
    binsize <- as.integer(binsize)
  }
  
  if (options != "majority" & options != "pq" & options != "segments") {
    stop(cat("Provide a valid option: majority, pq or segments."))
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
      
      assert_that("cell" %in% colnames(tr.tmp))
      assert_that("sample" %in% colnames(tr.tmp))
      
      if (tr.tmp[, length(unique(sample))] != 1) {
        stop("Your translocation file contains more than one sample name in the sample column.\n")
      }
      
      # correct format
      tr.tmp <- cbind(tr.tmp[, c("cell", "sample")], apply(tr.tmp[, -c("cell", "sample")], 2, function(x) {as.numeric(factor(x))}))
      sample.name <- tr.tmp[, unique(sample)]
      
      # name the translocation dt after the sample for later reference
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
    
    outlier.output <- file.path(output.folder.cursamp, "outliers")
    tryCatch({
      if (!dir.exists(file.path(output.folder.cursamp, "outliers"))) {
        dir.create(file.path(output.folder.cursamp, "outliers")) }
      }, warning = function(w) {
        stop("Could not create output folder for outliers.\n")
      })
    
    tr.output <- file.path(output.folder.cursamp, "translocations")
    tryCatch({
      if (!dir.exists(file.path(output.folder.cursamp, "translocations"))) {
        dir.create(file.path(output.folder.cursamp, "translocations")) }
    }, warning = function(w) {
      stop("Could not create output folder for translocations.\n")
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
    cat("\tRead phased file\n")
    
    ####################################################
    ## Analysis: majority + suspected translocation file
    ####################################################
    
    if (options == "majority") {
      
      cat("Starting analysis\n")
      phased[, len := end - start]
      # get the majority state by calculating which state covers the most DNA
      majority <- phased[,.(total.len = sum(len)), by = .(sample, cell, chrom, class)][, .(class = class[which.max(total.len)]), by = .(sample, cell, chrom)]
      
      maj.cast <- cast.majority(majority)
      
      if (translocations == T) {
        
        if (exists(paste0(cursample,".tr"))){
          tryCatch({
            maj.cast <- merge(maj.cast, get(paste0(cursample,".tr"))[, -"sample"], by = "cell", all.x=T)
          }, warning = function(w){
            stop("Check that the 'sample' column of your translocations file matches the 'sample' column used by MosaiCatcher.\n")
          })
        }
      }
      
      ####################################
      ## Analysis: majority + regions file
      ####################################
      
      # if a 'regions' file is provided, use those to make extra potential translocations
      if (use.regions == T) {
        
        regions.cursample <- regions.f[sample == cursample]
        
        # only use the file if regions for the current sample are provided
        if (dim(regions.cursample)[1] != 0){
          
          translos <- region.to.tr(regions.cursample, phased)
          
          # convert translocation file to factors
          names <- colnames(translos[, -"cell"])
          translos[, paste0(names, ".H1") := lapply(.SD, function(x) {as.numeric(factor(substring(x, 0, 1)))}), .SDcols = names]
          translos[, paste0(names, ".H2") := lapply(.SD, function(x) {as.numeric(factor(substring(x, 2)))}), .SDcols = names]
          
          maj.cast <- merge(maj.cast, translos[, .SD, .SDcols = names(translos) %like% ".H[1-2]" | names(translos) %like% "cell"], by = "cell")
        }
      }
      
      recurrent <- rec.seg(phased)
      p.values <- get.pvalue.dt(maj.cast)
      potential.tr <- p.values[str_extract(segA, "chr[0-9X-Y]+") != str_extract(segB, "chr[0-9X-Y]+")][pBH < cutoff][order(pBH)]
      
      # get all cells that do not conform to the translocation pattern
      if (dim(potential.tr)[1] != 0) {
        outlier.list <- lapply(1:nrow(potential.tr), get.outliers, phased = phased, hapmatrix = maj.cast, translocations = potential.tr)
      }
      
      ###########
      ## Output
      ###########
      
      cat("Writing output files...\n")
      for (x in 1:length(outlier.list)) {
        partnerA <- potential.tr[["segA"]][x]
        partnerB <- potential.tr[["segB"]][x]
        write.table(outlier.list[[x]], file.path(outlier.output, paste0(partnerA, "-", partnerB,"-outliers.txt")), quote=F, row.names = F)
      }
      
      write.table(p.values, file.path(tr.output, "pvalue-table.txt"), quote = F, row.names = F)
      write.table(potential.tr, file.path(tr.output, "translocations.txt"), quote = F, row.names = F)
      write.table(maj.cast, file.path(output.folder.cursamp, "haplotypes-per-arm.txt"), quote = F, row.names = F)
      write.table(recurrent, file.path(tr.output, "recurrent-segments.txt"), quote = F, row.names = F)
    }
    
    ##############################################
    ## Analysis: pq + suspected translocation file
    ##############################################
    
    if (options == "pq") {
      
      cat("Starting analysis\n")
      phased.pq <- pq.cast(phased)
      
      if (translocations == T) {
        
        if (exists(paste0(cursample,".tr"))){
          tryCatch({
            phased.pq <- merge(phased.pq, get(paste0(cursample,".tr"))[, -"sample"], by = "cell", all.x=T)
          }, warning = function(w){
            stop("Check that the 'sample' column of your translocations file matches the 'sample' column used by MosaiCatcher.\n")
          })
        }
      }
      
      ###############################
      ## Analysis: pq + regions file
      ###############################
      
      # if a 'regions' file is provided, use those to make extra potential translocations
      if (use.regions == T) {
        
        regions.cursample <- regions.f[sample == cursample]
        
        # only use the file if regions for the current sample are provided
        if (dim(regions.cursample)[1] != 0){
          
          translos <- region.to.tr(regions.cursample, phased)
          
          # convert translocation file to factors
          names <- colnames(translos[, -"cell"])
          translos[, paste0(names, ".H1") := lapply(.SD, function(x) {as.numeric(factor(substring(x, 0, 1)))}), .SDcols = names]
          translos[, paste0(names, ".H2") := lapply(.SD, function(x) {as.numeric(factor(substring(x, 2)))}), .SDcols = names]
          
          phased.pq <- merge(phased.pq, translos[, .SD, .SDcols = names(translos) %like% ".H[1-2]" | names(translos) %like% "cell"], by = "cell")
        }
      }
      
      recurrent <- rec.seg(phased)
      p.values <- get.pvalue.dt(phased.pq)
      potential.tr <- p.values[str_extract(segA, "chr[0-9X-Y]+") != str_extract(segB, "chr[0-9X-Y]+")][pBH < cutoff][order(pBH)]
      
      # get all cells that do not conform to the translocation pattern
      if (dim(potential.tr)[1] != 0) {
        outlier.list <- lapply(1:nrow(potential.tr), get.outliers, phased = phased, hapmatrix = phased.pq, translocations = potential.tr)
      }
      
      ###########
      ## Output
      ###########
      
      cat("Writing output files...\n")
      for (x in 1:length(outlier.list)) {
        partnerA <- potential.tr[["segA"]][x]
        partnerB <- potential.tr[["segB"]][x]
        write.table(outlier.list[[x]], file.path(outlier.output, paste0(partnerA, "-", partnerB,"-outliers.txt")), quote=F, row.names = F)
      }
      
      write.table(recurrent, file.path(tr.output, "recurrent-segments.txt"), quote = F, row.names = F)
      write.table(p.values, file.path(tr.output, "pvalue-table.txt"), quote = F, row.names = F)
      write.table(potential.tr, file.path(tr.output, "translocations.txt"), quote = F, row.names = F)
      write.table(phased.pq, file.path(output.folder.cursamp, "haplotypes-per-arm.txt"), quote = F, row.names = F)
    
    } # option = pq if statement
    
    if (options == "segments") {
      
      recurrent.segs <- rec.seg(phased)
      allsegs <- cast.haplotypes(phased)
      
      #################################
      ## Analysis: p-value calculations
      #################################
      
      p.values <- get.pvalue.dt(allsegs)
      potential.tr <- p.values[str_extract(segA, "chr[0-9X-Y]+") != str_extract(segB, "chr[0-9X-Y]+")][pBH < cutoff][order(pBH)]
      
      ###########
      ## Output
      ###########
      
      cat("Writing output files...\n")
      for (x in 1:length(outlier.list)) {
        partnerA <- potential.tr[["segA"]][x]
        partnerB <- potential.tr[["segB"]][x]
        write.table(outlier.list[[x]], file.path(outlier.output, paste0(partnerA, "-", partnerB,"-outliers.txt")), quote=F, row.names = F)
      }
      
      write.table(p.values, file.path(tr.output, "pvalue-table.txt"), quote = F, row.names = F)
      write.table(potential.tr, file.path(tr.output, "translocations.txt"), quote = F, row.names = F)
      write.table(phased.pq, file.path(output.folder.cursamp, "haplotypes-per-arm.txt"), quote = F, row.names = F)
      
    } # options = segments
    
    if (options == "make.segments") {
      
      ###############
      ## File checking and reading
      ###############
      
      # locate data directory
      translocator.data <- system.file(package = "TranslocatoR", "data", mustWork = T)
      
      # load blacklist
      if (blacklist == F) {
        cat("Proceeding without blacklist. This will affect performance negatively.\n")
      }
      else {
        blacklist <- fread(file.path(translocator.data, "blacklist.txt"))
      }
      
      # locate and read counts and segments files
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
      cat("\tRead counts file\n")
      
      tryCatch({
        segments <- fread(file.path(segments.folder, paste0(binsize, "_fixed_norm.txt")))
      }, warning = function(w) {
        stop(paste("Can't read the segments file", paste0(binsize, "_fixed_norm.txt"), "in", segments.folder, "\n"))
      })
      cat("\tRead segments file\n")
      
      #################################
      ## Analysis: cleaning and phasing
      #################################
      
      total.cells <- counts[, length(unique(cell))]
      cat(paste("Analyzing", total.cells, "cells in sample", cursample,"\n"))
      
      # remove all None bins and centromeric bins
      counts <- clean.bins(counts, blacklist, total.cells)
      
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
      
      final.segs <- unique(bp.segs)
      
      final.segs.H1H2 <- cast.haplotypes(final.segs)
      
      #################################
      ## Analysis: p-value calculations
      #################################
      
      p.values <- get.pvalue.dt(final.segs.H1H2)
      potential.tr <- p.values[str_extract(segA, "chr[0-9X-Y]+") != str_extract(segB, "chr[0-9X-Y]+")][pBH < cutoff][order(pBH)]
      
      #############
      ## Output
      #############
      
      write.table(potential.tr, file.path(tr.output, "translocations.txt"), quote = F, row.names = F)
      write.table(p.values, file.path(pval.output, "pvalue-table.txt"), quote = F, row.names = F)
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
