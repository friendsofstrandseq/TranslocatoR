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
  
  if (options == "segments" & (!is.null(trfile) | !is.null(regions))) {
    cat("The option \"segments\" does not incorporate trfile or region information. Proceeding without.")
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
          write.table(maj.cast, file.path(output.folder.cursamp, "haplotypes-per-chrom.txt"), quote = F, row.names = F)
        }
      }
      
      recurrent <- rec.seg(phased)
      p.values <- get.pvalue.dt(maj.cast)
      potential.tr <- p.values[str_extract(segA, "chr[0-9X-Y]+") != str_extract(segB, "chr[0-9X-Y]+")][pBH < cutoff][order(pBH)]
      
      # get all cells that do not conform to the translocation pattern
      if (dim(potential.tr)[1] != 0) {
        outlier.list <- lapply(1:nrow(potential.tr), get.outliers, phased = phased, hapmatrix = maj.cast, translocations = potential.tr)
      }
    } # options = pq if statement
    
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
          write.table(phased.pq, file.path(output.folder.cursamp, "haplotypes-per-arm.txt"), quote = F, row.names = F)
        }
      }
      
      # get recurrent segments, p-value matrix and potential translocations based on the corrected p-value (below cutoff)
      recurrent <- rec.seg(phased)
      p.values <- get.pvalue.dt(phased.pq)
      potential.tr <- p.values[str_extract(segA, "chr[0-9X-Y]+") != str_extract(segB, "chr[0-9X-Y]+")][pBH < cutoff][order(pBH)]
      
      # get all cells that do not conform to the translocation pattern
      outliers <- FALSE
      if (dim(potential.tr)[1] != 0) {
        outliers <- TRUE
        outlier.list <- lapply(1:nrow(potential.tr), get.outliers, phased = phased, hapmatrix = phased.pq, translocations = potential.tr)
      }
    
    } # option = pq if statement
    
    ###############################
    ## Analysis: segments
    ###############################
    
    if (options == "segments") {
      
      recurrent <- rec.seg(phased)
      allsegs <- cast.haplotypes(phased)
      write.table(allsegs, file.path(output.folder.cursamp, "haplotypes-per-segment.txt"), quote = F, row.names = F)
      
      p.values <- get.pvalue.dt(allsegs)
      potential.tr <- p.values[str_extract(segA, "chr[0-9X-Y]+") != str_extract(segB, "chr[0-9X-Y]+")][pBH < cutoff][order(pBH)]
      
      # get all cells that do not conform to the translocation pattern
      outliers <- FALSE
      if (dim(potential.tr)[1] != 0) {
        outliers <- TRUE
        outlier.list <- lapply(1:nrow(potential.tr), get.outliers, phased = phased, hapmatrix = allsegs, translocations = potential.tr)
      }
      
    } # options = segments
    
    ###########
    ## Output
    ###########
    
    cat("Writing output files...\n")
    if (outliers == TRUE) {
      for (x in 1:length(outlier.list)) {
        partnerA <- potential.tr[["segA"]][x]
        partnerB <- potential.tr[["segB"]][x]
        write.table(outlier.list[[x]], file.path(outlier.output, paste0(partnerA, "-", partnerB,"-outliers.txt")), quote=F, row.names = F)
      }
    }
    
    write.table(recurrent, file.path(tr.output, "recurrent-segments.txt"), quote = F, row.names = F)
    write.table(p.values, file.path(tr.output, "pvalue-table.txt"), quote = F, row.names = F)
    write.table(potential.tr, file.path(tr.output, "translocations.txt"), quote = F, row.names = F)
    
    #############
    ### Plotting
    #############
    
    ## sample loop 
  }
  
## end of function  
}  
