#' Find translocations
#' 
#' Find translocated chromosomal segments by correlating haplotype strand states.
#' 
#' @param data.folder MosaiCatcher output folder
#' @param samples on which samples in the MosaiCatcher output folder would you like to run TranslocatoR
#' @param trfile list of textfiles with manually determined translocation states if desired, must contain cell ID and approximate start/end
#' @param output.folder absolute path to output folder for TranslocatoR data
#' @return matrix
#' 
#' @author Alex van Vliet
#' 

translocatoR <- function(data.folder, samples, trfile, output.folder, binsize) {
  
  ### Libraries
  library(data.table)
  library(gtools)
  library(GenomicRanges)
  
  ### File checking, initializing
  
  # default binsize 100kb
  if(missing(binsize)) {
    binsize <- 100000
  }
  
  data.folder <- tools::file_path_as_absolute(data.folder)
  if(!dir.exists(data.folder)) {
    stop("Please check that you provided the correct path for the data folder. Stopping execution.")
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
      stop(paste("The folder with read counts for sample", sample, "is not in the data folder. Stopping execution."))
    }
    
    phased.folder <- file.path(data.folder, "strand_states", sample)
    if(!dir.exists(phased.folder)) {
      stop(paste("The folder with phased haplotypes for sample", sample, "is not in the data folder. Stopping execution."))
    }
    
    segments.folder <- file.path(data.folder, "segmentation", sample)
    if(!dir.exists(segments.folder)) {
      stop(paste("The folder with segments for sample", sample, "is not in the data folder. Stopping execution."))
    }
    
    # use normalized counts always
    counts <- fread(paste("zcat", file.path(counts.folder, paste0(binsize, "_fixed_norm.txt.gz")))
    
  ## sample loop 
  }
  
## end of function  
}  