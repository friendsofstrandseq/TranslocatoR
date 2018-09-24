#' Find translocations
#' 
#' Find translocated chromosomal segments by correlating haplotype strand states in Mosaicatcher output. Fully vectorized approach heavily reliant on data.table.
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
#' TODO: incorporate discreteMTP calc
#' TODO: option to look at one or more segments (instead of seeing which segments are characterized by strand changes, get a segment and look which state it has in each chrom)

### Libraries --> move to namespace
library(data.table)
library(gtools)
library(GenomicRanges) # necessary?
library(rprojroot)
library(dplyr)
library(ggplot2)
library(stringr)

translocatoR <- function(data.folder, samples, output.folder, binsize = 100000L, trfile = NULL, blacklist = T) {
  
  ###############################
  ### File checking, initializing
  ###############################
  
  # locate root directory
  root.dir <- rprojroot::find_root("translocator.R")
  # load helper functions
  source(file.path(root.dir, "utils/helpers.R"))
  
  if (blacklist == F) {
    cat("Proceeding without blacklist. This will affect performance negatively.")
  }
  else {
    blacklist <- fread(file.path(root.dir, "data/blacklist.txt"))
  }
  
  # default binsize 100kb, check if integer
  if(class(binsize) != "integer") {
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
    
    # identify stretches of consecutive states (taken from Venla's SCE detector) -- remove here?
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
    rm(switch.counts)
    
    ### Using MosaiCatcher segments to segment the data and get the haplotype-aware state for each segment
    
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
    newseg <- newseg[, .SD[,cbind(.(start = c(i.start, start, end)), .(end = c(start, end, i.end)))], by= .(chrom, start, end)]
    setnames(newseg, c("start", "end", "V1", "V2"), c("old.start", "old.end", "start", "end"))
    # now overlap: remove old segment(s) and replace new at the same time
    setkey(newseg, chrom, start, end)
    setkey(segments.many, chrom, start, end)
    segments.many <- foverlaps(segments.many, newseg)
    segments.many[start >= i.start & end <= i.end, c("i.start", "i.end"):=.(start, end)]
    # bug: if a new segments had the same start/end as an old one you get a segment of size zero i.e. start == end
    segments.many <- segments.many[,.(chrom, i.start, i.end)]
    setnames(segments.many, c("i.start", "i.end"), c("start", "end"))
    
    setkey(segments.many, chrom, start, end)
    setkey(counts, chrom, start, end)
    segment.counts <-foverlaps(counts, segments.many)
    segment.counts <- segment.counts[order(cell,chrom, start, end)]
    
    ## print how many cells had 2 or more choices?
    segment.class <- segment.counts[, .(nbins = max(.N)), by = .(sample, cell, chrom, start, end, class)][, .(class_ = class[which.max(nbins)]), by=.(sample, cell, chrom, start, end)]
    
    # now making own segments again... change!
    ## extract often-changing segments?
    segment.class <- segment.class[order(sample, cell, chrom, start, end)]
    segment.class$cnsc <- cumsum(segment.class[, .(consecutive = c(1,abs(diff(as.numeric(factor(class_)))))), by = .(sample, cell, chrom)]$consecutive)
    bp.segs <- segment.class[, .SD[c(1,.N)], by=cnsc]
    bp.segs[, c("start","end"):=.(min(start), max(end)), by=cnsc]
    bp.segs <- bp.segs[,.SD[1], by=.(sample, cell, chrom, start, end)]
    
    ## section may not be important ###
    bps <- bp.segs[,.N, by=.(chrom, start, end)][mixedorder(chrom)]
    bps[, whole.chrom := start==min(start) & end == max(end), by=chrom]
    bps <- bps[whole.chrom==F]
    bps[, size := end - start]
    bps.pot <- bps[N >= total.cells*0.1 & size > 1000000] # most important recurrent bps: are larger than 1 Mb, occur at least in 10% of cells
    ## remove above? ###
    
    unique.bps <- bp.segs[,unique(c(start, end)) ,by=.(cell, chrom)][, .(unique.bp = sapply(.SD, toString)), by=.(cell, chrom)]
    bp.segs <- merge(bp.segs, unique.bps)
    
    bp.segs[, cells.unique.bp := length(unique(cell)), by=.(chrom, unique.bp)]
    # select cells where bp is <= 2 (bp occurs max twice in sample)
    #bp.segs[cells.unique.bp<=2]
    bp.segs <- merge(bp.segs, bps, by=c("chrom", "start", "end"))
    bp.segs <- bp.segs[order(cell, chrom, start, end)]
    bp.segs[, len := end - start]
    
    bp.segs[N <= 2, c("start","end", "majority"):=.(min(start), max(end), .SD[which.max(len), class_]), by=.(cell, chrom)]
    final.segs <- unique(bp.segs[is.na(majority) | class_==majority, .(sample, cell, chrom, start, end, class_)])
    
    final.segs[, label:=paste0(chrom, ":", start, "-", end)]
    final.segs[,c("H1","H2"):=.(substring(class_, 0,1), substring(class_,2))]
    final.segs[,c("H1_factor","H2_factor"):=.(as.numeric(factor(H1)),as.numeric(factor(H2)))]
    final.segs.castH1 <- dcast(final.segs, cell~label,value.var = 'H1_factor', drop=F)
    setcolorder(final.segs.castH1, mixedorder(names(final.segs.castH1)))
    final.segs.castH2 <- dcast(final.segs, cell~label,value.var = 'H2_factor', drop=F)
    setcolorder(final.segs.castH2, mixedorder(names(final.segs.castH2)))
    
    setnames(final.segs.castH1, names(final.segs.castH1)[-1], paste0(names(final.segs.castH1)[-1], ".H1"))
    setnames(final.segs.castH2, names(final.segs.castH2)[-1], paste0(names(final.segs.castH2)[-1], ".H2"))
    final.segs.H1H2 <- cbind(final.segs.castH1, final.segs.castH2[,-1])
    
    spearman.dt <- getcor(final.segs.H1H2) # silence warnings
    test <- sapply(1:nrow(spearman.dt), function(x){if (x%%1000 == 0){cat(paste(x, "combinations out of", nrow(spearman.dt),"calculated\n"))}
      p.helper(final.segs.H1H2[,.(get(spearman.dt[["Var1"]][x]), get(spearman.dt[["Var2"]][x]))])})
    test2 <- as.data.table(matrix(test, nrow = nrow(spearman.dt), ncol = 2, byrow = T))
    setnames(test2, c("V1", "V2"), c("or", "p"))
    
    pvals <- cbind(spearman.dt, test2)
    pvals[,log10.p := -log10(p)]
    pvals[,BH.p := p.adjust(p, "BH")]
    pvals[,B.p := p.adjust(p, "bonferroni")]
    
    #############
    ### Plotting
    #############
    
    ## plot by haplotype (H1-H1, H1-H2, H2-H1, H2-H2)
    
    ## Interchromosomal translocation candidates
    tr.can <- pvals[BH.p < 0.01][str_extract(Var1, 'chr[0-9X-Y]+') != str_extract(Var2, 'chr[0-9X-Y]+')]
    
    ## Segments in consistent order, with lowest-number chrom first
    ## watch for chrX/Y
    tr.can[, c("chrom1", "chrom2") := .(as.numeric(str_extract(Var1, '\\d+')), as.numeric(str_extract(Var2, '\\d+')))]
    tr.can[chrom1 > chrom2, c("Var2", "Var1") := .(Var1, Var2)]
    tr.can <- tr.can[mixedorder(Var1)]
    tr.can[, c("chrom1", "chrom2") := .(as.numeric(str_extract(Var1, '\\d+')), as.numeric(str_extract(Var2, '\\d+')))]
    
    ## get segments back as numeric
    tr.can[, c("start.chrom1", "end.chrom1", "start.chrom2", "end.chrom2") := 
             .(as.numeric(str_extract(Var1, '(?<=:)\\d+')), as.numeric(str_extract(Var1, '(?<=-)\\d+')), 
               as.numeric(str_extract(Var2, '(?<=:)\\d+')), as.numeric(str_extract(Var2, '(?<=-)\\d+')))]
    tr.can[, plottogether := .GRP, by=.(chrom1, chrom2)]
    tr.can[, haplo :=paste0(str_extract(Var1, "H[1-2]"), "-",str_extract(Var2, "H[1-2]"))]
    tr.can[cor > 0, posneg := "positive"]
    tr.can[cor < 0, posneg := "negative"]
    ## ggplot takes only data frames
    plotdf <- setDF(tr.can)
    
    ## updating defaults is the only way to get the text to change colour
    update_geom_defaults("text", list(colour = "white"))
    for (i in unique(plotdf$plottogether)) {
      curpartners <- paste0(unique(plotdf[plotdf$plottogether == i,"chrom1"]),"-", unique(plotdf[plotdf$plottogether == i,"chrom2"]))
      ggplot(plotdf[plotdf$plottogether == i,]) + geom_tile(aes(Var1, Var2, fill = posneg))+ 
        geom_text(aes(Var1, Var2, label=formatC(as.numeric(BH.p), format = "e", digits = 2)))+
        scale_fill_manual(values = c("positive" = "darkgreen", "negative" = "darkred"))+
        facet_wrap(~haplo)+
        theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 8, hjust = 1))
      ggsave(filename = paste0(sample, "_co-segregation_", curpartners,".pdf"), plot = last_plot(), path = output.folder)
    }
    ## sample loop 
  }
  
## end of function  
}  
