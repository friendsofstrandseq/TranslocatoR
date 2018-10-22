#' Edit Mosaicatcher segments if necessary
#' 
#' @param mc.segments segments as determined by MosaiCatcher (with lowest sse)
#' @param rec.segments recurrent segemnts determined from count data
#' 
#' @return mc.segments with new segments added if necessary
#' @author Alex van Vliet
#' 
#' @export
#' 

# TODO: majority: >= 0.5 so must be in line with rest of code

segment.check <- function(mc.segments, rec.segments) {
  
  setkey(rec.segments, chrom, start, end)
  # by how many MosaiCatcher segments is each recurrent segment represented
  whichseg <- foverlaps(mc.segments, rec.segments)
  whichseg[, seg.id := paste0(i.start, "-", i.end)]
  howmanyseg <- na.omit(whichseg[, .SD[,uniqueN(seg.id)], by = .(chrom, start, end)])
  setnames(howmanyseg, "V1", "numseg")
  
  # if the recurrent segment falls in just one MosaiCatcher segment and is not the majority it can never be represented
  # so a final check is to see if that one segment covers the majority of the region
  # if not, make a new segment
  newseg <- howmanyseg[numseg == 1]
  newseg <- merge(newseg, whichseg, by=c("chrom", "start", "end"))
  newseg[, seg.len := i.end - i.start]
  newseg[, majority := len/seg.len >= 0.5]
  
  # if there are new segments that are not represented by MC segmentation: incorporate those
  if (newseg[, !all(majority)]) {
    newseg <- newseg[majority==F, .(start = min(start), end = max(end)),by=.(chrom, i.start, i.end)]
    newseg <- newseg[, .SD[,cbind(.(start = c(i.start, start, end)), .(end = c(start, end, i.end)))], by= .(chrom, start, end)]
    setnames(newseg, c("start", "end", "V1", "V2"), c("old.start", "old.end", "start", "end"))
    
    # now overlap: remove old segment(s) and replace new at the same time
    setkey(newseg, chrom, start, end)
    setkey(mc.segments, chrom, start, end)
    mc.segments <- foverlaps(mc.segments, newseg)
    mc.segments[start >= i.start & end <= i.end, c("i.start", "i.end"):=.(start, end)]
    
    # bug: if a new segments had the same start/end as an old one you get a segment of size zero i.e. start == end
    # will be resolved after overlaps in the next stage
    mc.segments <- unique(mc.segments[,.(chrom, i.start, i.end)])
    setnames(mc.segments, c("i.start", "i.end"), c("start", "end"))
    
    return(mc.segments)
  }
  
  else {
    return(mc.segments)
  }
  
}
