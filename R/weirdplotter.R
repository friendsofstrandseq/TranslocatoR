## plot by haplotype (H1-H1, H1-H2, H2-H1, H2-H2)

## Interchromosomal translocation candidates
stupidplot <- function(x){tr.can <- pvals[pBH < 0.01][str_extract(segA, 'chr[0-9X-Y]+') != str_extract(segB, 'chr[0-9X-Y]+')]

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
}