region_analysis <- function(contrast_coeff, gene_database, betacutoff, fdr){
  myAnnotation <- cpg.annotate(object = Mval, datatype = "array", what = "M", 
                               arraytype = c("EPICv1"), 
                               analysis.type = "differential", design = design, contrasts = T, cont.matrix = contrast,
                               fdr = fdr, coef = contrast_coeff)
  DMRs <- suppressMessages(dmrcate(myAnnotation, C=2, betacutoff = betacutoff))
  results.ranges <- extractRanges(DMRs, genome = "hg19")
  
  # plot DMRs
  if disease == disease
  cols <- c(1,2)[disease]
  names(cols) <-disease
  par(mfrow=c(1,1))
  DMR.plot(ranges=results.ranges, dmr=2, CpGs=beta_values, phen.col=cols, 
           what="Beta", arraytype="EPICv1", genome="hg19")
  gst.region <- goregion(results.ranges, all.cpg=rownames(Mval), 
                         collection= gene_database, array.type="EPIC", plot.bias = T)
  topGSA(gst.region)
}