process_methyl_microarray <- function(rgSet){
  
  # Detect p values
  p_values <- detectionP(rgSet, type = "m+u")
  
  mSet <- preprocessNoob(rgSet)
  
  # ensure probes are in the same order in the mSetSq and detP objects
  p_values <- p_values[match(featureNames(mSet),rownames(p_values)),]
  
  # remove any probes that have failed in one or more samples; this next line
  # checks for each row of p_values whether the number of values < 0.01 is equal
  # to the number of samples (TRUE) or not (FALSE)
  keep <- rowSums(p_values < 0.01) == ncol(mSet)
  table(keep)
  
  # Subset the GenomicRatioSet
  mSet <- mSet[keep,]
  
  # Remove probes with known SNPs
  gmSet <- gmSet <- dropLociWithSnps(mapToGenome(mSet))
  
  # Remove cross reactive probes
  xreactive_probes <- xreactive_probes(array_type="EPIC")
  keep <- !(featureNames(gmSet) %in% xreactive_probes)
  gmSet <- gmSet[keep,]
  
  # Remove probes mapping to sex chromosomes
  ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  keep <- !(featureNames(gmSet) %in% ann$Name[ann$chr %in% c("chrX", "chrY")])
  gmSet <- gmSet[keep, ]
  
  # make gmSet available in the global environment so M values and Beta values can be saved
  return(gmSet)
}

# Single nucleotide polymorphism (SNP) are biological markers of DNA variation. Removal reduces ambiguity, reduces false positives and increases power. However, it may lead to exclusion of important regulatory regions.

# Beta values, which range between [0,1], represent the proportion of methylated DNA copies at a specific CpG site compared to all of the DNA copies at that site. Beta values for the ith CpG site are mathematically defined as,

# Beta(i) = max(yi(methyl),0) / [max(yi(unmethyl), 0) + max(yi(methyl), 0) + alpha],

# where yi methyl and unmethyl are the intensities measured by the ith methylated and unmethylated probes, respectively. 0 ensures that negative signals that take place after background adjustment are not considered. Alpha, a constant (usually 100), serves as an offset to regularize beta when both methylated and unmethylated instensities are low. Beta values are intuitive, and hence useful for visual presentations of the reads. (Du et al., 2010)
# M values are log transformed beta values. Negative M values indicate higher levels of unmethylated DNA copies than methylated DNA copies. Postive M values indicate higher levels of methylated DNA copies than unmethylated. M values close to zero indicate similar intensities between methylated and unmethylated probes. M values are mathematically expressed as,

# M(i) = log2([max(yi(methyl))+alpha]/[max(yi(unmethyl))+alpha])

# where alpha is a constant to offset dramatic changes resulting from small estimation errors. M values can be expressed in terms of beta values,

# M(i) = log2(Beta(i)/[1 - Beta(i)])

# and therefore, beta values can be expressed in terms of M values,

# Beta(i) = 2\^M(i) / [2\^M(i) + 1]

# M values are more statistically sound than beta values.