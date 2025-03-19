# Illumina-EPICv1-DNA-Methylation-Analysis
Pre-processing, differential methylation analysis, differential variation analysis, and pathway analysis of Illumina EPICv1 microarray data using R.
The two R files (process_microarry and region_analysis) are user defined functions executed in the R markdown files.

Current state of the files are messy in explanation and execution. This is my first time analyzing DNA methylation data, so this is a work in progress. Updates and improvements will be made regularly.

The following sources have been essential in allowing me to complete these projects;

For cleaning and processing the raw data:
https://github.com/raulsanzr/DNA-Methylation/blob/main/R/human_epic_analysis.R

For differential analyses and pathway analyses:
https://www.bioconductor.org/packages/devel/bioc/vignettes/missMethyl/inst/doc/missMethyl.html#annotation-for-illumina-methylationepic-v2.0-beadchip-datasets

For differential analysis while adjusting for unwanted variation, when RUV was unsuccessful:

Leek JT and Storey JD. (2008) A general framework for multiple testing
dependence. Proceedings of the National Academy of Sciences , 105:
18718-18723.

LeekJTandStoreyJD.(2007)Capturingheterogeneityingeneexpression
studies by ‘Surrogate Variable Analysis’. PLoS Genetics, 3: e161.
