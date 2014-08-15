## Data cleaning/preprocessing before conducting PCA on the methylation dataset

CovariateFileName <- "~/methydata/covariates.txt"
covariates<-read.delim(CovariateFileName, header=TRUE, stringsAsFactors=TRUE)

overlap <- na.omit(covariates[which(intersect(covariates$Comments, phenoframe$SampleID)>0),])

final<-cbind(overlap[1], overlap[4], overlap[9:13], overlap[15:17])

# Example of PCA on 2 columns
prcomp(final[2:10])



source("http://bioconductor.org/biocLite.R")
biocLite("pcaMethods")
library(pcaMethods)
