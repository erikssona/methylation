source("http://bioconductor.org/biocLite.R")
biocLite("FDb.InfiniumMethylation.hg19")

#fix sampleName problem
sampleNames(lmdata)<-t(phenoframe[1])
lmdata.analysis<-MethyLumiM2GenoSet(lmdata)
mldata.mlumi.analysis<-MethyLumiM2GenoSet(mldata.mlumi)
# At this point lmdata.analysis and mldata.mlumi.analysis can be directly compared,
# with lmdata originating from lumi and mldata from watermelon.

## get sample type information
sampleType <- pData(mldata.mlumi.analysis)$SampleType
## try doing smoothing first
mldata.mlumi.analysis.sm <- smoothMethyData(mldata.mlumi.analysis, winSize = 250)

## Do differential test (this gives an error currently)
allResult <- detectDMR.slideWin(mldata.mlumi.analysis.sm, sampleType=sampleType, testmethod="Wilcox")

## get sample type information
sampleType <- pData(lmdata.analysis)$SampleType
## Do differential test
allResult <- detectDMR.slideWin(lmdata.analysis, sampleType=sampleType)
