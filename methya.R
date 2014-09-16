source("http://bioconductor.org/biocLite.R")
biocLite("FDb.InfiniumMethylation.hg19")
library(methyAnalysis)

#fix sampleName problem in lumi objects
sampleNames(lmdata)<-t(phenoframe[1])
## normalize lumi objects
library(lumi)
lmdata.quantile <- normalizeMethylation.quantile(lmdata)
lmdata.ssn <- normalizeMethylation.ssn(lmdata)


# At this point lmdata.analysis and mldata.mlumi.analysis can be directly compared,
# with lmdata originating from lumi and mldata from watermelon.

## Now we're going to do the DMR test on normalized data
## first with dasen (watermelon)
#fix sampleName problem
sampleNames(phenoData(mldata.dasen.pf))<-sampleNames(assayData(mldata.dasen.pf))
# Make a MethyLumiM object based on the mldata MethyLumiSet object (watermelon)
mldata.dasen.mlumi <- as(mldata.dasen.pf, 'MethyLumiM')
mldata.dasen.mlumi.analysis<-MethyLumiM2GenoSet(mldata.dasen.mlumi)
## get sample type information
mldata.dasen.mlumi.sampleType <- pData(mldata.dasen.mlumi.analysis)$status
## try doing smoothing first
mldata.dasen.mlumi.analysis.sm <- smoothMethyData(mldata.dasen.mlumi.analysis, winSize = 250)
## Do differential test
mldata.dasen.mlumi.Result <- detectDMR.slideWin(mldata.dasen.mlumi.analysis.sm, sampleType=mldata.dasen.mlumi.sampleType, testmethod="Wilcox")
mldata.dasen.mlumi.DMRinfo = identifySigDMR(mldata.dasen.mlumi.Result, topNum=100)

# Now with nasen (watermelon)
#fix sampleName problem
sampleNames(phenoData(mldata.nasen.pf))<-sampleNames(assayData(mldata.nasen.pf))
# Make a MethyLumiM object based on the mldata MethyLumiSet object (watermelon)
mldata.nasen.mlumi <- as(mldata.nasen.pf, 'MethyLumiM')
# Make a MethyGenoSet object from the MethyLumiM object (lumi)
mldata.nasen.mlumi.analysis<-MethyLumiM2GenoSet(mldata.nasen.mlumi)
## get sample type information
mldata.nasen.mlumi.sampleType <- pData(mldata.nasen.mlumi.analysis)$status
## try doing smoothing first
mldata.nasen.mlumi.analysis.sm <- smoothMethyData(mldata.nasen.mlumi.analysis, winSize = 250)
## Do differential test
mldata.nasen.mlumi.Result <- detectDMR.slideWin(mldata.nasen.mlumi.analysis.sm, sampleType=mldata.nasen.mlumi.sampleType, testmethod="Wilcox")
mldata.nasen.mlumi.DMRinfo = identifySigDMR(mldata.nasen.mlumi.Result, topNum=100)

# Now with lmdata (lumi) quantile normalized
# Make a MethyGenoSet object from the MethyLumiM object (lumi)
lmdata.quantile.analysis<-MethyLumiM2GenoSet(lmdata.quantile)
## get sample type information
lmdata.quantile.sampleType <- pData(lmdata.quantile.analysis)$status
## try doing smoothing first
lmdata.quantile.analysis.sm <- smoothMethyData(lmdata.quantile.analysis, winSize = 250)
## Do differential test
lmdata.quantile.Result <- detectDMR.slideWin(lmdata.quantile.analysis.sm, sampleType=lmdata.quantile.sampleType, testmethod="Wilcox")
lmdata.quantile.DMRinfo = identifySigDMR(lmdata.quantile.Result, topNum=100)

# Now with lmdata (lumi) ssn normalized
# Make a MethyGenoSet object from the MethyLumiM object (lumi)
lmdata.ssn.analysis<-MethyLumiM2GenoSet(lmdata.ssn)
## get sample type information
lmdata.ssn.sampleType <- pData(lmdata.ssn.analysis)$status
## try doing smoothing first
lmdata.ssn.analysis.sm <- smoothMethyData(lmdata.ssn.analysis, winSize = 250)
## Do differential test
lmdata.ssn.Result <- detectDMR.slideWin(lmdata.ssn.analysis.sm, sampleType=lmdata.ssn.sampleType, testmethod="Wilcox")
lmdata.ssn.DMRinfo = identifySigDMR(lmdata.ssn.Result, topNum=100)



### Next: Annotations
## Annotate significant DMR info
lmdata.quantile.DMRinfo <- annotateDMRInfo(lmdata.quantile.DMRinfo, 'TxDb.Hsapiens.UCSC.hg19.knownGene')

