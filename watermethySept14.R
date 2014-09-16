## Sept 2014 version of wateRmelon test script
## watermethy.R uses wateRmelon package (version 1.0) to import data into a MethyLumiSet object
## Also creates a MethyLumiM object using the lumi package (v. 1.1.0)
## and a exprmethy450 object using IMA (v. 2.0).

library(lumi)
library(methylumi)
library('wateRmelon')

## Original Files
MethyFileName<-"~/methydata/MethyFileName.txt"
PhenoFileName<-"~/methydata/PhenoFileName.txt"
## Edited Files
MethyFile<-"~/methydata/MethyFile.txt"
PhenoFile<-"~/methydata/PhenoFile.txt"

# Read in the files
methyframe<-read.delim(MethyFileName, header=TRUE, stringsAsFactors=FALSE)
phenoframe<-read.delim(PhenoFileName, header=TRUE, stringsAsFactors=FALSE)

## Map SampleID's to some other field
# get rid of duplicates
phenoframe$SampleLabel[76]<-"Ctrl2"
# Create the ID map
id_map<-data.frame(phenoframe$SampleID, phenoframe$SampleLabel)

# Change the IDs in methyframe
for (i in 1:96) {
  names(methyframe)<-gsub(id_map[,1][i], id_map[,2][i], names(methyframe))
}

#Get rid of unusable pheno samples
phenoframe<-phenoframe[which(grepl("[0123456789]{6}[A|C]", phenoframe$SampleLabel)),]
#Get rid of sample that doesn't have Stim_med data
phenoframe<-phenoframe[which(!grepl("101892C", phenoframe$SampleLabel)),]

#Get rid of unusable samples from methy file
methy<-names(methyframe)[which(grepl("X[0123456789]{6}[A|C].*", names(methyframe)))]
notmethy<-names(methyframe)[-which(grepl("X[0123456789]{6}[A|C].*", names(methyframe)))]
notmethy<-notmethy[which(grepl("X[^_]{8}", notmethy))]
submethy<-setdiff(notmethy,methy)
methyframe<-methyframe[setdiff(names(methyframe),submethy)]
#Get rid of sample that doesn't have Stim_med data
methyframe<-methyframe[which(!grepl("101892C", names(methyframe)))]

# (add X's and rearrange phenoframe columns)
for (i in 1:91) {
  try(phenoframe$SampleLabel[i]<-paste("X",phenoframe$SampleLabel[i],sep=""), silent=TRUE)
}
#There is an error that says: 'Error in `$<-.data.frame`(`*tmp*`, "SampleLabel", value = c("X107612C",  : 
#replacement has 91 rows, data has 90' - this doesn't change the output though so I silenced it

#change column labels
phenoframe<-data.frame(phenoframe$SampleLabel, phenoframe)
names(phenoframe)[2]<-"OldSampleID"
names(phenoframe)[1]<-"SampleID"

## Write to files
write.table(methyframe, file=MethyFile, sep="\t", quote=FALSE, row.names=FALSE)
write.table(phenoframe, file=PhenoFile, sep="\t", quote=FALSE, row.names=FALSE)

## create the exprmethy450 object using IMA
library(IMA)
imadata<-IMA.methy450R(fileName = MethyFile,
                       columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),
                       groupfile = PhenoFile)

## change column headers "Detection.Pval" to "Detection Pval"
names(methyframe)<-gsub("Detection.Pval", "Detection Pval", names(methyframe))
## Write to file
write.table(methyframe, file=MethyFile, sep="\t", quote=FALSE, row.names=FALSE)

## create the MethyLumiM object without pheno data
lmdata<-lumiMethyR(MethyFile)

## create the MethyLumiSet object without pheno data
mldata<-methylumiR(filename=MethyFile)

# change M/F to 0/1
phenoframe$Sex <- gsub("M", 0, phenoframe$Sex)
phenoframe$Sex <- gsub("F", 1, phenoframe$Sex)

# now add pheno data
pData(mldata)<-phenoframe
pData(lmdata)<-phenoframe

# filter and normalize 2 ways
mldata.pf<-pfilter(mldata)
mldata.dasen<-dasen(mldata.pf)
mldata.nasen<-nasen(mldata.pf)

#Check standard errors
dmrse_row(mldata.pf)
dmrse_row(mldata.dasen)
dmrse_row(mldata.nasen)

#this will give 3 values that will need to be averaged (wateRmelon people just take the mean)
genki(mldata.pf)
genki(mldata.dasen)
genki(mldata.nasen)

# Boxplots of methylated/unmethylated/betas reveal distribution of normalized values

# Can't calculate X-chromosome metrics on QC'd betas with seabi method with all male samples

## Beth says:
## After choosing the normalized betas, create a new IMA object;  this can be done if you need to change the phenofile for any reason
#newobj<-new("exprmethy450",bmatrix=betas(normalizeddata),detectP=pvals(normalizeddata),annot=annot.matrix,groupinfo=fromnormaliec data)
