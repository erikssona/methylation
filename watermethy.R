## watermethy.R uses wateRmelon package (version 1.0) to import data into a MethyLumiSet object
## Also creates a MethyLumiM object using the lumi package (v. 1.1.0)
## and a exprmethy450 object using IMA (v. 2.0).

## Original Files
MethyFileName<-"~/methydata/MethyFileName.txt"
PhenoFileName<-"~/methydata/PhenoFileName.txt"
## Edited Files
MethyFile<-"~/methydata/MethyFile.txt"
PhenoFile<-"~/methydata/PhenoFile.txt"

# Read in the files
methyframe<-read.delim(MethyFileName, header=TRUE, stringsAsFactors=FALSE)
phenoframe<-read.delim(PhenoFileName, header=TRUE, stringsAsFactors=FALSE)

# Map SampleID's to some other field
# get rid of duplicates
#levels(phenoframe$SampleLabel)<-c(levels(phenoframe$SampleLabel), "Ctrl2")
phenoframe$SampleLabel[76]<-"Ctrl2"
# Create the ID map
id_map<-data.frame(phenoframe$SampleID, phenoframe$SampleLabel)

# Change the IDs in methyframe
for (i in 1:96) {
  names(methyframe)<-gsub(id_map[,1][i], id_map[,2][i], names(methyframe))
}

#Get rid of unusable pheno samples
phenoframe<-phenoframe[which(grepl("[0123456789]{6}[A|C]", phenoframe$SampleLabel)),]

#Get rid of unusable samples from methy file
methy<-names(methyframe)[which(grepl("X[0123456789]{6}[A|C].*", names(methyframe)))]
notmethy<-names(methyframe)[-which(grepl("X[0123456789]{6}[A|C].*", names(methyframe)))]
notmethy<-notmethy[which(grepl("X[^_]{8}", notmethy))]
submethy<-setdiff(notmethy,methy)
methyframe<-methyframe[setdiff(names(methyframe),submethy)]


# (add X's and rearrange phenoframe columns)
#for (i in 1:96) {
#  levels(phenoframe$SampleLabel)[i]<-paste("X",levels(phenoframe$SampleLabel)[i],sep="")
#}
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
library(lumi)
lmdata<-lumiMethyR(MethyFile)

## create the MethyLumiSet object without pheno data
library(methylumi)
mldata<-methylumiR(filename=MethyFile)

### workspace saved 0728
# change M/F to 0/1
phenoframe$Sex <- gsub("M", 0, phenoframe$Sex)
phenoframe$Sex <- gsub("F", 1, phenoframe$Sex)

# now add pheno data
pData(mldata)<-phenoframe
pData(lmdata)<-phenoframe



library('wateRmelon')
mldata.pf<-pfilter(mldata)
mldata.dasen.pf<-dasen(mldata.pf)
mldata.nasen.pf<-nasen(mldata.pf)

#Check standard errors
dmrse_row(mldata.pf)
dmrse_row(mldata.dasen.pf)
dmrse_row(mldata.nasen.pf)

#this will give 3 values that will need to be averaged (wateRmelon people just take the mean)
genki(mldata.pf)
genki(mldata.dasen.pf)
genki(mldata.nasen.pf)

# Can't calculate X-chromosome metrics on QC'd betas with seabi method with all male samples


boxplot(log(methylated(mldata)), las=2, cex.axis=0.8, main="methylated, unfiltered" )
boxplot(log(methylated(mldata.pf)), las=2, cex.axis=0.8, main="methylated, filtered" )
boxplot(log(methylated(mldata.dasen.pf)), las=2, cex.axis=0.8, main="methylated, dasen" )
boxplot(log(methylated(mldata.nasen.pf)), las=2, cex.axis=0.8, main="methylated, nasen" )
boxplot(log(unmethylated(mldata)), las=2, cex.axis=0.8, main="unmethylated, unfiltered" )
boxplot(log(unmethylated(mldata.pf)), las=2, cex.axis=0.8, main="unmethylated, filtered" )
boxplot(log(unmethylated(mldata.dasen.pf)), las=2, cex.axis=0.8, main="unmethylated, dasen" )
boxplot(log(unmethylated(mldata.nasen.pf)), las=2, cex.axis=0.8, main="unmethylated, nasen" )

# fix name problem
sampleNames(phenoData(mldata))<-sampleNames(assayData(mldata))
# Make a MethyLumiM object based on the mldata MethyLumiSet object (watermelon)
mldata.mlumi <- as(mldata, 'MethyLumiM')

###


## Beth says:
## After choosing the normalized betas, create a new IMA object;  this can be done if you need to change the phenofile for any reason
#newobj<-new("exprmethy450",bmatrix=betas(normalizeddata),detectP=pvals(normalizeddata),annot=annot.matrix,groupinfo=fromnormaliec data)
