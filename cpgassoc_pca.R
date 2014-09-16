#######
## Is methylation (measured by betas) associated with
## phenotype status (ADHD/ctrl)?
## To answer this we use cpg.assoc to test the basic model
## methylation ~ phenostatus, then other models that include
## various covariates.
## One covariate is the principal component of PCA on the betas.
############################

source("http://bioconductor.org/biocLite.R")
biocLite("pcaMethods")
library(pcaMethods)## Data cleaning/preprocessing before conducting PCA on the methylation dataset
library("CpGassoc")

# Get the covariates
CovariateFileName <- "~/methydata/covariates.txt"
covariates<-read.delim(CovariateFileName, header=TRUE, stringsAsFactors=TRUE)
#add row names for merging
rownames(covariates)<-covariates$Comments

####### This was wrong, actually had to go back to original dataframes
# This was how I took care of the problem sample (must do this with all MethyLumiSet objects before/after normalizing)
#pData(mldata.dasen.pf)<-pData(mldata.pf)[pData(mldata.dasen.pf)[1]!=outcast,]
#betas(mldata.dasen.pf)<-betas(mldata.dasen.pf)[,colnames(betas(mldata.dasen.pf))!=outcast]

## Do PCA on betas
pca_samples<-prcomp(betas(mldata.dasen.pf),center=TRUE,scale=TRUE)
pca_loadings<-as.data.frame(pca_samples$rotation)
overlap <- merge(covariates, pca_loadings, by="row.names")
all_covariates <-cbind(overlap[11:14], overlap[18:102])
chipid <- as.vector(overlap$BeadChip)
status <- overlap$newstatus

correlations <- cor(all_covariates)

boxplot(correlations, main="All Correlations", las=2)
boxplot(correlations[,1:20], main="All Correlations, zoomed", las=2)

## Use cpgassoc to check for association between the beta values and the independent variable (ADHD),
#base model
cpgassoc_base<-cpg.assoc(betas(mldata.dasen.pf),indep=as.factor(status), chip.id=chipid)
## with PC1 as covariate in this case
cpgassoc_PC1<-cpg.assoc(betas(mldata.dasen.pf),indep=as.factor(status),covariates=all_covariates$PC1, chip.id=chipid, random=TRUE)
plot(cpgassoc_PC1)

### Some stats:
# Here's a plot of the percent of variance explained by each PC
# The red line shows a cutoff at 1.2 percent, the expected contribution
# of each PC if they were all equal.
sd <- pca_samples$sd
var<-sd^2
var.percent <- var/sum(var) * 100
barplot(var.percent, xlab="PC", ylab="Percent Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="gray")
abline(h=1/84*100, col="red")
#Here are the names of samples that contributed strongly to PC1:
row.names(pca_loadings)[which(pca_loadings$PC1>sqrt(1/84))]
