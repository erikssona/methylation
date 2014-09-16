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

## Do PCA on betas
pca<-prcomp(betas(mldata.dasen),center=TRUE,scale=TRUE)
# get the loadings or eigenvalues
pca_loadings<-as.data.frame(pca$rotation)
# get the standard deviations, variances and percent variance of the PCs
sd <- pca$sd
var<-sd^2
var.percent <- var/sum(var) * 100

### Some stats:
# Here's a plot of the percent of variance explained by each PC
# The red line shows a cutoff at 1.2 percent, the expected contribution
# of each PC if they were all equal.
barplot(var.percent, xlab="PC", ylab="Percent Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="gray")
abline(h=1/84*100, col="red")
#Here are the names of samples that contributed strongly to PC1:
#(strongly defined as more contribution than if all were equal contributors)
row.names(pca_loadings)[which(pca_loadings$PC1>sqrt(1/84))]

# get the samples for which we have covariates & add PCs
overlap <- merge(covariates, pca_loadings, by="row.names")
# isolate the numerical covariates we want to use
all_covariates <-cbind(overlap[11:14], overlap[18:102])
# chipid and status vectors
chipid <- as.vector(overlap$BeadChip)
status <- overlap$newstatus

#Look at correlations
correlations <- cor(all_covariates)

## Use cpgassoc to check for association between the beta values and the independent variable (ADHD),
#base model
cpgassoc_base<-cpg.assoc(betas(mldata.dasen),indep=as.factor(status), chip.id=chipid)
## with PC1 as covariate in this case
cpgassoc_PC1<-cpg.assoc(betas(mldata.dasen),indep=as.factor(status),covariates=all_covariates$PC1, chip.id=chipid, random=TRUE)

#plots
par(mfrow=c(2,2))
plot(cpgassoc_base)
plot(cpgassoc_PC1)
