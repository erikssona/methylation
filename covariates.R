
source("http://bioconductor.org/biocLite.R")
biocLite("pcaMethods")
library(pcaMethods)## Data cleaning/preprocessing before conducting PCA on the methylation dataset
library("CpGassoc")

# Get the covariates
CovariateFileName <- "~/methydata/covariates.txt"
covariates<-read.delim(CovariateFileName, header=TRUE, stringsAsFactors=TRUE)
#add row names for merging
rownames(covariates)<-covariates$Comments

#add pca_loadings (betas) to pheno matrix
#first figured out which record was causing a problem
#by having some kind of null value and screwing up the calculation
#Went back to watermethy script and removed that sample
#resulting in 84 samples for PCA and covariate analysis
row.names(pData(mldata.pf))<-pData(mldata.pf)$SampleID
overlap <- na.omit(merge(covariates,pData(mldata.pf)[1],by="row.names"))
overlap_no <- merge(covariates,pData(mldata.pf),by="row.names")
outcast <- setdiff(overlap_no.r[[1]], overlap.r[[1]])

#######
# This was how I took care of the problem sample (must do this with all MethyLumiSet objects before/after normalizing)
pData(mldata.pf)<-pData(mldata.pf)[pData(mldata.pf)[1]!=outcast,]
betas(mldata.pf)<-betas(mldata.pf)[,colnames(betas(mldata.pf))!=outcast]

## Do PCA on betas before normalizing
pca_samples.before<-prcomp(na.omit(betas(mldata.pf)),center=TRUE,scale=TRUE)
pca_loadings.before<-as.data.frame(pca_samples.before$rotation)
overlap.before <- merge(covariates, pca_loadings.before, by="row.names")
polymatrix.before <-cbind(overlap.before[10:14], overlap.before[16:102])
chipid.before <- overlap.before[5][[1]]

## Use cpgassoc to check for association between the beta values and the independent variable (ADHD),
## with PC1 as covariate in this case
cpgassoc_PC1.before<-cpg.assoc(betas(mldata.pf),indep=as.factor(polymatrix.before$newstatus),covariates=polymatrix.before$PC1, chip.id=chipid.before)



##Try the same thing all over again with nasen:

pca_samples<-prcomp(betas(mldata.nasen.pf),center=TRUE,scale=TRUE)
pca_loadings<-as.data.frame(pca_samples$rotation)

#add row names for merging
rownames(covariates)<-covariates$Comments

#add pca_loadings (betas) to pheno matrix
##FIRST figure out which record is causing the problem by having some kind of null value and screwing up the calculation
overlap <- na.omit(merge(covariates,pca_loadings,by="row.names"))
overlap_no <- merge(covariates,pca_loadings,by="row.names")
outcast <- setdiff(overlap_no[[1]], overlap[[1]])


### Before removing, try to figure out what's wrong with the sample outcast
which(betas(mldata.nasen.pf)[,colnames(betas(mldata.nasen.pf))==outcast])
which(!is.finite(as.matrix(pca_loadings)))
pData(mldata.dasen.pf)<-pData(mldata.dasen.pf)[pData(mldata.dasen.pf)[1]!=outcast,]
## Ha finally figured it out. The covariates were missing Stim_med column for that sample.


