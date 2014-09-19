#######
## Is methylation (measured by betas) associated with
## phenotype status (ADHD/ctrl)?
## To answer this we use cpg.assoc to test the basic model
## methylation ~ phenostatus, then other models that include
## various covariates.
## One covariate is the principal component of PCA on the betas.
############################
load("~/JustMethyNorm.RData")

library(methylumi) # needed for getting betas from methyLumiSet objects
library(pcaMethods) # Data cleaning/preprocessing before conducting PCA on the methylation dataset
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
all_covariates <-cbind(overlap[10:14], overlap[16:102])
# chipid and status vectors
chipid <- as.vector(overlap$BeadChip)
status <- overlap$newstatus

## Note there are 41 status 0 and 43 status 1'a
#The following plot looks at the question, is there a difference between the relationship between PC1 and PC2 based on status?
plot(pca_loadings$PC1[which(status==1)]/sum(pca_loadings$PC1), pca_loadings$PC2[which(status==1)]/sum(pca_loadings$PC2), col="blue")
points(pca_loadings$PC1[which(status==0)]/sum(pca_loadings$PC1), pca_loadings$PC2[which(status==0)]/sum(pca_loadings$PC2), col="red")
#And this plot explores whether there is a difference between case and control and the percentage of impact on PC1
plot(pca_loadings$PC1[which(status==0)]/sum(pca_loadings$PC1), col="red")
points(pca_loadings$PC1[which(status==1)]/sum(pca_loadings$PC1), col="blue")

#look at covariates
plot(all_covariates$Wave) #factor
plot(all_covariates$Stim_med) #factor
plot(all_covariates$MINORITY_STATUS) #factor
plot(all_covariates$FSIQ.y)
plot(all_covariates$Age_Months_byWave)
plot(all_covariates$newstatus) #factor
plot(all_covariates$newsubpheno) #factor
plot(all_covariates$rowcolumn2) #factor
plot(all_covariates$PC1)

summary(all_covariates$Wave) #factor
summary(all_covariates$Stim_med) #factor
summary(all_covariates$MINORITY_STATUS) #factor
summary(all_covariates$FSIQ.y)
summary(all_covariates$Age_Months_byWave)
summary(all_covariates$newstatus) #factor
summary(all_covariates$newsubpheno) #factor
summary(all_covariates$rowcolumn2) #factor
summary(all_covariates$PC1)

#Look at correlations
correlations_before <- cor(pca_loadings)
correlations <- cor(all_covariates)
correlations_before
correlations

## Use cpgassoc to check for association between the beta values and the independent variable (ADHD),
#base model
cpgassoc_base<-cpg.assoc(beta.val=betas(mldata.dasen),indep=as.factor(status), chip.id=chipid)
## with PC1 as covariate in this case
cpgassoc_PC1<-cpg.assoc(beta.val=betas(mldata.dasen),indep=as.factor(status),covariates=all_covariates$PC1, chip.id=chipid)
## with age as a covariate
cpgassoc_age<-cpg.assoc(beta.val=betas(mldata.dasen),indep=as.factor(status),covariates=all_covariates$Age_Months_byWave, chip.id=chipid)
## with both age and PC1 as covariates
cpgassoc_all<-cpg.assoc(beta.val=betas(mldata.dasen),indep=as.factor(status),covariates=all_covariates$PC1+all_covariates$Age_Months_byWave, chip.id=chipid)

###Also, I ran CpGassoc with meds, rowcolum, minority_status, chip as covariates as well.)
cpgassoc_meds<-cpg.assoc(beta.val=betas(mldata.dasen),indep=as.factor(status),covariates=as.factor(all_covariates$Stim_med), chip.id=chipid)
cpgassoc_rowcol<-cpg.assoc(beta.val=betas(mldata.dasen),indep=as.factor(status),covariates=as.factor(all_covariates$rowcolumn2), chip.id=chipid)
cpgassoc_minority<-cpg.assoc(beta.val=betas(mldata.dasen),indep=as.factor(status),covariates=as.factor(all_covariates$MINORITY_STATUS), chip.id=chipid)
cpgassoc_chip<-cpg.assoc(beta.val=betas(mldata.dasen),indep=as.factor(status), covariates=as.factor(chipid), chip.id=as.factor(chipid))

#plots
par(mfrow=c(2,2))
plot(cpgassoc_base, main="QQ plot for association between \nmethylation and status (base model)")
plot(cpgassoc_PC1, main="QQ plot for association between \nmethylation and status (PC1 as covariate)")
plot(cpgassoc_age, main="QQ plot for association between \nmethylation and status (age as covariate)")
plot(cpgassoc_all, main="QQ plot for association between \nmethylation and status (PC1 and age as covariates)")

plot(cpgassoc_meds, main="QQ plot of association between \nmethylation and status (meds as covariate)")
plot(cpgassoc_rowcol, main="QQ plot of association between \nmethylation and status (row/column as covariate)")
plot(cpgassoc_minority, main="QQ plot of association between \nmethylation and status (minority status as covariate)")
plot(cpgassoc_chip, main="QQ plot of association between \nmethylation and status (beadchip as covariate)")


#annotations for manhattan plots
annotation<-pData(featureData(mldata.dasen))

manhattan(cpgassoc_base,annotation$TargetID,annotation$CHR,annotation$MAPINFO, main="Association between methylation \nand status (base model)")
manhattan(cpgassoc_PC1,annotation$TargetID,annotation$CHR,annotation$MAPINFO, main="Association between methylation \nand status (PC1 as covariate)")
manhattan(cpgassoc_age,annotation$TargetID,annotation$CHR,annotation$MAPINFO, main="Association between methylation \nand status (age as covariate)")
manhattan(cpgassoc_all,annotation$TargetID,annotation$CHR,annotation$MAPINFO, main="Association between methylation \nand status (PC1 and age as covariates)")

manhattan(cpgassoc_meds,annotation$TargetID,annotation$CHR,annotation$MAPINFO, main="Association between methylation \nand status (meds as covariate)")
manhattan(cpgassoc_rowcol,annotation$TargetID,annotation$CHR,annotation$MAPINFO, main="Association between methylation \nand status (row/column as covariate)")
manhattan(cpgassoc_minority,annotation$TargetID,annotation$CHR,annotation$MAPINFO, main="Association between methylation \nand status (minority status as covariate)")
manhattan(cpgassoc_chip,annotation$TargetID,annotation$CHR,annotation$MAPINFO, main="Association between methylation \nand status (chip as covariate)")
