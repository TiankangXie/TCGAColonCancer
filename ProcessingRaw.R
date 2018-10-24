#setwd("..")
#load("PhaseI.RData")
##############################################################################################################
# minfi pre-processing and QC
# Script author: David Chen
# Date: 09/07/2017
# Notes:
###############################################################################################################
library(minfi)
library(matrixStats)
library(doParallel); registerDoParallel(detectCores() - 1)

#---------------------------------------IDAT file loading & outlier classification---------------------------------------

target_mo <- read.metharray.sheet(base = "~/BrockLab/OwnRes/Data",pattern="Oct16SampleSheetV3.csv");

mo_obj <- read.metharray.exp(targets=target_mo);

mo_obj@annotation[2] <- "ilmn12.hg19"; #update annotation; may not be required at a later time

## Extract sample ID and group label for plotting:

subjects <- pData(mo_obj)$Sample_Name;
groups <- pData(mo_obj)$tissue.definition;



somePDFPath = "~/BrockLab/OwnRes/Data/NewPlots/DiagnosisPlot.pdf"
pdf(file=somePDFPath)  
## Density plots:
densityPlot(mo_obj,sampGroups=groups,main="Raw intensities",cex=0.75);
densityBeanPlot(mo_obj, sampNames=subjects, sampGroups=groups, main="Raw intensities by sample");

## Outlier plot: Convert to a MethylSet
Mset <- preprocessRaw(mo_obj);
Mset <- minfiQC(Mset, fixOutliers=TRUE, verbose=TRUE);
plotQC(Mset$qc);
title("Poor-performing outlier identification"); 

## Bisulfite conversion check:
controlStripPlot(mo_obj, controls="BISULFITE CONVERSION I");
controlStripPlot(mo_obj, controls="BISULFITE CONVERSION II");

dev.off()
## Save initial workspace as backup: 
save(
  list = c("mo_obj","target_mo","Mset"),
  file = "~/BrockLab/OwnRes/Data/NewPlots/RawFiles.Rdata",
  compress = TRUE
);

#---------------------------------------Quality control & normalization---------------------------------------
## Step 1. Perform normalization (Funnorm) and background correction (methylumi.noob)
normalizedEPIC <- preprocessFunnorm(mo_obj);
dim(normalizedEPIC)

## Step 2. Remove probes failed to meet detection P-value threshold of 0.05 in >20% samples:
pvals <- detectionP(mo_obj);
failedP <- (pvals > 0.01); #Be careful with the inequality sign!
mean(failedP) #proportion of failed


fraction <- 0.20;
failedProbes <- rownames(failedP)[rowMeans(failedP) > fraction]; ##list of probes
sum(rowMeans(failedP) > fraction); 
mean(rowMeans(failedP) > fraction);

normalizedEPICfilter <- normalizedEPIC[! rownames(normalizedEPIC) %in% failedProbes];
dim(normalizedEPICfilter)

## Step 3. Remove non-CpGs, control SNP probes, and polymorphic SNP probes:
normEPIC.NoSNPs <- dropMethylationLoci(normalizedEPICfilter); #drop technical SNP probes
normEPIC.NoSNPs <- dropLociWithSnps(normEPIC.NoSNPs); #drop polymorphic SNPs; MAF set to 0 by default

dim(normEPIC.NoSNPs)

## Extract beta values:
betas_OCT22 <- getBeta(normEPIC.NoSNPs);

## Step 4. Remove sex probes (final step):
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot.850kb3 <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19);

sexProbes <- as.character(annot.850kb3$Name[annot.850kb3$chr %in% c("chrX", "chrY")]); #probes to REMOVE to get autosomal only

betas_OCT22 <- betas_OCT22[! rownames(betas_OCT22) %in% sexProbes, ]; 


#---------------------------------------Export beta-values with annnotations---------------------------------------
## Visualize densities of the final data set:
png("~/BrockLab/OwnRes/Data/NewPlots/FigureS1A.png", res=300, units="in", height=8.27, width=11.69);
densityPlot(betas_OCT22, sampGroups=groups, cex.lab=1.5);
dev.off();

png("~/BrockLab/OwnRes/Data/NewPlots/FigureS1B.png", res=300, units="in", height=8.27, width=11.69);
densityBeanPlot(betas_OCT22, sampNames=subjects, sampGroups=groups);
dev.off()

## Export (one-time only): 
save(
  list = c("betas_OCT22","target_mo"),
  file = "~/BrockLab/OwnRes/Data/NewPlots/OCT22Betas.Rdata",
  compress = TRUE
)