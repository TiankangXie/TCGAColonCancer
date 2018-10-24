################################
# Load Libraries
################################

################################
# Analysis
################################
#s <- c('low', 'high')
library(limma)
library(matrixStats)
library(doParallel); registerDoParallel(detectCores() - 1)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotCpGs <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

#t <- c('Basal', 'Her2', 'LumA', 'LumB', 'Normal')
low <- c("stage i","stage ia","stage ii","stage iia","stage iib","stage iic")
high <- c("stage iii","stage iiia","stage iiib","stage iiic","stage iv","stage iva","stage ivb")
load("F:\\Fall2018\\brock\\OwnRes\\CodesToUpload\\OCT22Betas.Rdata");
masterCovar <- read.csv("F:\\Fall2018\\brock\\OwnRes\\Datas\\Oct16SampleSheetV3.csv")
masterCovar$short <- NA
masterCovar$normal <- "Non"
#masterCovar$CMS_CBD <- as.character(masterCovar$CMS_CBD)

for(i in 1:nrow(masterCovar)){
  if(masterCovar$tissue.definition[i] == "Solid Tissue Normal"){
    masterCovar$normal[i] <- "normal"
  }
  if(masterCovar$stage[i] %in% low){
    masterCovar$short[i] <- "low"
  }
  else if (masterCovar$stage[i] %in% high){
    masterCovar$short[i] <- "high"
  }
}
s <- c('low', 'high')
t <- c("CMS1","CMS2","CMS3","CMS4","NOLBL")
#--------------------------------------Load data & switch out barcode with sample IDs--------------------------------------
## Annotation package:\

masterCovar <- masterCovar[!is.na(masterCovar$age),]

s<-c("high")
t <- c("NOLBL")
for(i in 1:length(s)){
  for(j in 1:length(t)){
    
    listofsamples <- as.character(masterCovar$Sample_Name[which(masterCovar$short == s[i] & masterCovar$CMS_CBD == t[j])])
    subBetas <- betas_OCT22[,listofsamples]
    sheet <- masterCovar[which(masterCovar$short == s[i] & masterCovar$CMS_CBD == t[j]),]
    
    k <- 1e4 #can choose a different number of CpGs
    sele <- order(rowVars(subBetas), decreasing=TRUE)[1:k]
    subBetas <- subBetas[sele, ]
    ## Transform to M-values:
    mVals <- minfi::logit2(subBetas) #can also use: boot::logit(betas)
    
    #----------------------------------Step 2. Assemble Design Matrix ----------------------------------
    ## Important checkpoint: Double check myDesign rows match with methylation data columns:
    stopifnot(identical(as.character(sheet$Sample_Name), colnames(mVals)))
    
    ## Assemble the design matrix:
    sheet$Status <- ifelse(sheet$normal=="Non", yes=1, no=0) #alternatively: sheet$Status <- 1*(sheet$Sample_Group=="Cases")
    myDesign <- model.matrix( ~ Status + sex + age, data=sheet) #include more/fewer covariates as needed
    myDesign #Make sure that column no.2 is named "Status1", not "Status0"
    
    #----------------------------------Step 3. LIMMA Comparisons----------------------------------
    ## Execute statistical tests:
    fit <- lmFit(mVals, design=myDesign)
    fit2 <- eBayes(fit)
    
    ## Subset 450K annotation & identify:
    ## Load Illumina annotation files via Bioconductor data packaage:

    ## Compile a list of CpGs with Illumina annotation:
    myCpGs <- annotCpGs[match(rownames(mVals),annotCpGs$Name), ]
    DMPs <- topTable(
      fit2,
      number = Inf,
      coef = "Status",
      genelist = myCpGs, 
      adjust.method = "fdr",
      sort.by = "p"
    )
    
    write.csv(DMPs, paste("F:\\Fall2018\\brock\\OwnRes\\CodesToUpload\\",s[i],t[j],".csv",sep = ""), row.names=FALSE, quote=FALSE)
    print("success")
    }
}
## Methylation data:
## Master covariate file:

#--------------------------------------Select most variable CpGs-------------------------------------