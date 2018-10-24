rm(list=ls());
library(RefFreeEWAS);
library(doParallel); registerDoParallel(detectCores() - 1)
Sample_sheet <- read.csv("~/BrockLab/OwnRes/Data/Oct16SampleSheetV3.csv")
#Sample_sheet <- read.csv("F:\\Fall2018\\brock\\OwnRes\\Datas\\Oct16SampleSheetV3.csv")
#betas <- readRDS("~/BrockLab/OwnRes/Data/NewPlots/OCT22Betas")

#betas <- readRDS("F:\\Fall2018\\brock\\OwnRes\\Datas\\betaWOSex.Rds")

load("~/BrockLab/OwnRes/Data/NewPlots/OCT22Betas.Rdata")

betas <- betas_OCT22
rm(betas_OCT22)

low <- c("stage i","stage ia","stage ii","stage iia","stage iib","stage iic")
high<- c("stage iii","stage iiia","stage iiib","stage iiic","stage iv","stage iva","stage ivb")
stable <- c("mss")
unstable <- c("msi")
subtypes <- c("CMS1","CMS2","CMS3","CMS4","NOLBL")


SubsetStage<- function(data, stage, stability,subtypes){
  if(stability == "stage"){
  tmpframe <- c()
  for(i in 1: length(stage)){
    tmp <- data[data$stage == stage[i],]
    tmp <- tmp[tmp$tissue.definition == "Primary solid Tumor",]
    tmp <- tmp[tmp$CMS_CBD == subtypes] 
    tmpframe <- rbind(tmpframe,tmp)
  }
  tmpframe <- rbind(tmpframe,data[data$tissue.definition == "Solid Tissue Normal",])
  return(tmpframe)
  }
  else if (stability == "stab"){
    tmpframe <- c()
    for(i in 1: length(stage)){
      tmp <- data[data$MSI == stage[i],]
      tmp <- tmp[tmp$tissue.definition == "Primary solid Tumor",]
      tmpframe <- rbind(tmpframe,tmp)
    }
    tmpframe <- rbind(tmpframe,data[data$tissue.definition == "Solid Tissue Normal",])
    return(tmpframe)
  }
}



stageCovL1 <- SubsetStage(Sample_sheet, low,"stage","CMS1")
stageCovL2 <- SubsetStage(Sample_sheet, low,"stage","CMS2")
stageCovL3 <- SubsetStage(Sample_sheet, low,"stage","CMS3")
stageCovL4 <- SubsetStage(Sample_sheet, low,"stage","CMS4")
stageCovLU <- SubsetStage(Sample_sheet, low,"stage","NOLBL")

stageCovH1 <- SubsetStage(Sample_sheet, high,"stage","CMS1")
stageCovH2 <- SubsetStage(Sample_sheet, high,"stage","CMS2")
stageCovH3 <- SubsetStage(Sample_sheet, high,"stage","CMS3")
stageCovH4 <- SubsetStage(Sample_sheet, high,"stage","CMS4")
stageCovHU <- SubsetStage(Sample_sheet, high,"stage","NOLBL")

stageCovStable <- SubsetStage(Sample_sheet,stable,"stab")
stageCovMSI <- SubsetStage(Sample_sheet,unstable,"stab")

#n <- nrow(stageCov[stageCov$tissue.definition == "Primary solid Tumor", ])

my.list <- list(stageCovL1,stageCovL2,stageCovL3,stageCovL4,stageCovLU,stageCovH1,stageCovH2,stageCovH3,stageCovH4,stageCovHU)


for(i in 1:length(my.list)){

stageCov <- my.list[[i]]

newBeta <- betas[,stageCov$Sample_Name]

newBeta <- data.matrix(newBeta)

DMGR_Var <- apply(newBeta,1,var)

rankvar <- rank(-DMGR_Var)

Y_shortened <- newBeta[rankvar<=10000,]

DMGR_RefFree_Array <- RefFreeCellMixArray(Y_shortened, Klist=2:10, iters=10)

Y_shortened2<-Y_shortened

sapply(DMGR_RefFree_Array,deviance,Y=Y_shortened)

DMGR_RefFree_Array2 <- RefFreeCellMixArray(newBeta, Klist=2:10, iters=10, Yfinal=Y_shortened2)

sapply(DMGR_RefFree_Array2,deviance,Y=Y_shortened2)

celprop<-DMGR_RefFree_Array2[[2]]$Omega

RefFree_DMGR_Boots = RefFreeCellMixArrayDevianceBoots(DMGR_RefFree_Array2, Y_shortened2, R=200, bootstrapIterations=10)

RefFree_DMGR_Boots

#RefFree_DMGR_Boots2 <- RefFree_DMGR_Boots

apply(RefFree_DMGR_Boots[-1,],2,mean,trim=0.25)

which.min(apply(RefFree_DMGR_Boots[-1,],2,mean,trim=0.25))

class(RefFree_DMGR_Boots)
# Save Results
save(DMGR_RefFree_Array, DMGR_RefFree_Array2, RefFree_DMGR_Boots,
     file=paste("~/BrockLab/OwnRes/Data/NewPlots/", i, "_", "RefFree2.RData", sep = ""), 
     compress=TRUE)
}
