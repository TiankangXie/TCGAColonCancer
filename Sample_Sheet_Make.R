library(GenomicDataCommons)
library(magrittr)
library(data.table)
library(readr)

#IMPORTANT!!!! 
#New method introduced here (primary method)
#--------------------------------------
#Deal with duplicated values
#---------------------------------
library(TCGAbiolinks)
#chekcing FFPE value
#FFPE <- read.csv("FFPEVals.csv")
#which(manifest$id %in% FFPE$UUID)


library(DT)


Garb1 <- read.table("F:\\Fall2018\\brock\\OwnRes\\CodesToUpload\\slide.tsv",sep = "\t",stringsAsFactors = F, header = T)

query <- GDCquery(project = c("TCGA-COAD"),
                  data.category = "Raw microarray data",
                  access = "open",
                  platform = "Illumina Human Methylation 450",
                  legacy = T)



query2 <- GDCquery(project = c("TCGA-READ"),
                   data.category = "Raw microarray data",
                   access = "open",
                   platform = "Illumina Human Methylation 450",
                   legacy = T)

# Only first 100 to make render faster

##Getting Sample sheet
tableResult <- getResults(query)
tableResult2 <- getResults(query2)

tableResult11 <- tableResult[-c(which(tableResult$data_type == "Normalized intensities")),]
tableResult21 <- tableResult2[-c(which(tableResult2$data_type == "Normalized intensities")),]

tableResult <- rbind(tableResult11,tableResult21)
rm(query,query2,tableResult11,tableResult21,tableResult2)

#Delete useless columns

deleting<- c("data_release","data_type","tags","submitter_id","file_size","state_comment","id","md5sum",
   "updated_datetime","data_format", "access" ,"platform", "version","data_category","type","experimental_strategy",
   "code","center_name","center_short_name","center_center_id","center_namespace","center_center_type" )


#delete = c(1,2,3,5,7,9,12,13,14,15,17,18,19,20,23,24,25,26,27)
tableResult = tableResult[,-which(colnames(tableResult)%in%deleting)]


#change some col names
colnames(tableResult)[1]<- "filename"

colnames(tableResult)[3] <- "FullSampleName"


#Sample_Name = c()
Sentrix_ID = c()
Sentrix_Position = c()
Sample_Name = c()
#colnames(RealSampleSheet)[2]
#x <- Intermediate$Sample_ID[1]
#Save a copy
firstMerge <- tableResult
#tableResult <- firstMerge
#firstMerge <- FinalSampSheet
for(i in 1:length(firstMerge$filename)){
  x <- tools::file_path_sans_ext(firstMerge$filename[i])
  NumPos <- gregexpr(pattern ='_',x)[[1]][2]
  pos <- strsplit(x,"_")
  Sample_Name[i] <- substring(x, 1,NumPos-1)
  Sentrix_ID[i] <- pos[[1]][1]
  Sentrix_Position[i] <- pos[[1]][2]
}
#Attach some necessary columns to the sample sheet
tableResult$Sentrix_ID <- Sentrix_ID
tableResult$Sentrix_Position <- Sentrix_Position
tableResult$Sample_Name <- Sample_Name
#FinalSampSheet$Sample_Name <- Sample_Name


Well_Position <- c()
New_SampleID <- c()

for(i in 1:length(tableResult$FullSampleName)){
  x <- as.character(tableResult$FullSampleName[i])
  New_SampleID[i] <- substring(x, 1,gregexpr(pattern ='-',x)[[1]][3]-1)
  Well_Position[i] <- substring(x, gregexpr(pattern ='-',x)[[1]][3]+1, gregexpr(pattern ='-',x)[[1]][4]-1)
}

tableResult$New_SampleID <- New_SampleID
tableResult$Well_Position <- Well_Position

#Download clinical data from GDC
clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
clinical2  <- GDCquery_clinic(project = "TCGA-READ", type = "clinical")

clinical2$submitter_id %in% clinical$submitter_id

clinical <- rbind(clinical,clinical2)




Hate = integer(0)
tableResult$gender <- NA
tableResult$age <- NA
tableResult$stage <- NA
tableResult$conciseSite <- NA


for(i in 1:nrow(tableResult)){
  indx <-  which(clinical$submitter_id == tableResult$New_SampleID[i])
  if(!identical(Hate,indx)){
    tableResult$gender[i] <- clinical$gender[indx]
    tableResult$age[i] <- as.numeric(clinical$age_at_diagnosis[indx]) %/% 365
    tableResult$stage[i] <- clinical$tumor_stage[indx]
    tableResult$conciseSite[i] <- clinical$site_of_resection_or_biopsy[indx]
  }
}



delete2 <- c(2)
SampleSheet <- tableResult[,-delete2] 
SampleSheet2 <- SampleSheet[-which(is.na(SampleSheet$gender)),]
SampleSheet3 <- SampleSheet2[-which(SampleSheet2$stage == "not reported"),]
SampleSheet <- SampleSheet3


#Merge Known Classification methods
#---------------------------------------------
NewSampleCBDSheet <- read.table("F:\\Fall2018\\brock\\OwnRes\\CodesToUpload\\ClinicalMode.txt",sep = "\t",stringsAsFactors = F, header = T)
NewSampleCBDSheet1 <- NewSampleCBDSheet[which(NewSampleCBDSheet$dataset == "tcga"),]



SampleSheet$cimp <- NA
SampleSheet$CMS_CBD <-NA
SampleSheet$MSI <- NA


for(i in 1:nrow(SampleSheet)){
  if(SampleSheet$New_SampleID[i] %in% NewSampleCBDSheet$sample){
    indx <-  which(NewSampleCBDSheet$sample == SampleSheet$New_SampleID[i])
    SampleSheet$cimp[i] <- NewSampleCBDSheet$cimp[indx]
    SampleSheet$MSI[i] <- NewSampleCBDSheet$msi[indx]
    SampleSheet$CMS_CBD[i] <- NewSampleCBDSheet$cms_label[indx]
  }else{
    SampleSheet$cimp[i] <- NA
    SampleSheet$MSI[i] <- NA
    SampleSheet$CMS_CBD[i] <- NA
  }
}
SampleSheet2 <- SampleSheet[-which(is.na(SampleSheet$CMS_CBD)),]



FinalSample <- SampleSheet2


#delete useless rows without gender/age/stage
#--------------------------

#12 to delete samples


Duptab <- FinalSample[which(duplicated(FinalSample$Sample_Name) == T),]

FinalSampSheet <- Duptab

FinalSampSheet$sex <- NA

for(i in 1:nrow(FinalSampSheet)){
  if(FinalSampSheet$gender[i] == "male"){
    FinalSampSheet$sex[i] <- "M"
  }
  else if (FinalSampSheet$gender[i] == "female"){
    FinalSampSheet$sex[i] <- "F"
  }
}

FinalSampSheet = FinalSampSheet[,-11]

write.csv(FinalSampSheet,"F:\\Fall2018\\brock\\OwnRes\\CodesToUpload\\Oct16SampleSheetV3.csv",row.names = F)



