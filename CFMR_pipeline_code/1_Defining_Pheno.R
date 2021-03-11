rm(list=ls());gc(TRUE)
options(stringsAsFactors = F)
set.seed(1)
require(feather);require(data.table) 
load("~/archive/START/pheno/PDB_2374_MoBa_Q1_Q4_MBRN_full.RData")
info=MoBa_Q1_Q4_MBRN_full
for(i in c("mother","father","child"))
{
  print(i)
  key_one=fread(paste("~/archive/START/pheno/2019_10_01_MoBaGenetics_",i,"_2374.csv",sep=""),
                data.table=F)  
  colnames(key_one)[2:7]=paste(substr(i,1,1),"_",colnames(key_one)[2:7],sep="")
  print(colnames(key_one)[2:7])
  print(table(key_one[,paste(substr(i,1,1),"_Role",sep="")]))
  info=merge(x=info,y=key_one,by=colnames(key_one)[1],sep="")
}
str(info)



momtemp <- read.csv("~/archive/START/pheno/2019_10_01_MoBaGenetics_mother_2374.csv", sep=";")
dadtemp <- read.csv("~/archive/START/pheno/2019_10_01_MoBaGenetics_father_2374.csv", sep=";")

mom_Sentrix_Preg_ID <-  merge ( momtemp, info, by = "M_ID_2374")





motherpheno <- mom_Sentrix_Preg_ID

library(dplyr)
#Some duplicate 
motherpheno <- motherpheno[ -which( duplicated(motherpheno$PREG_ID_2374)),]
#Removing mutiple pregnancy
temp <- motherpheno [ which(motherpheno$SENTRIXID %in% motherpheno$SENTRIXID[which( duplicated(motherpheno$SENTRIXID))]),]

motherpheno_out <- motherpheno[ - which(duplicated(motherpheno$SENTRIXID)), ]




motherpheno_out$BMI <- motherpheno_out$AA85/((motherpheno_out$AA87/100)^2)

#removing the duplicated individual
motherpheno_out$BMI[which(motherpheno_out$BMI>50)] <-NA
motherpheno_out$BMI[which(motherpheno_out$BMI<10)] <-NA

motherpheno_out <- motherpheno_out[, c("SENTRIXID","BMI")]

motherpheno_out$BMI[which(is.na(motherpheno_out$BMI ))] <- -9#plink  regression format


colnames(motherpheno_out)  <- c( "IID", "BMI")



head(motherpheno_out)
#Build pipeline for other trait 

write.table(motherpheno_out, file ="/mnt/work/william.denault/CFMR_Analysis/Data/Maternal_BMI.txt", quote=FALSE, row.names = FALSE)






#Performing sample splitting for 10 GWASs


set.seed(1)
ss <- sample(1:10,size=dim(motherpheno_out)[1],replace=TRUE,prob=rep(0.1,10))




motherpheno_split_1_out             <- motherpheno_out
motherpheno_split_1_out$BMI[which(ss==1) ] <- -9#removing from the GWAS the ind split 1, -9 correspond to NA in plink format
motherpheno_split_2_out             <- motherpheno_out
motherpheno_split_2_out$BMI[which(ss==2) ] <- -9 
motherpheno_split_3_out             <- motherpheno_out
motherpheno_split_3_out$BMI[which(ss==3) ] <- -9 
motherpheno_split_4_out             <- motherpheno_out
motherpheno_split_4_out$BMI[which(ss==4) ] <- -9 
motherpheno_split_5_out             <- motherpheno_out
motherpheno_split_5_out$BMI[which(ss==5) ] <- -9 
motherpheno_split_6_out             <- motherpheno_out
motherpheno_split_6_out$BMI[which(ss==6) ] <- -9 
motherpheno_split_7_out             <- motherpheno_out
motherpheno_split_7_out$BMI[which(ss==7) ] <- -9 
motherpheno_split_8_out             <- motherpheno_out
motherpheno_split_8_out$BMI[which(ss==8) ] <- -9 
motherpheno_split_9_out             <- motherpheno_out
motherpheno_split_9_out$BMI[which(ss==9) ] <- -9 
motherpheno_split_10_out             <- motherpheno_out
motherpheno_split_10_out$BMI[which(ss==10) ] <- -9 






write.table(motherpheno_split_1_out  , file ="/mnt/work/william.denault/CFMR_Analysis/Data/Split_1_Maternal_BMI.txt", quote=FALSE, row.names = FALSE)
write.table(motherpheno_split_2_out  , file ="/mnt/work/william.denault/CFMR_Analysis/Data/Split_2_Maternal_BMI.txt", quote=FALSE, row.names = FALSE)
write.table(motherpheno_split_3_out  , file ="/mnt/work/william.denault/CFMR_Analysis/Data/Split_3_Maternal_BMI.txt", quote=FALSE, row.names = FALSE)
write.table(motherpheno_split_4_out  , file ="/mnt/work/william.denault/CFMR_Analysis/Data/Split_4_Maternal_BMI.txt", quote=FALSE, row.names = FALSE)
write.table(motherpheno_split_5_out  , file ="/mnt/work/william.denault/CFMR_Analysis/Data/Split_5_Maternal_BMI.txt", quote=FALSE, row.names = FALSE)
write.table(motherpheno_split_6_out  , file ="/mnt/work/william.denault/CFMR_Analysis/Data/Split_6_Maternal_BMI.txt", quote=FALSE, row.names = FALSE)
write.table(motherpheno_split_7_out  , file ="/mnt/work/william.denault/CFMR_Analysis/Data/Split_7_Maternal_BMI.txt", quote=FALSE, row.names = FALSE)
write.table(motherpheno_split_8_out  , file ="/mnt/work/william.denault/CFMR_Analysis/Data/Split_8_Maternal_BMI.txt", quote=FALSE, row.names = FALSE)
write.table(motherpheno_split_9_out  , file ="/mnt/work/william.denault/CFMR_Analysis/Data/Split_9_Maternal_BMI.txt", quote=FALSE, row.names = FALSE)
write.table(motherpheno_split_10_out , file ="/mnt/work/william.denault/CFMR_Analysis/Data/Split_10_Maternal_BMI.txt", quote=FALSE, row.names = FALSE)

pheno_all <- cbind(motherpheno_split_1_out,
                   motherpheno_split_2_out[,2],
                   motherpheno_split_3_out[,2],
                   motherpheno_split_4_out[,2],
                   motherpheno_split_5_out[,2],
                   motherpheno_split_6_out[,2],
                   motherpheno_split_7_out[,2],
                   motherpheno_split_8_out[,2],
                   motherpheno_split_9_out[,2],
                   motherpheno_split_10_out[,2])

