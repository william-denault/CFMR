rm( list =ls())

library(glmnet)
library(hdm)
library(rms)
library(data.table)
library(SuperLearner)
library(xgboost)
library(glmnet)
library(ranger)
library(randomForest)


#Father and the 10 splits ----
library(data.table)

'%!in%' <- function(x,y)!('%in%'(x,y))
setwd(dir = "/mnt/work/william.denault/CFMR_Analysis/Data")
i=1
 
temp <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
print(dim(temp))
for ( i in 2:10)
{
   
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
    print(dim(tt))
  
  to_rm <- which( tt$ID %in% temp$ID)
  print(length(to_rm))
  if( length(to_rm)==0)
  {
    temp <- rbind(temp,tt)
  }  else{
    tt <- tt [-to_rm,]
    temp <- rbind(temp,tt)
  }
  
  
  
}
tt <-temp


#there is much more efficient for sure but this works
library(data.table)
library(dplyr)
gc(TRUE)
options(stringsAsFactors = F)

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
infom <-  merge ( momtemp, info[,c(1,3)], by = "M_ID_2374")



setwd("/mnt/work/william.denault/CFMR_Analysis/Data")
lf <- list.files()
lf <- lf[grep("^Extrac", lf)]
lf <- lf[grep(".vcf", lf)]

temp <- fread(lf[1])
var.names        <- temp$ID
temp             <- temp [ , - c( "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")]#removong all non SNP col
SENTRIXID        <- colnames(temp)
row.names(temp)  <- c()
temp             <- t(temp)#puting the matrix in line= in, row = SNP
colnames(temp)   <- var.names
temp2            <- temp[which(SENTRIXID %in% momtemp$SENTRIXID ),] #mother
temp2 <- cbind( temp2,SENTRIXID[which(SENTRIXID %in% momtemp$SENTRIXID )])
colnames(temp2)[dim(temp2)[2]] <- "SENTRIXID"
colnames(infom)
temp2 <- merge(temp2, infom[,c(1,2)], by="SENTRIXID")
temp2 <- temp2[,- dim(temp)[2]]

cleaning_file <- function(x)
{
  temp <- fread(lf[x])
  var.names        <- temp$ID
  temp             <- temp [ , - c( "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")]#removong all non SNP col
  SENTRIXID        <- colnames(temp)
  row.names(temp)  <- c()
  temp             <- t(temp)
  colnames(temp)   <- var.names
  temp2            <- temp[which(SENTRIXID %in% momtemp$SENTRIXID ),] #mother
  temp2 <- cbind( temp2,SENTRIXID[which(SENTRIXID %in% momtemp$SENTRIXID )])
  colnames(temp2)[dim(temp2)[2]] <- "SENTRIXID"
  colnames(infom)
  temp2 <- merge(temp2, infom[,c(1,2)], by="SENTRIXID")
  temp2 <- temp2[,- dim(temp)[2]]
  temp2
}



temp_l <-list()
temp_l[[1]] <- temp2
for ( i in 2: length(lf))
{
  temp_l[[i]] <- cleaning_file(i)
  print(i)
}

Reg_MAT <- do.call(cbind, temp_l)


Reg_MAT <- Reg_MAT[-which(duplicated(Reg_MAT[,1])),]
D_init <- read.csv("/mnt/work/william.denault/CFMR_Analysis/Data/Maternal_BMI.txt", sep="")

Reg_MAT <- Reg_MAT[which(Reg_MAT[,1]%in% D_init$IID),]
D_init <- D_init [ which( D_init$IID%in% Reg_MAT[,1]) ,]
Reg_MAT[1:10,1]
D_init[1:10,1]
D_init <- D_init[order(D_init$IID),]
Reg_MAT <- Reg_MAT[order(Reg_MAT[,1]),]
sum(D_init$IID==Reg_MAT[,1])


#Filtering on the SNPs selected in the GWAS

subdf <- as.matrix(Reg_MAT[,-1])

subdf[subdf=="0|0"] <-0
subdf[subdf=="0/0"] <-0
subdf[subdf=="0|1"] <-1
subdf[subdf=="1|0"] <-1
subdf[subdf=="0/1"] <-1
subdf[subdf=="1/0"] <-1
subdf[subdf=="1|1"] <-2
subdf[subdf=="1/1"] <-2

subdf <- subdf[,-grep( "SEN", colnames(subdf))]
subdf <- subdf[,-grep( "M_ID", colnames(subdf))]
tempcol <- colnames(subdf)
subdf <- t(apply(subdf, 1, as.numeric))
colnames(subdf ) <- tempcol
rownames(subdf) <- Reg_MAT[,1]




##Select split specific SNPs ----


##### PVthres 10-3----

set.seed(1)
options(mc.cores = 20)
#####Building IV1 ----
setwd(dir = paste("/mnt/work/william.denault/CFMR_Analysis/Data/", sep="") )
tt <- fread(paste(getwd(),"/Split_",1,"_SNP_for_extraction.txt",sep=""))
X <- subdf
 
colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
subtt <- tt[which(tt$P <10^(-3)),]
X <- X[ ,which(colnames(X) %in% subtt$ID)]




D_init$BMI[which(D_init$BMI ==-9)] <- NA
 

IV_all <- rep(NA,dim(D_init)[1])
#PREG ID in split1
Split_1_PREG_ID <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_1_Maternal_BMI.txt")
Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
#removing the one set as NA in the splitting procedure
sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
test_ID  <- (1:dim(D_init)[1])[-train_ID]

####Multiple NA in the observed data and glmnet does not handle them
y <- as.vector(D_init[train_ID,2])




x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
for( j in 1:dim( x)[2])
{
  x[,j] <- as.numeric(X[,j])
   
}
xtrain <- x[train_ID,]
xtest <- x[test_ID,]


Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1) 
save(Lasso_IV, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV1.RData")








#####Building IVs-----


#####Building the IVs ----
print("about to start buildin IV")
IVs_build <- function ( i)
  #for (i in 2:10)
{
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
   
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-3)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  sum( is.na(y))
   
   
  
   
  
  x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
  for( j in 1:dim( x)[2])
  {
    x[,j] <- as.numeric(X[,j])
     
  }
  xtrain <- x[train_ID,]
  xtest <- x[test_ID,]
  
  
  
  Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1)
  
  
  
   
  
   
  save(Lasso_IV , file = paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  print(paste("IV",i, "done"))
}

lapply(2:10, IVs_build)


tl <- list()
idx_list <-list()
mod_pred <- rep(NA, length(IV_all))
for( i in 1:10)
{
  
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
    
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-3)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  print(dim(X))
  
  length(which(!complete.cases(X)))
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  load(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  
  
  
  
  
  IV_all[test_ID] <- predict(Lasso_IV,X[test_ID,], s = "lambda.1se")
  tl[[i]] <- cbind(predict(Lasso_IV,X[train_ID,], s = "lambda.min"),y, rep(i,length(y)))
  
   
  idx_list[[i]] <- test_ID
   
  mod_pred[test_ID] <- rep( i, length(test_ID))
   
   
  Coef <- coef(Lasso_IV, s = "lambda.1se")
  print(length(Coef@i))
   
  
}


hist( IV_all)
IV_all[which(is.na(D_init$BMI))]  <- NA
plot( IV_all, D_init$BMI)

df_IV <- data.frame( IID= D_init$IID, IV=IV_all, BMI=D_init$BMI, GWAS=as.factor(mod_pred))
save(df_IV,file ="/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv3.RData")
ggplot()+
  geom_point(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),alpha=0.2 )+
  geom_smooth(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_IV,method="lm",se=FALSE,aes(x=IV, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the complementary data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_CFMRpv3.jpeg", device="jpeg",width = 4.35,height = 3.82)
ggplot(df_IV, aes(x = GWAS, y=(BMI-IV)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the complementary data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_box_CFMRpv3.jpeg", device="jpeg",width = 4.35,height = 3.82)
fit_test3 <- lm(D_init$BMI~IV_all)
summary(fit_test3)


df_train <- do.call(rbind, tl)
df_train <- data.frame(df_train)
colnames(df_train) <- c("Prediction","BMI", "GWAS")
df_train$GWAS <- as.factor(df_train$GWAS)


ggplot()+
  geom_point(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),alpha=0.1 )+
  geom_smooth(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_train,method="lm",se=FALSE,aes(x=Prediction, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the trainning data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")

ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_CFMRpv3.jpeg", device="jpeg",width = 4.35,height = 3.82)
ggplot(df_train, aes(x = GWAS, y=(BMI- Prediction)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the trainning data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")
ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_box_CFMRpv3.jpeg", device="jpeg",width = 4.35,height = 3.82)
fit_train3 <- lm(BMI~Prediction, df_train)
summary(fit_train3)


##### PVthres 10-4----

set.seed(1)
options(mc.cores = 20)
#####Building IV1 ----
setwd(dir = paste("/mnt/work/william.denault/CFMR_Analysis/Data/", sep="") )
tt <- fread(paste(getwd(),"/Split_",1,"_SNP_for_extraction.txt",sep=""))
X <- subdf
 
colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
subtt <- tt[which(tt$P <10^(-4)),]
X <- X[ ,which(colnames(X) %in% subtt$ID)]




D_init$BMI[which(D_init$BMI ==-9)] <- NA
 

IV_all <- rep(NA,dim(D_init)[1])
#PREG ID in split1
Split_1_PREG_ID <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_1_Maternal_BMI.txt")
Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
#removing the one set as NA in the splitting procedure
sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
test_ID  <- (1:dim(D_init)[1])[-train_ID]

####Multiple NA in the observed data and glmnet does not handle them
y <- as.vector(D_init[train_ID,2])




x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
for( j in 1:dim( x)[2])
{
  x[,j] <- as.numeric(X[,j])
   
}
xtrain <- x[train_ID,]
xtest <- x[test_ID,]


Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1) 
save(Lasso_IV, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV1.RData")








#####Building IVs-----


#####Building the IVs ----
print("about to start buildin IV")
IVs_build <- function ( i)
#for (i in 2:10)
{
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
   
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-4)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
y <- as.vector(D_init[train_ID,2])
sum( is.na(y))
 
 

 

x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
for( j in 1:dim( x)[2])
{
  x[,j] <- as.numeric(X[,j])
   
}
xtrain <- x[train_ID,]
xtest <- x[test_ID,]



Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1)

  
  
 

   
  save(Lasso_IV , file = paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  print(paste("IV",i, "done"))
}

lapply(2:10, IVs_build)


tl <- list()
idx_list <-list()
mod_pred <- rep(NA, length(IV_all))
for( i in 1:10)
{
  
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
    
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-4)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  print(dim(X))
  
  length(which(!complete.cases(X)))
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  load(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  
  

  
  
   IV_all[test_ID] <- predict(Lasso_IV,X[test_ID,], s = "lambda.1se")
   tl[[i]] <- cbind(predict(Lasso_IV,X[train_ID,], s = "lambda.min"),y, rep(i,length(y)))
   
    
  idx_list[[i]] <- test_ID
   
  mod_pred[test_ID] <- rep( i, length(test_ID))
  
   
  Coef <- coef(Lasso_IV, s = "lambda.1se")
  print(length(Coef@i))
   
  
}


hist( IV_all)
IV_all[which(is.na(D_init$BMI))]  <- NA
plot( IV_all, D_init$BMI)

df_IV <- data.frame( IID= D_init$IID, IV=IV_all, BMI=D_init$BMI, GWAS=as.factor(mod_pred))
save(df_IV,file ="/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv4.RData")
ggplot()+
  geom_point(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),alpha=0.2 )+
  geom_smooth(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_IV,method="lm",se=FALSE,aes(x=IV, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the complementary data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_CFMRpv4.jpeg", device="jpeg",width = 4.35,height = 3.82)
ggplot(df_IV, aes(x = GWAS, y=(BMI-IV)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the complementary data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_box_CFMRpv4.jpeg", device="jpeg",width = 4.35,height = 3.82)
fit_test4 <- lm(D_init$BMI~IV_all)
summary(fit_test4)


df_train <- do.call(rbind, tl)
df_train <- data.frame(df_train)
colnames(df_train) <- c("Prediction","BMI", "GWAS")
df_train$GWAS <- as.factor(df_train$GWAS)


ggplot()+
  geom_point(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),alpha=0.1 )+
  geom_smooth(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_train,method="lm",se=FALSE,aes(x=Prediction, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the trainning data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")

ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_CFMRpv4.jpeg", device="jpeg",width = 4.35,height = 3.82)
ggplot(df_train, aes(x = GWAS, y=(BMI- Prediction)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the trainning data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")
ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_box_CFMRpv4.jpeg", device="jpeg",width = 4.35,height = 3.82)
fit_train4 <- lm(BMI~Prediction, df_train)
summary(fit_train4)
##### PVthres 10-5-----------

set.seed(1)
options(mc.cores = 20)
#####Building IV1 ----
setwd(dir = paste("/mnt/work/william.denault/CFMR_Analysis/Data/", sep="") )
tt <- fread(paste(getwd(),"/Split_",1,"_SNP_for_extraction.txt",sep=""))
X <- subdf
 
colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
subtt <- tt[which(tt$P <10^(-5)),]
X <- X[ ,which(colnames(X) %in% subtt$ID)]




D_init$BMI[which(D_init$BMI ==-9)] <- NA
 

IV_all <- rep(NA,dim(D_init)[1])
#PREG ID in split1
Split_1_PREG_ID <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_1_Maternal_BMI.txt")
Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
#removing the one set as NA in the splitting procedure
sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
test_ID  <- (1:dim(D_init)[1])[-train_ID]

####Multiple NA in the observed data and glmnet does not handle them
y <- as.vector(D_init[train_ID,2])




x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
for( j in 1:dim( x)[2])
{
  x[,j] <- as.numeric(X[,j])
   
}
xtrain <- x[train_ID,]
xtest <- x[test_ID,]


Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1) 
save(Lasso_IV, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV1.RData")








#####Building IVs-----


#####Building the IVs ----
print("about to start buildin IV")
IVs_build <- function ( i)
  #for (i in 2:10)
{
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
   
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-5)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  sum( is.na(y))
   
   
  
   
  
  x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
  for( j in 1:dim( x)[2])
  {
    x[,j] <- as.numeric(X[,j])
     
  }
  xtrain <- x[train_ID,]
  xtest <- x[test_ID,]
  
  
  
  Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1)
  
  
  
   
  
   
  save(Lasso_IV , file = paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  print(paste("IV",i, "done"))
}

lapply(2:10, IVs_build)


tl <- list()
idx_list <-list()
mod_pred <- rep(NA, length(IV_all))
for( i in 1:10)
{
  
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
    
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-5)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  print(dim(X))
  
  length(which(!complete.cases(X)))
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  load(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  
  
  
  
  
  IV_all[test_ID] <- predict(Lasso_IV,X[test_ID,], s = "lambda.1se")
  tl[[i]] <- cbind(predict(Lasso_IV,X[train_ID,], s = "lambda.min"),y, rep(i,length(y)))
  
   
  idx_list[[i]] <- test_ID
   
  mod_pred[test_ID] <- rep( i, length(test_ID))
   
   
  Coef <- coef(Lasso_IV, s = "lambda.1se")
  print(length(Coef@i))
   
  
}


hist( IV_all)
IV_all[which(is.na(D_init$BMI))]  <- NA
plot( IV_all, D_init$BMI)

df_IV <- data.frame( IID= D_init$IID, IV=IV_all, BMI=D_init$BMI, GWAS=as.factor(mod_pred))
save(df_IV,file ="/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv5.RData")
ggplot()+
  geom_point(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),alpha=0.2 )+
  geom_smooth(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_IV,method="lm",se=FALSE,aes(x=IV, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the complementary data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_CFMRpv5.jpeg", device="jpeg",width = 4.35,height = 3.82)

ggplot(df_IV, aes(x = GWAS, y=(BMI-IV)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the complementary data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_box_CFMRpv5.jpeg", device="jpeg",width = 4.35,height = 3.82)
fit_test5 <- lm(D_init$BMI~IV_all)
summary(fit_test5)


df_train <- do.call(rbind, tl)
df_train <- data.frame(df_train)
colnames(df_train) <- c("Prediction","BMI", "GWAS")
df_train$GWAS <- as.factor(df_train$GWAS)


ggplot()+
  geom_point(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),alpha=0.1 )+
  geom_smooth(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_train,method="lm",se=FALSE,aes(x=Prediction, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the trainning data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")
ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_CFMRpv5.jpeg", device="jpeg",width = 4.35,height = 3.82)

ggplot(df_train, aes(x = GWAS, y=(BMI- Prediction)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the trainning data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")
ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_box_CFMRpv5.jpeg", device="jpeg",width = 4.35,height = 3.82)
fit_train5 <- lm(BMI~Prediction, df_train)
summary(fit_train5)

##### PVthres 10-6--------

set.seed(1)
options(mc.cores = 20)
#####Building IV1 ----
setwd(dir = paste("/mnt/work/william.denault/CFMR_Analysis/Data/", sep="") )
tt <- fread(paste(getwd(),"/Split_",1,"_SNP_for_extraction.txt",sep=""))
X <- subdf
 
colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
subtt <- tt[which(tt$P <10^(-6)),]
X <- X[ ,which(colnames(X) %in% subtt$ID)]




D_init$BMI[which(D_init$BMI ==-9)] <- NA
 

IV_all <- rep(NA,dim(D_init)[1])
#PREG ID in split1
Split_1_PREG_ID <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_1_Maternal_BMI.txt")
Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
#removing the one set as NA in the splitting procedure
sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
test_ID  <- (1:dim(D_init)[1])[-train_ID]

####Multiple NA in the observed data and glmnet does not handle them
y <- as.vector(D_init[train_ID,2])




x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
for( j in 1:dim( x)[2])
{
  x[,j] <- as.numeric(X[,j])
   
}
xtrain <- x[train_ID,]
xtest <- x[test_ID,]


Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1) 
save(Lasso_IV, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV1.RData")








#####Building IVs-----


#####Building the IVs ----
print("about to start buildin IV")
IVs_build <- function ( i)
  #for (i in 2:10)
{
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
   
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-6)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  sum( is.na(y))
   
   
  
   
  
  x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
  for( j in 1:dim( x)[2])
  {
    x[,j] <- as.numeric(X[,j])
     
  }
  xtrain <- x[train_ID,]
  xtest <- x[test_ID,]
  
  
  
  Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1)
  
  
  
   
  
   
  save(Lasso_IV , file = paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  print(paste("IV",i, "done"))
}

lapply(2:10, IVs_build)


tl <- list()
idx_list <-list()
mod_pred <- rep(NA, length(IV_all))
for( i in 1:10)
{
  
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
    
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-6)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  print(dim(X))
  
  length(which(!complete.cases(X)))
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  load(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  
  
  
  
  
  IV_all[test_ID] <- predict(Lasso_IV,X[test_ID,], s = "lambda.1se")
  tl[[i]] <- cbind(predict(Lasso_IV,X[train_ID,], s = "lambda.min"),y, rep(i,length(y)))
  
   
  idx_list[[i]] <- test_ID
   
  mod_pred[test_ID] <- rep( i, length(test_ID))
   
   
  Coef <- coef(Lasso_IV, s = "lambda.1se")
  print(length(Coef@i))
   
  
}


hist( IV_all)
IV_all[which(is.na(D_init$BMI))]  <- NA
plot( IV_all, D_init$BMI)

df_IV <- data.frame( IID= D_init$IID, IV=IV_all, BMI=D_init$BMI, GWAS=as.factor(mod_pred))
save(df_IV,file ="/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv6.RData")
ggplot()+
  geom_point(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),alpha=0.2 )+
  geom_smooth(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_IV,method="lm",se=FALSE,aes(x=IV, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the complementary data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_CFMRpv6.jpeg", device="jpeg",width = 4.35,height = 3.82)
ggplot(df_IV, aes(x = GWAS, y=(BMI-IV)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the complementary data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_box_CFMRpv6.jpeg", device="jpeg",width = 4.35,height = 3.82)
fit_test6 <- lm(D_init$BMI~IV_all)
summary(fit_test6)


df_train <- do.call(rbind, tl)
df_train <- data.frame(df_train)
colnames(df_train) <- c("Prediction","BMI", "GWAS")
df_train$GWAS <- as.factor(df_train$GWAS)


ggplot()+
  geom_point(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),alpha=0.1 )+
  geom_smooth(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_train,method="lm",se=FALSE,aes(x=Prediction, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the trainning data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")
ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_CFMRpv6.jpeg", device="jpeg",width = 4.35,height = 3.82)

ggplot(df_train, aes(x = GWAS, y=(BMI- Prediction)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the trainning data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")
ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_box_CFMRpv6.jpeg", device="jpeg",width = 4.35,height = 3.82)
fit_train6 <- lm(BMI~Prediction, df_train)
summary(fit_train6)

##### PVthres 10-7--------

set.seed(1)
options(mc.cores = 20)
#####Building IV1 ----
setwd(dir = paste("/mnt/work/william.denault/CFMR_Analysis/Data/", sep="") )
tt <- fread(paste(getwd(),"/Split_",1,"_SNP_for_extraction.txt",sep=""))
X <- subdf
 
colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
subtt <- tt[which(tt$P <10^(-7)),]
X <- X[ ,which(colnames(X) %in% subtt$ID)]




D_init$BMI[which(D_init$BMI ==-9)] <- NA
 

IV_all <- rep(NA,dim(D_init)[1])
#PREG ID in split1
Split_1_PREG_ID <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_1_Maternal_BMI.txt")
Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
#removing the one set as NA in the splitting procedure
sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
test_ID  <- (1:dim(D_init)[1])[-train_ID]

####Multiple NA in the observed data and glmnet does not handle them
y <- as.vector(D_init[train_ID,2])




x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
for( j in 1:dim( x)[2])
{
  x[,j] <- as.numeric(X[,j])
   
}
xtrain <- x[train_ID,]
xtest <- x[test_ID,]


Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1) 
save(Lasso_IV, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV1.RData")








#####Building IVs-----


#####Building the IVs ----
print("about to start buildin IV")
IVs_build <- function ( i)
  #for (i in 2:10)
{
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
   
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-7)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  sum( is.na(y))
   
   
  
   
  
  x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
  for( j in 1:dim( x)[2])
  {
    x[,j] <- as.numeric(X[,j])
     
  }
  xtrain <- x[train_ID,]
  xtest <- x[test_ID,]
  
  
  
  Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1)
  
  
  
   
  
   
  save(Lasso_IV , file = paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  print(paste("IV",i, "done"))
}

lapply(2:10, IVs_build)


tl <- list()
idx_list <-list()
mod_pred <- rep(NA, length(IV_all))
for( i in 1:10)
{
  
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
    
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-7)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  print(dim(X))
  
  length(which(!complete.cases(X)))
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  load(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  
  
  
  
  
  IV_all[test_ID] <- predict(Lasso_IV,X[test_ID,], s = "lambda.1se")
  tl[[i]] <- cbind(predict(Lasso_IV,X[train_ID,], s = "lambda.min"),y, rep(i,length(y)))
  
   
  idx_list[[i]] <- test_ID
   
  mod_pred[test_ID] <- rep( i, length(test_ID))
   
   
  Coef <- coef(Lasso_IV, s = "lambda.1se")
  print(length(Coef@i))
   
  
}


hist( IV_all)
IV_all[which(is.na(D_init$BMI))]  <- NA
plot( IV_all, D_init$BMI)

df_IV <- data.frame( IID= D_init$IID, IV=IV_all, BMI=D_init$BMI, GWAS=as.factor(mod_pred))
save(df_IV,file ="/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv7.RData")
ggplot()+
  geom_point(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),alpha=0.2 )+
  geom_smooth(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_IV,method="lm",se=FALSE,aes(x=IV, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the complementary data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_CFMRpv7.jpeg", device="jpeg",width = 4.35,height = 3.82)
ggplot(df_IV, aes(x = GWAS, y=(BMI-IV)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the complementary data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_box_CFMRpv7.jpeg", device="jpeg",width = 4.35,height = 3.82)
fit_test7 <- lm(D_init$BMI~IV_all)
summary(fit_test7)


df_train <- do.call(rbind, tl)
df_train <- data.frame(df_train)
colnames(df_train) <- c("Prediction","BMI", "GWAS")
df_train$GWAS <- as.factor(df_train$GWAS)


ggplot()+
  geom_point(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),alpha=0.1 )+
  geom_smooth(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_train,method="lm",se=FALSE,aes(x=Prediction, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the trainning data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")
ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_CFMRpv7.jpeg", device="jpeg",width = 4.35,height = 3.82)

ggplot(df_train, aes(x = GWAS, y=(BMI- Prediction)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the trainning data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")
ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_box_CFMRpv7.jpeg", device="jpeg",width = 4.35,height = 3.82)
fit_train7 <- lm(BMI~Prediction, df_train)
summary(fit_train7)

##### PVthres 10-8--------

set.seed(1)
options(mc.cores = 20)
#####Building IV1 ----
setwd(dir = paste("/mnt/work/william.denault/CFMR_Analysis/Data/", sep="") )
tt <- fread(paste(getwd(),"/Split_",1,"_SNP_for_extraction.txt",sep=""))
X <- subdf
 
colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
subtt <- tt[which(tt$P <10^(-8)),]
X <- X[ ,which(colnames(X) %in% subtt$ID)]




D_init$BMI[which(D_init$BMI ==-9)] <- NA
 

IV_all <- rep(NA,dim(D_init)[1])
#PREG ID in split1
Split_1_PREG_ID <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_1_Maternal_BMI.txt")
Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
#removing the one set as NA in the splitting procedure
sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
test_ID  <- (1:dim(D_init)[1])[-train_ID]

####Multiple NA in the observed data and glmnet does not handle them
y <- as.vector(D_init[train_ID,2])




x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
for( j in 1:dim( x)[2])
{
  x[,j] <- as.numeric(X[,j])
   
}
xtrain <- x[train_ID,]
xtest <- x[test_ID,]


Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1) 
save(Lasso_IV, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV1.RData")








#####Building IVs-----


#####Building the IVs ----
print("about to start buildin IV")
IVs_build <- function ( i)
  #for (i in 2:10)
{
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
   
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-8)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  sum( is.na(y))
   
   
  
   
  
  x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
  for( j in 1:dim( x)[2])
  {
    x[,j] <- as.numeric(X[,j])
     
  }
  xtrain <- x[train_ID,]
  xtest <- x[test_ID,]
  
  
  
  Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1)
  
  
  
   
  
   
  save(Lasso_IV , file = paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  print(paste("IV",i, "done"))
}

lapply(2:10, IVs_build)


tl <- list()
idx_list <-list()
mod_pred <- rep(NA, length(IV_all))
for( i in 1:10)
{
  
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
    
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-8)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  print(dim(X))
  
  length(which(!complete.cases(X)))
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  load(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  
  
  
  
  
  IV_all[test_ID] <- predict(Lasso_IV,X[test_ID,], s = "lambda.1se")
  tl[[i]] <- cbind(predict(Lasso_IV,X[train_ID,], s = "lambda.min"),y, rep(i,length(y)))
  
   
  idx_list[[i]] <- test_ID
   
  mod_pred[test_ID] <- rep( i, length(test_ID))
   
   
  Coef <- coef(Lasso_IV, s = "lambda.1se")
  print(length(Coef@i))
   
  
}


hist( IV_all)
IV_all[which(is.na(D_init$BMI))]  <- NA
plot( IV_all, D_init$BMI)

df_IV <- data.frame( IID= D_init$IID, IV=IV_all, BMI=D_init$BMI, GWAS=as.factor(mod_pred))
save(df_IV,file ="/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv8.RData")
ggplot()+
  geom_point(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),alpha=0.2 )+
  geom_smooth(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_IV,method="lm",se=FALSE,aes(x=IV, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the complementary data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_CFMRpv8.jpeg", device="jpeg",width = 4.35,height = 3.82)
ggplot(df_IV, aes(x = GWAS, y=(BMI-IV)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the complementary data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")

fit_test8 <- lm(D_init$BMI~IV_all)
summary(fit_test8)
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_box_CFMRpv8.jpeg", device="jpeg",width = 4.35,height = 3.82)

df_train <- do.call(rbind, tl)
df_train <- data.frame(df_train)
colnames(df_train) <- c("Prediction","BMI", "GWAS")
df_train$GWAS <- as.factor(df_train$GWAS)


ggplot()+
  geom_point(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),alpha=0.1 )+
  geom_smooth(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_train,method="lm",se=FALSE,aes(x=Prediction, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the trainning data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")
ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_CFMRpv8.jpeg", device="jpeg",width = 4.35,height = 3.82)

ggplot(df_train, aes(x = GWAS, y=(BMI- Prediction)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the trainning data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")
ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_box_CFMRpv8.jpeg", device="jpeg",width = 4.35,height = 3.82)
fit_train8 <- lm(BMI~Prediction, df_train)
summary(fit_train8)


##### PVthres 10-9--------

set.seed(1)
options(mc.cores = 20)
#####Building IV1 ----
setwd(dir = paste("/mnt/work/william.denault/CFMR_Analysis/Data/", sep="") )
tt <- fread(paste(getwd(),"/Split_",1,"_SNP_for_extraction.txt",sep=""))
X <- subdf
 
colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
subtt <- tt[which(tt$P <10^(-9)),]
X <- X[ ,which(colnames(X) %in% subtt$ID)]




D_init$BMI[which(D_init$BMI ==-9)] <- NA
 

IV_all <- rep(NA,dim(D_init)[1])
#PREG ID in split1
Split_1_PREG_ID <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_1_Maternal_BMI.txt")
Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
#removing the one set as NA in the splitting procedure
sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
test_ID  <- (1:dim(D_init)[1])[-train_ID]

####Multiple NA in the observed data and glmnet does not handle them
y <- as.vector(D_init[train_ID,2])




x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
for( j in 1:dim( x)[2])
{
  x[,j] <- as.numeric(X[,j])
   
}
xtrain <- x[train_ID,]
xtest <- x[test_ID,]


Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1) 
save(Lasso_IV, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV1.RData")








#####Building IVs-----


#####Building the IVs ----
print("about to start buildin IV")
IVs_build <- function ( i)
  #for (i in 2:10)
{
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
   
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-9)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  sum( is.na(y))
   
   
  
   
  
  x =matrix(NA,nrow=dim(X)[1], ncol=dim(X)[2])
  for( j in 1:dim( x)[2])
  {
    x[,j] <- as.numeric(X[,j])
     
  }
  xtrain <- x[train_ID,]
  xtest <- x[test_ID,]
  
  
  
  Lasso_IV <- cv.glmnet(x=xtrain, y=y, alpha=1)
  
  
  
   
  
   
  save(Lasso_IV , file = paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  print(paste("IV",i, "done"))
}

lapply(2:10, IVs_build)


tl <- list()
idx_list <-list()
mod_pred <- rep(NA, length(IV_all))
for( i in 1:10)
{
  
  tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
   
  X <- subdf
    
  colnames(X) <- colnames(Reg_MAT[,-c(grep("SEN", colnames(Reg_MAT)), grep("M_ID",colnames(Reg_MAT)))])
  
  subtt <- tt[which(tt$P <10^(-9)),]
  X <- X[ ,which(colnames(X) %in% subtt$ID)]
  print(dim(X))
  
  length(which(!complete.cases(X)))
  Split_1_PREG_ID <- fread(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Split_",i,"_Maternal_BMI.txt",sep=""))
  Split_1_PREG_ID <- Split_1_PREG_ID[order(Split_1_PREG_ID$IID),]
  train_ID <- which(D_init$IID %in% Split_1_PREG_ID$IID & !is.na(D_init[,2]) )
  #removing the one set as NA in the splitting procedure
  sub_indx <- Split_1_PREG_ID$IID[which(Split_1_PREG_ID$BMI==-9)]
  train_ID <- train_ID [ -which(D_init$IID[train_ID] %in% sub_indx)]
  test_ID  <- (1:dim(D_init)[1])[-train_ID]
  
  y <- as.vector(D_init[train_ID,2])
  load(paste("/mnt/work/william.denault/CFMR_Analysis/Data/Lasso_IV",i,".RData", sep=""))
  
  
  
  
  
  IV_all[test_ID] <- predict(Lasso_IV,X[test_ID,], s = "lambda.1se")
  tl[[i]] <- cbind(predict(Lasso_IV,X[train_ID,], s = "lambda.min"),y, rep(i,length(y)))
  
   
  idx_list[[i]] <- test_ID
   
  mod_pred[test_ID] <- rep( i, length(test_ID))
   
   
  Coef <- coef(Lasso_IV, s = "lambda.1se")
  print(length(Coef@i))
   
  
}


hist( IV_all)
IV_all[which(is.na(D_init$BMI))]  <- NA
plot( IV_all, D_init$BMI)

df_IV <- data.frame( IID= D_init$IID, IV=IV_all, BMI=D_init$BMI, GWAS=as.factor(mod_pred))
save(df_IV,file ="/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv9.RData")
ggplot()+
  geom_point(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),alpha=0.2 )+
  geom_smooth(data=df_IV, aes(x = IV, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_IV,method="lm",se=FALSE,aes(x=IV, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the complementary data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_CFMRpv9.jpeg", device="jpeg",width = 4.35,height = 3.82)
ggplot(df_IV, aes(x = GWAS, y=(BMI-IV)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the complementary data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")

fit_test9 <- lm(D_init$BMI~IV_all)
summary(fit_test9)
ggsave("/mnt/cargo/william.denault/graph_CFMR/test_perf_box_CFMRpv9.jpeg", device="jpeg",width = 4.35,height = 3.82)

df_train <- do.call(rbind, tl)
df_train <- data.frame(df_train)
colnames(df_train) <- c("Prediction","BMI", "GWAS")
df_train$GWAS <- as.factor(df_train$GWAS)


ggplot()+
  geom_point(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),alpha=0.1 )+
  geom_smooth(data=df_train, aes(x = Prediction, y=BMI,colour=GWAS),method="lm",se=FALSE)+
  geom_smooth(data=df_train,method="lm",se=FALSE,aes(x=Prediction, y=BMI),linetype = "dashed",colour="black", size=2)+
  
  ggtitle("Observed pre-pregnancy BMI against predicted pre-pregnancy BMI \n predictions based on the trainning data of each GWAS")+
  xlab("predicted pre-pregnancy BMI")+
  ylab("pre-pregnancy BMI")
ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_CFMRpv9.jpeg", device="jpeg",width = 4.35,height = 3.82)

ggplot(df_train, aes(x = GWAS, y=(BMI- Prediction)))+
  geom_boxplot()+
  ggtitle("Boxplot of the prediction error of each predictor\n predictions based on the trainning data of each GWAS ")+
  xlab("GWAS")+
  ylab("Error in prediction")
ggsave("/mnt/cargo/william.denault/graph_CFMR/training_perf_box_CFMRpv9.jpeg", device="jpeg",width = 4.35,height = 3.82)
fit_train9 <- lm(BMI~Prediction, df_train)
summary(fit_train9)


save(
     fit_train4,
     fit_train5,
     fit_train6,
     fit_train7,
     fit_train8,fit_train9,
     
     file="test_IV_train.Rdata")
save(
  
    fit_test4,
    fit_test5,
    fit_test6,
    fit_test7,
    fit_test8,
    fit_test9,file="test_IV_test.Rdata")
