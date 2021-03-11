rm(list=ls())
load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf.RData")
setwd(dir = "/home/william.denault/Causal_ART")
source("utils.R")
library(glmnet)
library(hdm)
library(rms)
library(data.table)

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

IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374


load("~/archive/START/pheno/PDB_2374_MoBa_Q1_Q4_MBRN_full.RData")
MoBa_Q1_Q4_MBRN_full <- MoBa_Q1_Q4_MBRN_full[order(MoBa_Q1_Q4_MBRN_full$PREG_ID_2374),]
phenofile <- MoBa_Q1_Q4_MBRN_full[which(MoBa_Q1_Q4_MBRN_full$PREG_ID_2374 %in% preg_id_analysis),]
phenofile <- phenofile[-which(duplicated(phenofile$PREG_ID_2374 )),]


phenofile$MSmoking <- rep( 0,dim(phenofile)[1])
phenofile$MSmoking[which(phenofile$AA1348 =="Yes")] <-1
### Estimating effect for different p thrsh -----
res<- list()
res_norm <- list()
res_counf <- list()
l <- 0

##### -3 -----
load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv3.RData")
l <- l+1

IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374
IV_df <- IV_df[order(IV_df$PREG_ID_2374),]
Y <- phenofile$VEKT
Y[which(Y>7000)] <-NA
Y[which(Y<400)] <-NA
GA <- phenofile$SVLEN_DG
SEX <- as.factor(phenofile$KJONN)
BMI <- IV_df$BMI
IV <- IV_df$IV
MSmoking <- phenofile$MSmoking
dfreg <- data.frame( BMI = BMI,
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)


fm <- ivreg(Y~BMI+GA+SEX+MSmoking|GA+SEX+MSmoking+IV, data=dfreg)
summary(fm)
confint(fm)[2,]
f_ob <-lm(Y~BMI,data=dfreg)

res[[l]]   <- c(-3,
                summary(fm)$coefficients[2,],
                confint(fm)[2,]
)

###Normalized version---
sd(BMI, na.rm=TRUE)
sd(BMI, na.rm=TRUE)
dfreg <- data.frame( BMI = BMI/sd(BMI, na.rm=TRUE),
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)
fm <- ivreg(Y~BMI+GA+SEX+MSmoking|GA+SEX+MSmoking+IV, data=dfreg)
summary(fm)
confint(fm)[2,]


res_norm[[l]]   <- c(-3,
                     summary(fm)$coefficients[2,],
                     confint(fm)[2,]
)

sage <- summary(lm(phenofile$MORS_ALDER~IV_df$IV))
sga <-  summary(lm(phenofile$SVLEN ~IV_df$IV))
Medu1 <- rep(0,length(IV))#check with 4 or more years of study
Medu1[which(phenofile$AA1124 %in% c("Higher education (university/college), over 4 years","Higher education (university/college), up to and including 4 years"))] <-1
smedu <- summary(glm(Medu1~IV, family = "binomial"))

MSmoking <- rep( 0,length(IV))
MSmoking[which(phenofile$AA1348 =="Yes")] <-1#smoking prior to pregnancy

ssmo <- summary(glm(MSmoking~IV, family = "binomial"))

res_counf[[l]] <- 
  
  rbind( c( -3,"Mother_Age",sage$coefficients[2,]),
         c( -3,"Mother_edu",smedu$coefficients[2,]),
         c( -3,"Gestational Age",sga$coefficients[2,]),
         c( -3,"Mother_Smoking",ssmo$coef[2,])
  )



##### -4 -----
load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv4.RData")
l <- l+1

IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374
IV_df <- IV_df[order(IV_df$PREG_ID_2374),]
Y <- phenofile$VEKT
Y[which(Y>7000)] <-NA
Y[which(Y<400)] <-NA
GA <- phenofile$SVLEN_DG
SEX <- as.factor(phenofile$KJONN)
BMI <- IV_df$BMI
IV <- IV_df$IV
MSmoking <- phenofile$MSmoking
dfreg <- data.frame( BMI = BMI,
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)


fm <- ivreg(Y~BMI+GA+SEX+MSmoking|GA+SEX+MSmoking+IV, data=dfreg)
summary(fm)
confint(fm)[2,]
f_ob <-lm(Y~BMI,data=dfreg)

res[[l]]   <- c(-4,
                summary(fm)$coefficients[2,],
                confint(fm)[2,]
)

###Normalized version---
sd(BMI, na.rm=TRUE)
sd(BMI, na.rm=TRUE)
dfreg <- data.frame( BMI = BMI/sd(BMI, na.rm=TRUE),
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)
fm <- ivreg(Y~BMI+GA+SEX+MSmoking|GA+SEX+MSmoking+IV, data=dfreg)
summary(fm)
confint(fm)[2,]


res_norm[[l]]   <- c(-4,
                     summary(fm)$coefficients[2,],
                     confint(fm)[2,]
)

sage <- summary(lm(phenofile$MORS_ALDER~IV_df$IV))
sga <-  summary(lm(phenofile$SVLEN ~IV_df$IV))
Medu1 <- rep(0,length(IV))#check with 4 or more years of study
Medu1[which(phenofile$AA1124 %in% c("Higher education (university/college), over 4 years","Higher education (university/college), up to and including 4 years"))] <-1
smedu <- summary(glm(Medu1~IV, family = "binomial"))

MSmoking <- rep( 0,length(IV))
MSmoking[which(phenofile$AA1348 =="Yes")] <-1#smoking prior to pregnancy

ssmo <- summary(glm(MSmoking~IV, family = "binomial"))

res_counf[[l]] <- 
  
  rbind( c( -4,"Mother_Age",sage$coefficients[2,]),
         c( -4,"Mother_edu",smedu$coefficients[2,]),
         c( -4,"Gestational Age",sga$coefficients[2,]),
         c( -4,"Mother_Smoking",ssmo$coef[2,])
  )



##### -5------
load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv5.RData")
l <- l+1
IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374

IV_df <- IV_df[order(IV_df$PREG_ID_2374),]
Y <- phenofile$VEKT
Y[which(Y>7000)] <-NA
Y[which(Y<400)] <-NA
GA <- phenofile$SVLEN_DG
SEX <- as.factor(phenofile$KJONN)
BMI <- IV_df$BMI
IV <- IV_df$IV
dfreg <- data.frame( BMI = BMI,
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)



fm <- ivreg(Y~BMI+GA+SEX+MSmoking|GA+SEX+MSmoking+IV, data=dfreg)
summary(fm)
confint(fm)[2,]
f_ob <-lm(Y~BMI,data=dfreg)

res[[l]]   <- c(-5,
                summary(fm)$coefficients[2,],
                confint(fm)[2,]
)

###Normalized version---
sd(BMI, na.rm=TRUE)
dfreg <- data.frame( BMI = BMI/sd(BMI, na.rm=TRUE),
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)
fm <- ivreg(Y~BMI+GA+SEX+MSmoking|GA+SEX+MSmoking+IV, data=dfreg)
summary(fm)
confint(fm)[2,]


res_norm[[l]]   <- c(-5,
                     summary(fm)$coefficients[2,],
                     confint(fm)[2,]
)


sage <- summary(lm(phenofile$MORS_ALDER~IV_df$IV))
sga <-  summary(lm(phenofile$SVLEN ~IV_df$IV))
Medu1 <- rep(0,length(IV))#check with 4 or more years of study
Medu1[which(phenofile$AA1124 %in% c("Higher education (university/college), over 4 years","Higher education (university/college), up to and including 4 years"))] <-1
smedu <- summary(glm(Medu1~IV, family = "binomial"))

MSmoking <- rep( 0,length(IV))
MSmoking[which(phenofile$AA1348 =="Yes")] <-1#smoking prior to pregnancy

ssmo <- summary(glm(MSmoking~IV, family = "binomial"))

res_counf[[l]] <- 
  
  rbind( c( -5,"Mother_Age",sage$coefficients[2,]),
         c( -5,"Mother_edu",smedu$coefficients[2,]),
         c( -5,"Gestational Age",sga$coefficients[2,]),
         c( -5,"Mother_Smoking",ssmo$coef[2,])
  )



###### -6 -------
load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv6.RData")
l <- l+1
IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374

IV_df <- IV_df[order(IV_df$PREG_ID_2374),]
Y <- phenofile$VEKT
Y[which(Y>7000)] <-NA
Y[which(Y<400)] <-NA
GA <- phenofile$SVLEN_DG
SEX <- as.factor(phenofile$KJONN)
BMI <- IV_df$BMI
IV <- IV_df$IV
dfreg <- data.frame( BMI = BMI,
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)



fm <- ivreg(Y~BMI+GA+SEX+MSmoking|GA+SEX+MSmoking+IV, data=dfreg)
summary(fm)
confint(fm)[2,]
f_ob <-lm(Y~BMI,data=dfreg)

res[[l]]   <- c(-6,
                summary(fm)$coefficients[2,],
                confint(fm)[2,]
)

###Normalized version---
sd(BMI, na.rm=TRUE)
sd(BMI, na.rm=TRUE)
dfreg <- data.frame( BMI = BMI/sd(BMI, na.rm=TRUE),
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)
fm <- ivreg(Y~BMI+GA+SEX+MSmoking|GA+SEX+MSmoking+IV, data=dfreg)
summary(fm)
confint(fm)[2,]


res_norm[[l]]   <- c(-6,
                     summary(fm)$coefficients[2,],
                     confint(fm)[2,]
)
sage <- summary(lm(phenofile$MORS_ALDER~IV_df$IV))
sga <-  summary(lm(phenofile$SVLEN ~IV_df$IV))
Medu1 <- rep(0,length(IV))#check with 4 or more years of study
Medu1[which(phenofile$AA1124 %in% c("Higher education (university/college), over 4 years","Higher education (university/college), up to and including 4 years"))] <-1
smedu <- summary(glm(Medu1~IV, family = "binomial"))

MSmoking <- rep( 0,length(IV))
MSmoking[which(phenofile$AA1348 =="Yes")] <-1#smoking prior to pregnancy

ssmo <- summary(glm(MSmoking~IV, family = "binomial"))

res_counf[[l]] <- 
  
  rbind( c( -6,"Mother_Age",sage$coefficients[2,]),
         c( -6,"Mother_edu",smedu$coefficients[2,]),
         c( -6,"Gestational Age",sga$coefficients[2,]),
         c( -6,"Mother_Smoking",ssmo$coef[2,])
  )


#### -7 -------
load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv7.RData")
l <- l+1
IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374

IV_df <- IV_df[order(IV_df$PREG_ID_2374),]
Y <- phenofile$VEKT
Y[which(Y>7000)] <-NA
Y[which(Y<400)] <-NA
GA <- phenofile$SVLEN_DG
SEX <- as.factor(phenofile$KJONN)
BMI <- IV_df$BMI
IV <- IV_df$IV
dfreg <- data.frame( BMI = BMI,
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)



fm <- ivreg(Y~BMI+GA+SEX+MSmoking|GA+SEX+MSmoking+IV, data=dfreg)
summary(fm)
confint(fm)[2,]
f_ob <-lm(Y~BMI,data=dfreg)

res[[l]]   <- c(-7,
                summary(fm)$coefficients[2,],
                confint(fm)[2,]
)

###Normalized version---
sd(BMI, na.rm=TRUE)
dfreg <- data.frame( BMI = BMI/sd(BMI, na.rm=TRUE),
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)
fm <- ivreg(Y~BMI+GA+SEX+MSmoking|GA+SEX+MSmoking+IV, data=dfreg)
summary(fm)
confint(fm)[2,]


res_norm[[l]]   <- c(-7,
                     summary(fm)$coefficients[2,],
                     confint(fm)[2,]
)

sage <- summary(lm(phenofile$MORS_ALDER~IV_df$IV))
sga <-  summary(lm(phenofile$SVLEN ~IV_df$IV))
Medu1 <- rep(0,length(IV))#check with 4 or more years of study
Medu1[which(phenofile$AA1124 %in% c("Higher education (university/college), over 4 years","Higher education (university/college), up to and including 4 years"))] <-1
smedu <- summary(glm(Medu1~IV, family = "binomial"))

MSmoking <- rep( 0,length(IV))
MSmoking[which(phenofile$AA1348 =="Yes")] <-1#smoking prior to pregnancy

ssmo <- summary(glm(MSmoking~IV, family = "binomial"))

res_counf[[l]] <- 
  
  rbind( c( -7,"Mother_Age",sage$coefficients[2,]),
         c( -7,"Mother_edu",smedu$coefficients[2,]),
         c( -7,"Gestational Age",sga$coefficients[2,]),
         c( -7,"Mother_Smoking",ssmo$coef[2,])
  )


##### -8 -----
load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv8.RData")
l <- l+1
IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374

IV_df <- IV_df[order(IV_df$PREG_ID_2374),]
Y <- phenofile$VEKT
Y[which(Y>7000)] <-NA
Y[which(Y<400)] <-NA
GA <- phenofile$SVLEN_DG
SEX <- as.factor(phenofile$KJONN)
BMI <- IV_df$BMI
IV <- IV_df$IV
dfreg <- data.frame( BMI = BMI,
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)



fm <- ivreg(Y~BMI+GA+SEX+MSmoking|GA+SEX+MSmoking+IV, data=dfreg)
summary(fm)
confint(fm)[2,]
f_ob <-lm(Y~BMI,data=dfreg)

res[[l]]   <- c(-8,
                summary(fm)$coefficients[2,],
                confint(fm)[2,]
)

###Normalized version---
sd(BMI, na.rm=TRUE)
dfreg <- data.frame( BMI = BMI/sd(BMI, na.rm=TRUE),
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)
fm <- ivreg(Y~BMI+GA+SEX+MSmoking|GA+SEX+MSmoking+IV, data=dfreg)
summary(fm)
confint(fm)[2,]


res_norm[[l]]   <- c(-8,
                     summary(fm)$coefficients[2,],
                     confint(fm)[2,]
)




f_ob <-lm(Y~BMI,data=dfreg)
summary(fm)
confint(fm)[2,]




#Highest educational qualification attained
#Gestational/existing diabetes
#Smoking (current smoker vs non-smoker)
#Age


sage <- summary(lm(phenofile$MORS_ALDER~IV_df$IV))
sga <-  summary(lm(phenofile$SVLEN ~IV_df$IV))
Medu1 <- rep(0,length(IV))#check with 4 or more years of study
Medu1[which(phenofile$AA1124 %in% c("Higher education (university/college), over 4 years","Higher education (university/college), up to and including 4 years"))] <-1
smedu <- summary(glm(Medu1~IV, family = "binomial"))

MSmoking <- rep( 0,length(IV))
MSmoking[which(phenofile$AA1348 =="Yes")] <-1#smoking prior to pregnancy

ssmo <- summary(glm(MSmoking~IV, family = "binomial"))

res_counf[[l]] <- 
  
  rbind( c( -8,"Mother_Age",sage$coefficients[2,]),
         c( -8,"Mother_edu",smedu$coefficients[2,]),
         c( -8,"Gestational Age",sga$coefficients[2,]),
         c( -8,"Mother_Smoking",ssmo$coef[2,])
  )

####Plots ----
res <- do.call(rbind,res )
res_norm <- do.call(rbind,res_norm )
res_counf <- do.call(rbind, res_counf)

res <- data.frame(res)
colnames(res)[c(1,6,7)] <- c("pv","low", "up")
res$pv <- as.factor(res$pv)

res_norm <- data.frame(res_norm)
colnames(res_norm)[c(1,6,7)] <- c("pv","low", "up")
res_norm$pv <- as.factor(res_norm$pv)

ggplot(res,aes(x=pv, y=Estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=low, ymax=up))

res_normT <- res_norm[, c(1,2,6,7)]

rbind( res_normT, c( "Tyrrel et al." ))
P1 <- ggplot(res_norm,aes(x=pv, y=Estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=low, ymax=up))+
  ggtitle("CFMR estimates of pre-pregancy BMI effect\n on birth weight with 95%confidence intervals ")+
  xlab("-log10 p-value")+
  ylab("CFMR estimates")
P1 
ggsave(P1, file="/mnt/cargo/william.denault/graph_CFMR/CFMR_est.png",width = 4.35,height = 3.82)
weight <- (res[,7]-res[,6])
weighted.mean(res$Estimate, w=weight)
weight <- (res_norm[,7]-res_norm[,6])
weighted.mean(res_norm$Estimate, w=weight)


##### Var explained -3 ----
l <- 1
var_exp <- list()
load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv3.RData")
IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374

IV_df <- IV_df[order(IV_df$PREG_ID_2374),]
Y <- phenofile$VEKT
Y[which(Y>7000)] <-NA
Y[which(Y<400)] <-NA
GA <- phenofile$SVLEN_DG
SEX <- as.factor(phenofile$KJONN)
BMI <- IV_df$BMI
IV <- IV_df$IV
dfreg <- data.frame( BMI = BMI,
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)

lol <-summary(lm(BMI~IV,dfreg))
var_exp[[l]] <- c(-3,lol$adj.r.squared)
##### Var explained -4 ----
l <- 1+l

load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv4.RData")
IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374

IV_df <- IV_df[order(IV_df$PREG_ID_2374),]
Y <- phenofile$VEKT
Y[which(Y>7000)] <-NA
Y[which(Y<400)] <-NA
GA <- phenofile$SVLEN_DG
SEX <- as.factor(phenofile$KJONN)
BMI <- IV_df$BMI
IV <- IV_df$IV
dfreg <- data.frame( BMI = BMI,
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)

lol <-summary(lm(BMI~IV,dfreg))
var_exp[[l]] <- c(-4,lol$adj.r.squared)
##### Var explained -5 ----
l <- 1+l

load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv5.RData")
IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374

IV_df <- IV_df[order(IV_df$PREG_ID_2374),]
Y <- phenofile$VEKT
Y[which(Y>7000)] <-NA
Y[which(Y<400)] <-NA
GA <- phenofile$SVLEN_DG
SEX <- as.factor(phenofile$KJONN)
BMI <- IV_df$BMI
IV <- IV_df$IV
dfreg <- data.frame( BMI = BMI,
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)

lol <-summary(lm(BMI~IV,dfreg))
var_exp[[l]] <- c(-5,lol$adj.r.squared)
##### Var explained -6 ----
l <- 1+l

load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv6.RData")
IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374

IV_df <- IV_df[order(IV_df$PREG_ID_2374),]
Y <- phenofile$VEKT
Y[which(Y>7000)] <-NA
Y[which(Y<400)] <-NA
GA <- phenofile$SVLEN_DG
SEX <- as.factor(phenofile$KJONN)
BMI <- IV_df$BMI
IV <- IV_df$IV
dfreg <- data.frame( BMI = BMI,
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)

lol <-summary(lm(BMI~IV,dfreg))
var_exp[[l]] <- c(-6,lol$adj.r.squared)
##### Var explained -7 ----
l <- 1+l

load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv7.RData")
IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374

IV_df <- IV_df[order(IV_df$PREG_ID_2374),]
Y <- phenofile$VEKT
Y[which(Y>7000)] <-NA
Y[which(Y<400)] <-NA
GA <- phenofile$SVLEN_DG
SEX <- as.factor(phenofile$KJONN)
BMI <- IV_df$BMI
IV <- IV_df$IV
dfreg <- data.frame( BMI = BMI,
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)

lol <-summary(lm(BMI~IV,dfreg))
var_exp[[l]] <- c(-7,lol$adj.r.squared)
##### Var explained -8 ----
l <- 1+l

load("/mnt/work/william.denault/CFMR_Analysis/Data/IVdf_pv8.RData")
IV_df <- merge(df_IV, infom, by.x="IID", by.y="SENTRIXID")
IV_df <- IV_df[-which(duplicated(IV_df$PREG_ID_2374)),]
preg_id_analysis <- IV_df$PREG_ID_2374

IV_df <- IV_df[order(IV_df$PREG_ID_2374),]
Y <- phenofile$VEKT
Y[which(Y>7000)] <-NA
Y[which(Y<400)] <-NA
GA <- phenofile$SVLEN_DG
SEX <- as.factor(phenofile$KJONN)
BMI <- IV_df$BMI
IV <- IV_df$IV
dfreg <- data.frame( BMI = BMI,
                     IV= IV, 
                     GA=GA,
                     SEX=SEX,
                     Y=Y
)

lol <-summary(lm(BMI~IV,dfreg))
var_exp[[l]] <- c(-8,lol$adj.r.squared)



var_exp <- do.call(rbind, var_exp)
save(res,res_norm, res_counf, var_exp, file="/mnt/cargo/william.denault/CFMR_estimate.RData")
res
res_norm
res_counf
