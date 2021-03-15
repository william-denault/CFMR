library(data.table)
##Cleanning Results ----
rm(list=ls())
library(data.table)

maf <- 0.05
setwd("/home/william.denault/archive/MERGE/markerinfo")
info_marker<- fread("1-markerinfo")
info_marker <- info_marker[ which(info_marker$`[8]RefPanelAF` >maf ),c(3,7,8)]
for( i in 2:22)
{
  temp <- fread(paste (i,"-markerinfo",sep=""))
  temp <- temp[ which(temp$`[8]RefPanelAF` >maf ),c(3,7,8)]
  info_marker <- rbind(info_marker,temp)
  print(i)
}

colnames(info_marker) <- c("ID", "INFO","RefPanelAF")
library(data.table)

maf <- 0.05
setwd("/home/william.denault/archive/MERGE/markerinfo")
info_marker<- fread("1-markerinfo")
info_marker <- info_marker[ which(info_marker$`[8]RefPanelAF` >maf ),c(3,7,8)]
for( i in 2:22)
{
  temp <- fread(paste (i,"-markerinfo",sep=""))
  temp <- temp[ which(temp$`[8]RefPanelAF` >maf ),c(3,7,8)]
  info_marker <- rbind(info_marker,temp)
  print(i)
}

colnames(info_marker) <- c("ID", "INFO","RefPanelAF")
#Split 1 ----


setwd("/mnt/work/william.denault/CFMR_Analysis/Data/Split_1_GWAS_Maternal")
lf <- list.files()
lf <- lf[grep("glm", lf)] #keep only the results files

cleaning_f <- function(x)
{
  temp <- fread(x)
  print(length(is.na(temp$P)))
  temp <- temp[ - which( is.na(temp$P)),]
  return(temp)
}


GWAS_mother <- do.call(rbind,lapply(lf, cleaning_f))
GWAS_mother <- merge( GWAS_mother,info_marker, by="ID") #removing low frequency variant
GWAS_mother <- GWAS_mother[order(GWAS_mother$`#CHROM`,GWAS_mother$POS),]
write.table(GWAS_mother, file ="Split_1_GWAS_mother.txt", quote=FALSE, row.names = FALSE)

#Split 2 ----


setwd("/mnt/work/william.denault/CFMR_Analysis/Data/Split_2_GWAS_Maternal")
lf <- list.files()
lf <- lf[grep("glm", lf)] #keep only the results files

cleaning_f <- function(x)
{
  temp <- fread(x)
  print(length(is.na(temp$P)))
  temp <- temp[ - which( is.na(temp$P)),]
  return(temp)
}


GWAS_mother <- do.call(rbind,lapply(lf, cleaning_f))
GWAS_mother <- merge( GWAS_mother,info_marker, by="ID") #removing low frequency variant
GWAS_mother <- GWAS_mother[order(GWAS_mother$`#CHROM`,GWAS_mother$POS),]
write.table(GWAS_mother, file ="Split_2_GWAS_mother.txt", quote=FALSE, row.names = FALSE)



#Split 3 ----


setwd("/mnt/work/william.denault/CFMR_Analysis/Data/Split_3_GWAS_Maternal")
lf <- list.files()
lf <- lf[grep("glm", lf)] #keep only the results files

cleaning_f <- function(x)
{
  temp <- fread(x)
  print(length(is.na(temp$P)))
  temp <- temp[ - which( is.na(temp$P)),]
  return(temp)
}


GWAS_mother <- do.call(rbind,lapply(lf, cleaning_f))
GWAS_mother <- merge( GWAS_mother,info_marker, by="ID") #removing low frequency variant
GWAS_mother <- GWAS_mother[order(GWAS_mother$`#CHROM`,GWAS_mother$POS),]
write.table(GWAS_mother, file ="Split_3_GWAS_mother.txt", quote=FALSE, row.names = FALSE)


#Split 4 ----


setwd("/mnt/work/william.denault/CFMR_Analysis/Data/Split_4_GWAS_Maternal")
lf <- list.files()
lf <- lf[grep("glm", lf)] #keep only the results files

cleaning_f <- function(x)
{
  temp <- fread(x)
  print(length(is.na(temp$P)))
  temp <- temp[ - which( is.na(temp$P)),]
  return(temp)
}


GWAS_mother <- do.call(rbind,lapply(lf, cleaning_f))
GWAS_mother <- merge( GWAS_mother,info_marker, by="ID") #removing low frequency variant
GWAS_mother <- GWAS_mother[order(GWAS_mother$`#CHROM`,GWAS_mother$POS),]
write.table(GWAS_mother, file ="/mnt/work/william.denault/CFMR_Analysis/Data/Split_4_GWAS_Maternal/Split_4_GWAS_mother.txt", quote=FALSE, row.names = FALSE)



#Split 5 ----


setwd("/mnt/work/william.denault/CFMR_Analysis/Data/Split_5_GWAS_Maternal")
lf <- list.files("/mnt/work/william.denault/CFMR_Analysis/Data/Split_5_GWAS_Maternal")
lf <- lf[grep("glm", lf)] #keep only the results files

cleaning_f <- function(x)
{
  temp <- fread(x)
  print(length(is.na(temp$P)))
  temp <- temp[ - which( is.na(temp$P)),]
  return(temp)
}


GWAS_mother <- do.call(rbind,lapply(lf, cleaning_f))
GWAS_mother <- merge( GWAS_mother,info_marker, by="ID") #removing low frequency variant
GWAS_mother <- GWAS_mother[order(GWAS_mother$`#CHROM`,GWAS_mother$POS),]
write.table(GWAS_mother, file ="Split_5_GWAS_mother.txt", quote=FALSE, row.names = FALSE)



#Split 6 ----


setwd("/mnt/work/william.denault/CFMR_Analysis/Data/Split_6_GWAS_Maternal")
lf <- list.files()
lf <- lf[grep("glm", lf)] #keep only the results files

cleaning_f <- function(x)
{
  temp <- fread(x)
  print(length(is.na(temp$P)))
  temp <- temp[ - which( is.na(temp$P)),]
  return(temp)
}


GWAS_mother <- do.call(rbind,lapply(lf, cleaning_f))
GWAS_mother <- merge( GWAS_mother,info_marker, by="ID") #removing low frequency variant
GWAS_mother <- GWAS_mother[order(GWAS_mother$`#CHROM`,GWAS_mother$POS),]
write.table(GWAS_mother, file ="Split_6_GWAS_mother.txt", quote=FALSE, row.names = FALSE)



#Split 7 ----


setwd("/mnt/work/william.denault/CFMR_Analysis/Data/Split_7_GWAS_Maternal")
lf <- list.files()
lf <- lf[grep("glm", lf)] #keep only the results files

cleaning_f <- function(x)
{
  temp <- fread(x)
  print(length(is.na(temp$P)))
  temp <- temp[ - which( is.na(temp$P)),]
  return(temp)
}


GWAS_mother <- do.call(rbind,lapply(lf, cleaning_f))
GWAS_mother <- merge( GWAS_mother,info_marker, by="ID") #removing low frequency variant
GWAS_mother <- GWAS_mother[order(GWAS_mother$`#CHROM`,GWAS_mother$POS),]
write.table(GWAS_mother, file ="Split_7_GWAS_mother.txt", quote=FALSE, row.names = FALSE)



#Split 8 ----


setwd("/mnt/work/william.denault/CFMR_Analysis/Data/Split_8_GWAS_Maternal")
lf <- list.files()
lf <- lf[grep("glm", lf)] #keep only the results files

cleaning_f <- function(x)
{
  temp <- fread(x)
  print(length(is.na(temp$P)))
  temp <- temp[ - which( is.na(temp$P)),]
  return(temp)
}


GWAS_mother <- do.call(rbind,lapply(lf, cleaning_f))
GWAS_mother <- merge( GWAS_mother,info_marker, by="ID") #removing low frequency variant
GWAS_mother <- GWAS_mother[order(GWAS_mother$`#CHROM`,GWAS_mother$POS),]
write.table(GWAS_mother, file ="Split_8_GWAS_mother.txt", quote=FALSE, row.names = FALSE)



#Split 9 ----


setwd("/mnt/work/william.denault/CFMR_Analysis/Data/Split_9_GWAS_Maternal")
lf <- list.files()
lf <- lf[grep("glm", lf)] #keep only the results files

cleaning_f <- function(x)
{
  temp <- fread(x)
  print(length(is.na(temp$P)))
  temp <- temp[ - which( is.na(temp$P)),]
  return(temp)
}


GWAS_mother <- do.call(rbind,lapply(lf, cleaning_f))
GWAS_mother <- merge( GWAS_mother,info_marker, by="ID") #removing low frequency variant
GWAS_mother <- GWAS_mother[order(GWAS_mother$`#CHROM`,GWAS_mother$POS),]
write.table(GWAS_mother, file ="Split_9_GWAS_mother.txt", quote=FALSE, row.names = FALSE)

#Split 10 ----


setwd("/mnt/work/william.denault/CFMR_Analysis/Data/Split_10_GWAS_Maternal")
lf <- list.files()
lf <- lf[grep("glm", lf)] #keep only the results files

cleaning_f <- function(x)
{
  temp <- fread(x)
  temp <- temp[ - which( is.na(temp$P)),]
  return(temp)
}


GWAS_mother <- do.call(rbind,lapply(lf, cleaning_f))
GWAS_mother <- merge( GWAS_mother,info_marker, by="ID") #removing low frequency variant
GWAS_mother <- GWAS_mother[order(GWAS_mother$`#CHROM`,GWAS_mother$POS),]
write.table(GWAS_mother, file ="Split_10_GWAS_mother.txt", quote=FALSE, row.names = FALSE)









##Clumping Split 1 -----
#keeping only SNP with pv below 10^-3
#SNP s for propensity score (PS)

rm(list=ls())
library(data.table)
source("/home/william.denault/Causal_ART/utils.R")



Split_1_GWAS_mother <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_1_GWAS_Maternal/Split_1_GWAS_mother.txt")


jpeg(file="/mnt/cargo/william.denault/graph_CFMR/mh1.jpg", width = 850, height = 650)
manhattan.plot(chr= Split_1_GWAS_mother$`#CHROM`,
               pos= Split_1_GWAS_mother$POS,
               pvalue= Split_1_GWAS_mother$P)
abline(a=-log10(5*10^(-8)), b=0)
dev.off()
SNP_PS <- Split_1_GWAS_mother[ -which(Split_1_GWAS_mother$P  >10^(-2)),]

if( length(which(SNP_PS$ID == ".")>0))#to deal without ID
{
  for( i in which(SNP_PS$ID == "."))
  {
    SNP_PS$ID[i] <-  paste(SNP_PS$`#CHROM`[i],SNP_PS$POS[i], sep=":" )
  }
}

tt <-  SNP_PS[order(SNP_PS$P),]
print("Starting Clumping")

thresh <- 500000#max distancez between two SNPs in the PS
maxit  <- dim(tt)[1]


for( i in 1:maxit)
{
  temp <- which(tt$`#CHROM` == tt$`#CHROM`[i] )
  if( length(temp)>1)
  {
    tempd     <- abs( tt$POS - tt$POS[i] )
    tempd     <- which(tempd < thresh) #index SNPS to close
    to_remove <- temp[which(temp%in%tempd)]#the SNP on the same chromosome being to close to the lead SNP
    if ( length(to_remove)>1)#1 because the leading SNP is included
    {
      to_remove <- to_remove[ -which(to_remove ==i)]
      tt <- tt[-to_remove]
    }
  }
  if(dim(tt)[1]==i)
  {
    break()
  }
}

SNP_PS <- tt
head(SNP_PS)
dim(SNP_PS)#for 2 splits 1210 and 16
# 
write.table(SNP_PS, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Split_1_SNP_for_extraction.txt",row.names = FALSE , col.names = TRUE, quote= FALSE)






##Clumping Split 2 -----
#keeping only SNP with pv below 10^-3
#SNP s for propensity score (PS)
#rm(list=ls())
library(data.table)

Split_2_GWAS_mother <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_2_GWAS_Maternal/Split_2_GWAS_mother.txt")

Split_2_GWAS_mother$P <- as.numeric(Split_2_GWAS_mother$P)
jpeg(file="/mnt/cargo/william.denault/graph_CFMR/mh2.jpg", width = 850, height = 650)
manhattan.plot(chr= Split_2_GWAS_mother$`#CHROM`,
               pos= Split_2_GWAS_mother$POS,
               pvalue= Split_2_GWAS_mother$P)
abline(a=-log10(5*10^(-8)), b=0)
dev.off()
SNP_PS <- Split_2_GWAS_mother[ -which(Split_2_GWAS_mother$P  >10^(-2)),]

if( length(which(SNP_PS$ID == ".")>0))#to deal without ID
{
  for( i in which(SNP_PS$ID == "."))
  {
    SNP_PS$ID[i] <-  paste(SNP_PS$`#CHROM`[i],SNP_PS$POS[i], sep=":" )
  }
}

tt <-  SNP_PS[order(SNP_PS$P),]
print("Starting Clumping")

thresh <- 500000#max distancez between two SNPs in the PS
maxit  <- dim(tt)[1]


for( i in 1:maxit)
{
  temp <- which(tt$`#CHROM` == tt$`#CHROM`[i] )
  if( length(temp)>1)
  {
    tempd     <- abs( tt$POS - tt$POS[i] )
    tempd     <- which(tempd < thresh) #index SNPS to close
    to_remove <- temp[which(temp%in%tempd)]#the SNP on the same chromosome being to close to the lead SNP
    if ( length(to_remove)>1)#1 because the leading SNP is included
    {
      to_remove <- to_remove[ -which(to_remove ==i)]
      tt <- tt[-to_remove]
    }
  }
  if(dim(tt)[1]==i)
  {
    break()
  }
}

SNP_PS <- tt
head(SNP_PS)
dim(SNP_PS) #for two splits 1518 SNPS
##plot( -log10(SNP_PS$P), main=2)
write.table(SNP_PS, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Split_2_SNP_for_extraction.txt",row.names = FALSE , col.names = TRUE, quote= FALSE)




##Clumping Split 3 -----
#keeping only SNP with pv below 10^-3
#SNP s for propensity score (PS)
##rm(list=ls())
library(data.table)

Split_3_GWAS_mother <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_3_GWAS_Maternal/Split_3_GWAS_mother.txt")

jpeg(file="/mnt/cargo/william.denault/graph_CFMR/mh3.jpg", width = 850, height = 650)
manhattan.plot(chr= Split_3_GWAS_mother$`#CHROM`,
               pos= Split_3_GWAS_mother$POS,
               pvalue= Split_3_GWAS_mother$P)
abline(a=-log10(5*10^(-8)), b=0)
dev.off()
SNP_PS <- Split_3_GWAS_mother[ -which(Split_3_GWAS_mother$P  >10^(-2)),]

if( length(which(SNP_PS$ID == ".")>0))#to deal without ID
{
  for( i in which(SNP_PS$ID == "."))
  {
    SNP_PS$ID[i] <-  paste(SNP_PS$`#CHROM`[i],SNP_PS$POS[i], sep=":" )
  }
}

tt <-  SNP_PS[order(SNP_PS$P),]
print("Starting Clumping")

thresh <- 500000#max distancez between two SNPs in the PS
maxit  <- dim(tt)[1]


for( i in 1:maxit)
{
  temp <- which(tt$`#CHROM` == tt$`#CHROM`[i] )
  if( length(temp)>1)
  {
    tempd     <- abs( tt$POS - tt$POS[i] )
    tempd     <- which(tempd < thresh) #index SNPS to close
    to_remove <- temp[which(temp%in%tempd)]#the SNP on the same chromosome being to close to the lead SNP
    if ( length(to_remove)>1)#1 because the leading SNP is included
    {
      to_remove <- to_remove[ -which(to_remove ==i)]
      tt <- tt[-to_remove]
    }
  }
  if(dim(tt)[1]==i)
  {
    break()
  }
}

SNP_PS <- tt
head(SNP_PS)
dim(SNP_PS) #for two splits 1518 SNPS
#plot( -log10(SNP_PS$P), main=3)
write.table(SNP_PS, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Split_3_SNP_for_extraction.txt",row.names = FALSE , col.names = TRUE, quote= FALSE)



##Clumping Split 4 -----
#keeping only SNP with pv below 10^-3
#SNP s for propensity score (PS)
#rm(list=ls())
library(data.table)

Split_4_GWAS_mother <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_4_GWAS_Maternal/Split_4_GWAS_mother.txt")
jpeg(file="/mnt/cargo/william.denault/graph_CFMR/mh4.jpg", width = 850, height = 650)
manhattan.plot(chr= Split_4_GWAS_mother$`#CHROM`,
               pos= Split_4_GWAS_mother$POS,
               pvalue= Split_4_GWAS_mother$P)
abline(a=-log10(5*10^(-8)), b=0)
dev.off()
SNP_PS <- Split_4_GWAS_mother[ -which(Split_4_GWAS_mother$P  >10^(-2)),]

if( length(which(SNP_PS$ID == ".")>0))#to deal without ID
{
  for( i in which(SNP_PS$ID == "."))
  {
    SNP_PS$ID[i] <-  paste(SNP_PS$`#CHROM`[i],SNP_PS$POS[i], sep=":" )
  }
}

tt <-  SNP_PS[order(SNP_PS$P),]
print("Starting Clumping")

thresh <- 500000#max distancez between two SNPs in the PS
maxit  <- dim(tt)[1]


for( i in 1:maxit)
{
  temp <- which(tt$`#CHROM` == tt$`#CHROM`[i] )
  if( length(temp)>1)
  {
    tempd     <- abs( tt$POS - tt$POS[i] )
    tempd     <- which(tempd < thresh) #index SNPS to close
    to_remove <- temp[which(temp%in%tempd)]#the SNP on the same chromosome being to close to the lead SNP
    if ( length(to_remove)>1)#1 because the leading SNP is included
    {
      to_remove <- to_remove[ -which(to_remove ==i)]
      tt <- tt[-to_remove]
    }
  }
  if(dim(tt)[1]==i)
  {
    break()
  }
}

SNP_PS <- tt
head(SNP_PS)
dim(SNP_PS) #for two splits 1518 SNPS


write.table(SNP_PS, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Split_4_SNP_for_extraction.txt",row.names = FALSE , col.names = TRUE, quote= FALSE)



##Clumping Split 5 -----
#keeping only SNP with pv below 10^-3
#SNP s for propensity score (PS)
#rm(list=ls())
library(data.table)

Split_5_GWAS_mother <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_5_GWAS_Maternal/Split_5_GWAS_mother.txt")
jpeg(file="/mnt/cargo/william.denault/graph_CFMR/mh5.jpg", width = 850, height = 650)
manhattan.plot(chr= Split_5_GWAS_mother$`#CHROM`,
               pos= Split_5_GWAS_mother$POS,
               pvalue= Split_5_GWAS_mother$P)
abline(a=-log10(5*10^(-8)), b=0)
dev.off()
SNP_PS <- Split_5_GWAS_mother[ -which(Split_5_GWAS_mother$P  >10^(-2)),]

if( length(which(SNP_PS$ID == ".")>0))#to deal without ID
{
  for( i in which(SNP_PS$ID == "."))
  {
    SNP_PS$ID[i] <-  paste(SNP_PS$`#CHROM`[i],SNP_PS$POS[i], sep=":" )
  }
}

tt <-  SNP_PS[order(SNP_PS$P),]
print("Starting Clumping")

thresh <- 500000#max distancez between two SNPs in the PS
maxit  <- dim(tt)[1]


for( i in 1:maxit)
{
  temp <- which(tt$`#CHROM` == tt$`#CHROM`[i] )
  if( length(temp)>1)
  {
    tempd     <- abs( tt$POS - tt$POS[i] )
    tempd     <- which(tempd < thresh) #index SNPS to close
    to_remove <- temp[which(temp%in%tempd)]#the SNP on the same chromosome being to close to the lead SNP
    if ( length(to_remove)>1)#1 because the leading SNP is included
    {
      to_remove <- to_remove[ -which(to_remove ==i)]
      tt <- tt[-to_remove]
    }
  }
  if(dim(tt)[1]==i)
  {
    break()
  }
}

SNP_PS <- tt
head(SNP_PS)
dim(SNP_PS) #for two splits 1518 SNPS


write.table(SNP_PS, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Split_5_SNP_for_extraction.txt",row.names = FALSE , col.names = TRUE, quote= FALSE)




##Clumping Split 6 -----
#keeping only SNP with pv below 10^-3
#SNP s for propensity score (PS)
#rm(list=ls())
library(data.table)

Split_6_GWAS_mother <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_6_GWAS_Maternal/Split_6_GWAS_mother.txt")
jpeg(file="/mnt/cargo/william.denault/graph_CFMR/mh6.jpg", width = 850, height = 650)
manhattan.plot(chr= Split_6_GWAS_mother$`#CHROM`,
               pos= Split_6_GWAS_mother$POS,
               pvalue= Split_6_GWAS_mother$P)
abline(a=-log10(5*10^(-8)), b=0)
dev.off()
SNP_PS <- Split_6_GWAS_mother[ -which(Split_6_GWAS_mother$P  >10^(-2)),]

if( length(which(SNP_PS$ID == ".")>0))#to deal without ID
{
  for( i in which(SNP_PS$ID == "."))
  {
    SNP_PS$ID[i] <-  paste(SNP_PS$`#CHROM`[i],SNP_PS$POS[i], sep=":" )
  }
}

tt <-  SNP_PS[order(SNP_PS$P),]
print("Starting Clumping")

thresh <- 500000#max distancez between two SNPs in the PS
maxit  <- dim(tt)[1]


for( i in 1:maxit)
{
  temp <- which(tt$`#CHROM` == tt$`#CHROM`[i] )
  if( length(temp)>1)
  {
    tempd     <- abs( tt$POS - tt$POS[i] )
    tempd     <- which(tempd < thresh) #index SNPS to close
    to_remove <- temp[which(temp%in%tempd)]#the SNP on the same chromosome being to close to the lead SNP
    if ( length(to_remove)>1)#1 because the leading SNP is included
    {
      to_remove <- to_remove[ -which(to_remove ==i)]
      tt <- tt[-to_remove]
    }
  }
  if(dim(tt)[1]==i)
  {
    break()
  }
}

SNP_PS <- tt
head(SNP_PS)
dim(SNP_PS) #for two splits 1518 SNPS

write.table(SNP_PS, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Split_6_SNP_for_extraction.txt",row.names = FALSE , col.names = TRUE, quote= FALSE)




##Clumping Split 7 -----
#keeping only SNP with pv below 10^-3
#SNP s for propensity score (PS)
#rm(list=ls())
library(data.table)

Split_7_GWAS_mother <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_7_GWAS_Maternal/Split_7_GWAS_mother.txt")

jpeg(file="/mnt/cargo/william.denault/graph_CFMR/mh7.jpg", width = 850, height = 650)
manhattan.plot(chr= Split_7_GWAS_mother$`#CHROM`,
               pos= Split_7_GWAS_mother$POS,
               pvalue= Split_7_GWAS_mother$P)
abline(a=-log10(5*10^(-8)), b=0)
dev.off()

SNP_PS <- Split_7_GWAS_mother[ -which(Split_7_GWAS_mother$P  >10^(-2)),]

if( length(which(SNP_PS$ID == ".")>0))#to deal without ID
{
  for( i in which(SNP_PS$ID == "."))
  {
    SNP_PS$ID[i] <-  paste(SNP_PS$`#CHROM`[i],SNP_PS$POS[i], sep=":" )
  }
}

tt <-  SNP_PS[order(SNP_PS$P),]
print("Starting Clumping")

thresh <- 500000#max distancez between two SNPs in the PS
maxit  <- dim(tt)[1]


for( i in 1:maxit)
{
  temp <- which(tt$`#CHROM` == tt$`#CHROM`[i] )
  if( length(temp)>1)
  {
    tempd     <- abs( tt$POS - tt$POS[i] )
    tempd     <- which(tempd < thresh) #index SNPS to close
    to_remove <- temp[which(temp%in%tempd)]#the SNP on the same chromosome being to close to the lead SNP
    if ( length(to_remove)>1)#1 because the leading SNP is included
    {
      to_remove <- to_remove[ -which(to_remove ==i)]
      tt <- tt[-to_remove]
    }
  }
  if(dim(tt)[1]==i)
  {
    break()
  }
}

SNP_PS <- tt
head(SNP_PS)
dim(SNP_PS) #for two splits 1518 SNPS

write.table(SNP_PS, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Split_7_SNP_for_extraction.txt",row.names = FALSE , col.names = TRUE, quote= FALSE)




##Clumping Split 8 -----
#keeping only SNP with pv below 10^-3
#SNP s for propensity score (PS)
#rm(list=ls())
library(data.table)

Split_8_GWAS_mother <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_8_GWAS_Maternal/Split_8_GWAS_mother.txt")
jpeg(file="/mnt/cargo/william.denault/graph_CFMR/mh8.jpg", width = 850, height = 650)
manhattan.plot(chr= Split_8_GWAS_mother$`#CHROM`,
               pos= Split_8_GWAS_mother$POS,
               pvalue= Split_8_GWAS_mother$P)
abline(a=-log10(5*10^(-8)), b=0)
dev.off()
SNP_PS <- Split_8_GWAS_mother[ -which(Split_8_GWAS_mother$P  >10^(-2)),]

if( length(which(SNP_PS$ID == ".")>0))#to deal without ID
{
  for( i in which(SNP_PS$ID == "."))
  {
    SNP_PS$ID[i] <-  paste(SNP_PS$`#CHROM`[i],SNP_PS$POS[i], sep=":" )
  }
}

tt <-  SNP_PS[order(SNP_PS$P),]
print("Starting Clumping")

thresh <- 500000#max distancez between two SNPs in the PS
maxit  <- dim(tt)[1]


for( i in 1:maxit)
{
  temp <- which(tt$`#CHROM` == tt$`#CHROM`[i] )
  if( length(temp)>1)
  {
    tempd     <- abs( tt$POS - tt$POS[i] )
    tempd     <- which(tempd < thresh) #index SNPS to close
    to_remove <- temp[which(temp%in%tempd)]#the SNP on the same chromosome being to close to the lead SNP
    if ( length(to_remove)>1)#1 because the leading SNP is included
    {
      to_remove <- to_remove[ -which(to_remove ==i)]
      tt <- tt[-to_remove]
    }
  }
  if(dim(tt)[1]==i)
  {
    break()
  }
}

SNP_PS <- tt
head(SNP_PS)
dim(SNP_PS) #for two splits 1518 SNPS

write.table(SNP_PS, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Split_8_SNP_for_extraction.txt",row.names = FALSE , col.names = TRUE, quote= FALSE)




##Clumping Split 9 -----
#keeping only SNP with pv below 10^-3
#SNP s for propensity score (PS)
#rm(list=ls())
library(data.table)

Split_9_GWAS_mother <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_9_GWAS_Maternal/Split_9_GWAS_mother.txt")
jpeg(file="/mnt/cargo/william.denault/graph_CFMR/mh9.jpg", width = 850, height = 650)
manhattan.plot(chr= Split_9_GWAS_mother$`#CHROM`,
               pos= Split_9_GWAS_mother$POS,
               pvalue= Split_9_GWAS_mother$P)
abline(a=-log10(5*10^(-8)), b=0)
dev.off()
SNP_PS <- Split_9_GWAS_mother[ -which(Split_9_GWAS_mother$P  >10^(-2)),]

if( length(which(SNP_PS$ID == ".")>0))#to deal without ID
{
  for( i in which(SNP_PS$ID == "."))
  {
    SNP_PS$ID[i] <-  paste(SNP_PS$`#CHROM`[i],SNP_PS$POS[i], sep=":" )
  }
}

tt <-  SNP_PS[order(SNP_PS$P),]
print("Starting Clumping")

thresh <- 500000#max distancez between two SNPs in the PS
maxit  <- dim(tt)[1]


for( i in 1:maxit)
{
  temp <- which(tt$`#CHROM` == tt$`#CHROM`[i] )
  if( length(temp)>1)
  {
    tempd     <- abs( tt$POS - tt$POS[i] )
    tempd     <- which(tempd < thresh) #index SNPS to close
    to_remove <- temp[which(temp%in%tempd)]#the SNP on the same chromosome being to close to the lead SNP
    if ( length(to_remove)>1)#1 because the leading SNP is included
    {
      to_remove <- to_remove[ -which(to_remove ==i)]
      tt <- tt[-to_remove]
    }
  }
  if(dim(tt)[1]==i)
  {
    break()
  }
}

SNP_PS <- tt
head(SNP_PS)
dim(SNP_PS) #for two splits 1518 SNPS

write.table(SNP_PS, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Split_9_SNP_for_extraction.txt",row.names = FALSE , col.names = TRUE, quote= FALSE)





##Clumping Split 10 -----
#keeping only SNP with pv below 10^-3
#SNP s for propensity score (PS)
#rm(list=ls())
library(data.table)

Split_10_GWAS_mother <- fread("/mnt/work/william.denault/CFMR_Analysis/Data/Split_10_GWAS_Maternal/Split_10_GWAS_mother.txt")
jpeg(file="/mnt/cargo/william.denault/graph_CFMR/mh10.jpg", width = 850, height = 650)
manhattan.plot(chr= Split_10_GWAS_mother$`#CHROM`,
               pos= Split_10_GWAS_mother$POS,
               pvalue= Split_10_GWAS_mother$P)
abline(a=-log10(5*10^(-8)), b=0)
dev.off()
SNP_PS <- Split_10_GWAS_mother[ -which(Split_10_GWAS_mother$P  >10^(-2)),]

if( length(which(SNP_PS$ID == ".")>0))#to deal without ID
{
  for( i in which(SNP_PS$ID == "."))
  {
    SNP_PS$ID[i] <-  paste(SNP_PS$`#CHROM`[i],SNP_PS$POS[i], sep=":" )
  }
}

tt <-  SNP_PS[order(SNP_PS$P),]
print("Starting Clumping")

thresh <- 500000#max distancez between two SNPs in the PS
maxit  <- dim(tt)[1]


for( i in 1:maxit)
{
  temp <- which(tt$`#CHROM` == tt$`#CHROM`[i] )
  if( length(temp)>1)
  {
    tempd     <- abs( tt$POS - tt$POS[i] )
    tempd     <- which(tempd < thresh) #index SNPS to close
    to_remove <- temp[which(temp%in%tempd)]#the SNP on the same chromosome being to close to the lead SNP
    if ( length(to_remove)>1)#1 because the leading SNP is included
    {
      to_remove <- to_remove[ -which(to_remove ==i)]
      tt <- tt[-to_remove]
    }
  }
  if(dim(tt)[1]==i)
  {
    break()
  }
}

SNP_PS <- tt
head(SNP_PS)
dim(SNP_PS) #for two splits 1518 SNPS

write.table(SNP_PS, file= "/mnt/work/william.denault/CFMR_Analysis/Data/Split_10_SNP_for_extraction.txt",row.names = FALSE , col.names = TRUE, quote= FALSE)


#to continue the pipe
#setwd(dir = "/mnt/work/william.denault/Causal_ART")