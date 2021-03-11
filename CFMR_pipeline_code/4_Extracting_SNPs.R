rm(list=ls())

#Mother and the 10 splits ----
library(data.table)
setwd(dir = "/mnt/work/william.denault/CFMR_Analysis/Data")
i=1
temp <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))

for ( i in 2:10)
{
   if( i==1)
  {
    
    tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
    print(dim(tt))
  }else{
    
    tt <- fread(paste(getwd(),"/Split_",i,"_SNP_for_extraction.txt",sep=""))
    print(dim(tt))
  }
  to_rm <- which( tt$ID %in% temp$ID)
  if( length(to_rm)==0)
  {
    temp <- rbind(temp,tt)
  }  else{
    tt <- tt [-to_rm,]
    temp <- rbind(temp,tt)
  }
  
  
  
}


#extracting all the SNPs for all the split at the same time is faster than doing it split 
#by split
#We will write the split specific matrix afterwards



for ( i in 1:22)
{
  if( length(which(temp$`#CHROM`==i))>0)
    write.table(temp$ID[which(temp$`#CHROM`==i)],
                file ="/mnt/work/william.denault/CFMR_Analysis/Data/temp.txt",
                quote=FALSE,
                row.names = FALSE,
                col.names = FALSE)
  system(
    paste(
      "plink2 --vcf  /home/william.denault/archive/MERGE/vcf/",i,".vcf.gz",
      " --extract /mnt/work/william.denault/CFMR_Analysis/Data/temp.txt",
      " --export vcf",
      " --out /mnt/work/william.denault/CFMR_Analysis/Data/Extracted_SNPs_CHR",i,
      sep=""
    )
  )
}

