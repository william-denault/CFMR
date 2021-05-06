########Split 1 ----
for ( i in 1:22)
{
  break
  system(
    paste( "plink2 --vcf /home/william.denault/archive/MERGE/vcf/", i ,".vcf.gz",
           " --pheno /mnt/work/william.denault/CFMR_Analysis/Data/Split_1_Maternal_BMI.txt",
           " --linear /mnt/work/william.denault/CFMR_Analysis/Data/PC_1", 
           " --out /mnt/work/william.denault/CFMR_Analysis/Data/Split_1_GWAS_Maternal/Split_1_Maternal_Chr_",
           i,
           sep=""
    )
  )
  print(paste ("Split 1 Father CHR ",i ))
}

########Split 2 ----
for ( i in 1:22)
{
  break
  system(
    paste( "plink2 --vcf /home/william.denault/archive/MERGE/vcf/", i ,".vcf.gz",
           " --pheno /mnt/work/william.denault/CFMR_Analysis/Data/Split_2_Maternal_BMI.txt",
           " --linear /mnt/work/william.denault/CFMR_Analysis/Data/PC_2", 
           " --out /mnt/work/william.denault/CFMR_Analysis/Data/Split_2_GWAS_Maternal/Split_2_Maternal_Chr_",
           i,
           sep=""
    )
  )
  print(paste ("Split 2 Father CHR ",i ))
}



########Split 3 ----
for ( i in 1:22)
{
  break
  system(
    paste( "plink2 --vcf /home/william.denault/archive/MERGE/vcf/", i ,".vcf.gz",
           " --pheno /mnt/work/william.denault/CFMR_Analysis/Data/Split_3_Maternal_BMI.txt",
           " --linear /mnt/work/william.denault/CFMR_Analysis/Data/PC_3", 
           " --out /mnt/work/william.denault/CFMR_Analysis/Data/Split_3_GWAS_Maternal/Split_3_Maternal_Chr_",
           i,
           sep=""
    )
  )
  print(paste ("Split 3 Father CHR ",i ))
}



########Split 4 ----
for ( i in 1:22)
{
  break
  system(
    paste( "plink2 --vcf /home/william.denault/archive/MERGE/vcf/", i ,".vcf.gz",
           " --pheno /mnt/work/william.denault/CFMR_Analysis/Data/Split_4_Maternal_BMI.txt",
           " --linear /mnt/work/william.denault/CFMR_Analysis/Data/PC_4", 
           " --out /mnt/work/william.denault/CFMR_Analysis/Data/Split_4_GWAS_Maternal/Split_4_Maternal_Chr_",
           i,
           sep=""
    )
  )
  print(paste ("Split 4 Father CHR ",i ))
}

#######Split 5 ----
for ( i in 1:22)
{
  system(
    paste( "plink2 --vcf /home/william.denault/archive/MERGE/vcf/", i ,".vcf.gz",
           " --pheno /mnt/work/william.denault/CFMR_Analysis/Data/Split_5_Maternal_BMI.txt",
           " --linear /mnt/work/william.denault/CFMR_Analysis/Data/PC_5",
           " --out /mnt/work/william.denault/CFMR_Analysis/Data/Split_5_GWAS_Maternal/Split_5_Maternal_Chr_",
           i,
           sep=""
    )
  )
  print(paste ("Split 5 Father CHR ",i ))
}


########Split 6 ----
for ( i in 1:22)
{
  system(
    paste( "plink2 --vcf /home/william.denault/archive/MERGE/vcf/", i ,".vcf.gz",
           " --pheno /mnt/work/william.denault/CFMR_Analysis/Data/Split_6_Maternal_BMI.txt",
           " --linear /mnt/work/william.denault/CFMR_Analysis/Data/PC_6", 
           " --out /mnt/work/william.denault/CFMR_Analysis/Data/Split_6_GWAS_Maternal/Split_6_Maternal_Chr_",
           i,
           sep=""
    )
  )
  print(paste ("Split 6 Father CHR ",i ))
}


########Split 7 ----
for ( i in 1:22)
{
  system(
    paste( "plink2 --vcf /home/william.denault/archive/MERGE/vcf/", i ,".vcf.gz",
           " --pheno /mnt/work/william.denault/CFMR_Analysis/Data/Split_7_Maternal_BMI.txt",
           " --linear /mnt/work/william.denault/CFMR_Analysis/Data/PC_7", 
           " --out /mnt/work/william.denault/CFMR_Analysis/Data/Split_7_GWAS_Maternal/Split_7_Maternal_Chr_",
           i,
           sep=""
    )
  )
  print(paste ("Split 7 Father CHR ",i ))
}

########Split 8 ----
for ( i in 1:22)
{
  system(
    paste( "plink2 --vcf /home/william.denault/archive/MERGE/vcf/", i ,".vcf.gz",
           " --pheno /mnt/work/william.denault/CFMR_Analysis/Data/Split_8_Maternal_BMI.txt",
           " --linear /mnt/work/william.denault/CFMR_Analysis/Data/PC_8", 
           " --out /mnt/work/william.denault/CFMR_Analysis/Data/Split_8_GWAS_Maternal/Split_8_Maternal_Chr_",
           i,
           sep=""
    )
  )
  print(paste ("Split 8 Father CHR ",i ))
}

########Split 9 ----
for ( i in 1:22)
{
  system(
    paste( "plink2 --vcf /home/william.denault/archive/MERGE/vcf/", i ,".vcf.gz",
           " --pheno /mnt/work/william.denault/CFMR_Analysis/Data/Split_9_Maternal_BMI.txt",
           " --linear /mnt/work/william.denault/CFMR_Analysis/Data/PC_9",
           " --out /mnt/work/william.denault/CFMR_Analysis/Data/Split_9_GWAS_Maternal/Split_9_Maternal_Chr_",
           i,
           sep=""
    )
  )
  print(paste ("Split 9 Father CHR ",i ))
}


########Split 10 ----
for ( i in 1:22)
{
  system(
    paste( "plink2 --vcf /home/william.denault/archive/MERGE/vcf/", i ,".vcf.gz",
           " --pheno /mnt/work/william.denault/CFMR_Analysis/Data/Split_10_Maternal_BMI.txt",
           " --linear /mnt/work/william.denault/CFMR_Analysis/Data/PC_10",
           " --out /mnt/work/william.denault/CFMR_Analysis/Data/Split_10_GWAS_Maternal/Split_10_Maternal_Chr_",
           i,
           sep=""
    )
  )
  print(paste ("Split 10 Father CHR ",i ))
}


 
