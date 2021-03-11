rm(list=ls())
library(glmnet)
library(parallel)
library( mvtnorm)
library(AER)

n <- 10000
#Number of selected SNPs
p <- 300
#Number of SNPs affecting the exposure
np_act <- 5 #sample at 0.3 MAF HWE
sub_split <- 10 #number of split

set.seed(1)


simu_cross_fit_IV_est <- function(n,beta0,h2 )
{
 
  #Genotype data
  G <- matrix(sample(c(0,1,2),
                     size= (n*p),
                     prob = c(0.7^2, 2*0.7*0.3, 0.3^2),
                     replace = TRUE),
              nrow= n)
  
  #Hidden confounding
  H <- rnorm(n)
  #correlated noise 
  noise <-rmvnorm(n, 
                  mean = rep(0, 2), 
                  sigma = matrix(nrow = 2,
                                 byrow = TRUE,
                                 c(1,0.2,0.2,1)
                  )
  )
  


  V <- noise[,2]
  
  U <- noise[,1]
  

  #rescaling for insuring h2 = h2 on average
  tt <- apply( G[, (1:np_act)],1,sum)#genotype effect
  if( h2==0)
  {
    resc <- 0
  }else{resc <- ((var(H)+var(V))*(h2/(1-h2)))/var(tt)
  
  
  }
  
  ZPI <- sqrt(resc)*tt 
  
  
  #X is compound of ZPI , H hidden confounding, and a noise that correlate with the noise of Y
  
  X <- ZPI + H + V
  
  #Y the exposure alos affected by H
  
  Y <- X + H + U
  
  build_IV_sub_sample <- function(i) #i in 0:(sub_split-1)
  {
    indx <- (1+i*(n/sub_split)):((i+1)*(n/sub_split))
    
    
    Gtemp <- G[-indx,]
    Xtemp <- X[-indx]
    
    lasso <-  cv.glmnet(x=Gtemp, y=Xtemp,alpha=1)
    return(predict(lasso, G[indx,], s = "lambda.min"))
  }
  #Cross fitting LASSO
  #res  <- mclapply(0:(sub_split-1), build_IV_sub_sample, mc.cores = 6)
  res  <-  lapply(0:(sub_split-1), build_IV_sub_sample )
  
  IV   <- do.call(c,res)#cross fitted instrument
  res2 <- summary(ivreg(Y~X|IV))$coef[2,]
  
  
  out <- c(res2,n,beta0, h2)
  return(out)
}

simu_cross_fit_IV_est(n=5000, beta=1,h2=0.15)


beta0 <- c(0,0.05,-0.05,0.08,-0.08)
h2 <- c(0.2)
n <- c( 1000, 2000, 3000, 4000, 5000, 6000,7000, 8000, 9000, 10000)
N         <- rep(n, each=5000)  
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N))
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_Comparison_Deng_etal.RData")


beta0 <- c(0,0.05,0.08)
h2 <- c(0.2,0.1,0.05,0.01,0)
n <- c( 1000,  5000,  10000)
N         <- rep(n, each=15000)  
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_weak_IV_consistency.RData")


set.seed(2)
 
beta0 <- c(0,0.05,0.08)
h2 <- c(0.2,0.1,0.05,0.01,0)
n <- c(    100000)
N         <- rep(n, each=1000)
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 4)
save(res, file = "simulation_CFMR_weak_IV_consistency_100K2.RData")
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 4)
save(res, file = "simulation_CFMR_weak_IV_consistency_100K3.RData")
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 4)
save(res, file = "simulation_CFMR_weak_IV_consistency_100K4.RData")
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 4)
save(res, file = "simulation_CFMR_weak_IV_consistency_100K5.RData")
set.seed(2)
#
beta0 <- c(0,0.05,0.08)
h2 <- c(0.2,0.1,0.05,0.01,0)
n <- c(    500000)
N         <- rep(n, each=1000)
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_weak_IV_consistency_500K1.RData")
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_weak_IV_consistency_500K2.RData")
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_weak_IV_consistency_500K3.RData")
set.seed(3)
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_weak_IV_consistency_500K4.RData")
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_weak_IV_consistency_500K5.RData")

beta0 <- c(0,0.05,0.08)
h2 <-  c(0.00001, 0.0001,0.001,0.01)
n <- c( 1000,  5000,  10000)
N         <- rep(n, each=12000)  
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
table(N,B,H2)
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_very_weak_IV1k5K10K.RData")



beta0 <- c(0,0.05,0.08)
h2 <-  c(0.00001, 0.0001,0.001,0.01)
n <- c( 50000)
N         <- rep(n, each=12000)    
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
table(N,B,H2)
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_very_weak_IV1k50K.RData")





beta0 <- c(0,0.05,0.08)
h2 <-  c(0.00001, 0.0001,0.001,0.01)
n <- c( 50000)
N         <- rep(n, each=3600)   
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
table(N,B,H2)
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_very_weak_IV50K.RData")



beta0 <- c(0,0.05,0.08)
h2 <-  c(0.00001, 0.0001,0.001,0.01)
n <- c( 100000)
N         <- rep(n, each=3600)  
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
table(N,B,H2)
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_very_weak_IV100K.RData")




beta0 <- c(0,0.05,0.08)
h2 <-  c(0.00001, 0.0001,0.001,0.01)
n <- c( 100000)
N         <- rep(n, each=6000)  
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
table(N,B,H2)
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_very_weak_IV100K1.RData")
beta0 <- c(0,0.05,0.08)
h2 <-  c(0.00001, 0.0001,0.001,0.01)
n <- c( 100000)
N         <- rep(n, each=6000)  
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
table(N,B,H2)
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_very_weak_IV100K2.RData")
beta0 <- c(0,0.05,0.08)
h2 <-  c(0.00001, 0.0001,0.001,0.01)
n <- c( 100000)
N         <- rep(n, each=6000)  
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
table(N,B,H2)
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_very_weak_IV100K3.RData")

beta0 <- c(0,0.05,0.08)
h2 <-  c(0.00001, 0.0001,0.001,0.01)
n <- c( 500000)
N         <- rep(n, each=6000)  
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
table(N,B,H2)
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_very_weak_IV500K1.RData")
beta0 <- c(0,0.05,0.08)
h2 <-  c(0.00001, 0.0001,0.001,0.01)
n <- c( 500000)
N         <- rep(n, each=6000)  
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
table(N,B,H2)
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_very_weak_IV500K2.RData")
beta0 <- c(0,0.05,0.08)
h2 <-  c(0.00001, 0.0001,0.001,0.01)
n <- c( 500000)
N         <- rep(n, each=6000)  
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N)/length(h2))
table(N,B,H2)
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 5)
save(res, file = "simulation_CFMR_very_weak_IV500K3.RData")