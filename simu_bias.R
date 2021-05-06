
simu_cross_fit_IV_est <- function(n,p,np_act, h2, beta, cov=0.3 )
{
  if( missing(h2))
  {
    h2 <- 0.2
  }
  
  if( missing(beta))
  {
    beta <- 1
  }
  if(np_act >p)
  {
    print( "np_act should be smaller than p")
    break
  }
  n <- n
  #Number of selected SNPs
  p <- p
  #Number of SNPs affecting the exposure
  np_act <- np_act
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
                                 c(1,cov,cov,1)
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
  
  
  #X is compound of ZPI , H hidden confounding, and a noise that correlates with the noise of Y
  
  X <- ZPI + H + V
  
  #Y the exposure also affected by H
  
  Y <- beta*X + 40*H + U
  
  build_IV_sub_sample <- function(i) 
  {
    indx <- (1+i*(n/sub_split)):((i+1)*(n/sub_split))
    
    
    Gtemp <- G[-indx,]
    Xtemp <- X[-indx]
    
    lasso <-  cv.glmnet(x=Gtemp, y=Xtemp,alpha=1)
    return(predict(lasso, G[indx,], s = "lambda.min"))
  }
  
  
  
  #LM estimate
  res1 <- summary(lm(Y~X))$coef[2,]
  #Cross fitting LASSO
  res  <- lapply(0:(sub_split-1), build_IV_sub_sample)
  IV   <- do.call(c,res)
  res2 <- summary(ivreg(Y~X|IV))$coef[2,]
  
  #One sample lasso IV
  lasso <-  cv.glmnet(x=G, y=X,alpha=1)
  IV <- predict(lasso, G, s = "lambda.min")
  res3 <- summary(ivreg(Y~X|IV))$coef[2,]
  
  
  out  <- data.frame(rbind(res1,res2,res3))
  out$type <- c("LM","CFI_LASSO","IV_LASSO")
  out <- c(out[,1],n,p,beta,h2,np_act, cov)
  names(out ) <- c("LM","CFI_LASSO","IV_LASSO", "n","p","beta","h2","np_act","cov")
  return(out)
}

sub_split=5
library(mvtnorm)
library(glmnet)
library(AER)

my_n      <- rep(2000,50)
my_p      <- rep(300 ,50)
my_beta   <- rep(-0.8 ,50)
my_h2     <- rep(0.1 ,50)
my_cov    <- rep(0.9,50)
my_np_act <- rep(5   ,50)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

mean(out[,2])
mean(out[,3])
my_n      <- rep(1000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.5 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.3 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
              simu_cross_fit_IV_est, 
              n=my_n,
              p=my_p, 
              beta=my_beta,
              h2= my_h2, 
              np_act=my_np_act,
              cov=my_cov
              ))

save("check_bias_1SMR_1K_cov0.3")


my_n      <- rep(1000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.5 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_1K_cov0.5")

my_n      <- rep(1000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.7 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_1K_cov0.7")


my_n      <- rep(1000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.9 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_1K_cov0.9")

my_n      <- rep(5000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.3 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_5K_cov0.3")


my_n      <- rep(5000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.5 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_5K_cov0.5")

my_n      <- rep(5000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.7 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_5K_cov0.7")


my_n      <- rep(5000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.9 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_5K_cov0.9")



my_n      <- rep(10000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.3 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_10K_cov0.3")


my_n      <- rep(10000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.5 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_10K_cov0.5")

my_n      <- rep(10000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.7 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_10K_cov0.7")


my_n      <- rep(10000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.9 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_10K_cov0.9")




my_n      <- rep(50000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.3 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_50K_cov0.3")


my_n      <- rep(50000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.5 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_50K_cov0.5")

my_n      <- rep(50000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.7 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_50K_cov0.7")


my_n      <- rep(50000,3000)
my_p      <- rep(100 ,3000)
my_beta   <- rep(0.8 ,3000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
my_cov    <- rep(0.9 ,3000)
my_np_act <- rep(5   ,3000)
out <- t(mapply(
  simu_cross_fit_IV_est, 
  n=my_n,
  p=my_p, 
  beta=my_beta,
  h2= my_h2, 
  np_act=my_np_act,
  cov=my_cov
))

save("check_bias_1SMR_50K_cov0.9")




my_h2     <- rep(0.1 ,1000)
my_h2     <- c(rep(0.01 ,1000),rep(0.05,1000), rep(0.1,1000))
out
mean(out[2,])
mean(out[3,])


out
mean(out[2,])
mean(out[3,])



