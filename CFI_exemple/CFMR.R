rm(list=ls())
library(glmnet)
library(AER)
library( mvtnorm)

#function to perform the cross fitting
build_IV_sub_sample <- function(i) #i in 0:(sub_split-1)
{
  indx <- (1+i*(n/sub_split)):((i+1)*(n/sub_split))# here the splits are not random
  Gtemp <- G[-indx,]
  Xtemp <- X[-indx]
  lasso <-  cv.glmnet(x=Gtemp, y=Xtemp,alpha=1)
  return(predict(lasso, G[indx,], s = "lambda.min")
  )
}


set.seed(1)
#number of individuals
n <- 5000
#Number of selected SNPs
p <- 100
#Number of SNPs affecting the exposure
np_act <-5
sub_split <- 5#10 #number of split



#Y and X are sample with some correlated noise
noise <-rmvnorm(n, 
                mean = rep(0, 2), 
                sigma = matrix(nrow = 2,
                               byrow = TRUE,
                               c(1,0.7,0.7,1)
                )
)


G <- matrix(sample(c(0,1,2),
                   size= (n*p),
                   replace = TRUE),
            nrow= n
)
#Hidden confounding, increase sd to increase confouding
H <- rnorm(n,
           sd=1
          )
#ZPI is the part of X explained by the SNPs
ZPI <- apply( G[, (1:np_act)],
              1,
              mean
              )
#X is compound of ZPI , H hidden confounding, and a noise that correlate with the noise of Y
X <- ZPI+H+ noise[,2]
#Y the exposure alos affected by H
Y <- X+H+noise[,1] #Effect of X on Y is 1

cor(H,Y)



###Standard linear model
res0 <- lm(Y~X)
summary(res0) #obviously confouded effect, true effect is 1

#Identification and estimation of the insturement within the sample
lasso <-  cv.glmnet(x=G,
                    y=X,
                    alpha=1
)
IV <-  predict(lasso,
               G,
               s = "lambda.min"
)



res1 <- ivreg(Y~X|IV)
summary(res1) #biased estimate around 1.1

#Below CFMR estimation

res  <- lapply(0:(sub_split-1), build_IV_sub_sample)#performing cross-fitting


IV <- do.call(c,res)

#cross fitted IV, CFMR2
res2 <- ivreg(Y~X|IV)


summary(res2)#unbiased estimate

