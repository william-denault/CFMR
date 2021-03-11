rm(list=ls())
library(ggplot2)
load("~/DMLMR/simulation_CFMR_weak_IV_consistency.RData")
res1 <- res
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_100K.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_100K2.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_100K3.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_100K4.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_very_weak_IV100K1.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_very_weak_IV100K2.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_very_weak_IV100K3.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_100K5.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_50K1.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_50K2.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_50K3.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_50K4.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_50K5.RData")
res1 <- cbind(res1,res)
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_500K1.RData")


res <-  do.call( rbind,res[- which(do.call(c,(lapply(res,is.character))))])

res1 <- cbind(res1,t(res))
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_500K2.RData")


res <-  do.call( rbind,res[- which(do.call(c,(lapply(res,is.character))))])

res1 <- cbind(res1,t(res))
load("~/DMLMR/simulation_CFMR_weak_IV_consistency_500K3.RData")


res <-  do.call( rbind,res[- which(do.call(c,(lapply(res,is.character))))])


res1 <- cbind(res1,t(res))
res1 <- data.frame( t(res1))
res <- res1
colnames(res)[5:7] <- c("N", "beta","h2")
str(res)
res$power <- ifelse(res$Pr...t.. >0.05,0,1)

temp <- res[ which(res$beta==0.00),]
temp$N <- as.factor(temp$N)
temp <- temp[ -which(abs(temp$Estimate)>10),]
ggplot(temp[-which(temp$h2==0),], aes(x=Estimate, col=N))+ 
  geom_density()+
  geom_vline(xintercept = 0)+
  xlim(c(-0.2,0.2))+
  facet_wrap(~h2,scales = "free")+
  ggtitle("Density of the estimation of beta0 under H0\n for different  values of h2\n the vertical line correspond to beta=0")

ggsave("/mnt/cargo/william.denault/graph_CFMR/pvdensH0.pdf")
ggplot(temp[which(temp$h2==0),], aes(x=Estimate, col=N))+ 
  geom_density()+
  geom_vline(xintercept = 0)+
  ggtitle("Density of the estimation of beta0 under H0\n for different  h2=0")

ggsave("/mnt/cargo/william.denault/graph_CFMR/pvdensH0_h20.pdf")


temp <- res[ which(res$beta==0.05),]
temp$N <- as.factor(temp$N)
temp <- temp[ -which(abs(temp$Estimate)>10),]
ggplot(temp[-which(temp$h2==0),], aes(x=Estimate, col=N))+ 
  geom_density()+
  xlim(c(-0.25,0.35))+
  facet_wrap(~h2)+
  geom_vline(xintercept = 0.05)+
  ggtitle("Density of the estimation of beta0=0.05\n for different  values of h2")
ggsave("/mnt/cargo/william.denault/graph_CFMR/pvdensBeta005.pdf")
ggplot(temp[which(temp$h2==0),], aes(x=Estimate, col=N))+ 
  geom_density()+
  geom_vline(xintercept = 0.05)+
  ggtitle("Density of the estimation of beta0=0.05\n for of h2=0")
ggsave("/mnt/cargo/william.denault/graph_CFMR/pvdensBeta005h20.pdf")






temp <- res[ which(res$beta==0.08),]
temp$N <- as.factor(temp$N)
temp <- temp[ -which(abs(temp$Estimate)>10),]
ggplot(temp[-which(temp$h2==0),], aes(x=Estimate, col=N))+ 
  geom_density()+
  xlim(c(-0.12,0.22))+
  geom_vline(xintercept = 0.08)+
  facet_wrap(~h2)+
  ggtitle("Density of the estimation of beta0=0.08\n for different  values of h2")
ggsave("/mnt/cargo/william.denault/graph_CFMR/pvdensBeta008.pdf")

ggplot(temp[which(temp$h2==0),], aes(x=Estimate, col=N))+ 
  geom_density()+
 # xlim(c(-0.12,0.22))+
  geom_vline(xintercept = 0.08)+
  #facet_wrap(~h2)+
  ggtitle("Density of the estimation of beta0=0.08\n for h2 =0")
ggsave("/mnt/cargo/william.denault/graph_CFMR/pvdensBeta008h20.pdf")





library(dplyr)
out <- res %>% 
  group_by(h2,N,beta)%>%
  summarise(Beta_Est = median(Estimate ),sd_emp=sd(Estimate ), sd_theo = mean(Std..Error), power=mean(power),n=n() )
out <- data.frame(out)
print(xtable(out[-which(out$beta ==0),], digits = c(0,5,0,2,4,2,2,3,3)),include.rownames = FALSE)
      
      
out$rel_bias <- (out$beta-out$Beta_Est)/out$beta

is.na(out) <- sapply(out, is.infinite)

plot( out$N, out$sd_emp)
eq = function(x){+2/(sqrt(x))}
lines(y=eq(1000:10000), x=c(1000:10000))



lol <- data.frame(out)
lol <- lol[which(lol$h2 >0.0),]
lol$h2 <- as.factor(lol$h2)
lol$beta<- as.factor(lol$beta)
ggplot( lol, aes(N, sd_emp,shape=beta, col=h2))+
  geom_point(size=2)+
  ylim(c(0,0.152))+
  stat_function(fun=eq)+
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                           colour = "grey"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',
                                        colour = "grey"))+
  scale_x_log10(limits=c(800,100000))+
  ylab("Mean of the estimated standard deviations")+
  ggtitle("Mean of the estimated standard deviations of CFMR \n for different value of beta, h2 and N \n each configuration was simulated 1000 times \n the solid line is the curve n= sigma/sqrt(n)")

ggsave("/mnt/cargo/william.denault/graph_CFMR/conv_sd.pdf")
  
ggplot( lol, aes(Beta_Est, beta, col=as.factor(h2)))+
  geom_point()+
  xlim(c(-0.01,0.09))+
  geom_abline(slope = 1, intercept = 0)
  
  


is.na(out) <- sapply(out, is.infinite)


out_temp <- out
out_temp$N <- as.factor(out_temp$N)
out_temp$beta <-as.factor(out_temp$beta)
out_temp <- out_temp[-which(out_temp$h2==0),]
ggplot( out_temp, aes(h2, rel_bias, shape=beta,col=N))+
  geom_point()+
geom_abline(intercept = 0, slope = 0)+
  ylim(c(-0.3,0.3))
  
out_temp$N <-as.numeric(out_temp$N)
ggplot( out_temp, aes(N, rel_bias, shape=beta,col=as.factor(h2)))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE)+
  geom_abline(intercept = 0, slope = 0)+
  ylim(c(-0.3,0.3))

library(ggplot2)
temp <- out
temp$N <- as.factor(temp$N)
temp$h2 <- as.factor(temp$h2)
temp$beta<- as.factor(temp$beta)
ggplot( temp, aes(log(sd_emp),log(sd_theo),  colour= h2, shape= beta))+
  geom_point(size=4)+
  geom_abline(slope = 1, intercept = 0)+
  ylab("Log mean of the estimated standard deviations")+
  xlab("Log empricial standard deviations ")+
  facet_wrap(~N)+
  ggtitle("Empricial standard deviation of CFMR against mean of the estimated standard deviations \n of CFMR for different value of beta, h2 and N \n each configuration was simulated 1000 times")

ggsave("/mnt/cargo/william.denault/graph_CFMR/sd_to_large_panel.pdf")


temp <- out
temp$N <- as.factor(temp$N)
temp$h2 <- as.factor(temp$h2)
ggplot( temp, aes(beta, Beta_Est, colour= h2))+
  geom_point(size=4)+
  facet_wrap(~N)+
  xlab("Theoretical beta")+
  ylab("Mean estimate of beta")+
  geom_abline(slope = 1, intercept = 0)+
  ggtitle("Mean estimate of beta by CFMR against true beta \n for different value of beta, h2 and N")

ggsave("/mnt/cargo/william.denault/graph_CFMR/mean_est_panel.pdf")




library(xtable)


source("Deng_power.R")
out <- data.frame(out)


xtable(out)
out$type <- rep("CFMR", dim(out)[1])
Deng_Po$type <- rep("TSMR", dim(Deng_Po)[1])

ggplot( data.frame(rbind(out[, c("N","beta","power","type")],Deng_Po)),
        aes(N, power, colour= as.factor(beta), linetype=as.factor(type)))+
  geom_line(size=1.5)




res$power0.05 <- ifelse(res$Pr...t.. >0.05,0,1)
res$power0.01 <- ifelse(res$Pr...t.. >0.01,0,1)
res$power0.001 <- ifelse(res$Pr...t.. >0.001,0,1)

library(dplyr)
outb <- res %>% 
  group_by(h2,N,beta)%>%
  summarise( power0.05=mean(power0.05),power0.01=mean(power0.01),power0.001=mean(power0.001),n=n() )
data.frame(outb)


typeI_error <- outb[which(outb$beta==0),]
typeI_error
xtable(typeI_error,digits=3)

print(xtable(typeI_error[,-3],digits=c(0,5,0,4,4,4,0)), include.rownames = FALSE)
