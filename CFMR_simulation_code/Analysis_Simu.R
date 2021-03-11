rm(list = ls())
load("~/DMLMR/simulation_CFMR_Comparison_Deng_etal.RData")
str(res)
res <- data.frame(t(res))
colnames(res)[5:7] <- c("N", "beta","h2")
str(res)
res$power <- ifelse(res$Pr...t.. >0.05,0,1)
library(ggplot2)
temp <- res[ which(res$beta==0),]
temp$N <- as.factor(temp$N)
ggplot(temp, aes(x=Estimate, col=as.factor(N)))+#, col= as.factor(res$N+res$beta)))+
  geom_density()+
  ggtitle("Density of estimation of beta0=0 by CFMR \n for h2=20% and different value of N")+
  


  


library(dplyr)
out <- res %>% 
  group_by(N,beta)%>%
  summarise(Beta_Est = mean(Estimate ),sd_emp=sd(Estimate ), sd_theo = mean(Std..Error), power=mean(power),n=n() )
data.frame(out)
plot( out$N, out$sd_emp)
eq = function(x){2/(sqrt(x))}
lines(y=eq(1000:10000), x=c(1000:10000))

xtable(out[-which(out$beta==0),])
library(ggplot2)
ggplot( out, aes(N, sd_emp))+
  geom_point()+
  stat_function(fun=eq)+
  xlim(c(900,10000))+
  geom_point(out,mapping = aes(N, sd_theo, colour="red"))+
  ylab("Empirical standard deviation")+
  guides(color=FALSE)+
  ggtitle("In black empirical standard deviation of the estimates per condition \n in red the mean theoretical standard error per condition\n (h2=20%)")
ggsave("/mnt/cargo/william.denault/graph_CFMR/cv_sd_h2_0.2.pdf")


library(xtable)
xtable(out[-which(out$beta==0),], digits=3)


out <- data.frame(out)
out <- data.frame(out)
out$rel_bias <- (out$beta-out$Beta_Est)/out$beta

out

xtable(out)
source("Deng_power.R")
out$type <- rep("CFMR", dim(out)[1])
Deng_Po$type <- rep("TSMR", dim(Deng_Po)[1])

temp <- data.frame(rbind(out[, c("N","beta","power","type")],Deng_Po))
temp$beta <- as.factor(temp$beta)

temp$type <- as.factor(temp$type)
ggplot( temp, aes(N, power, colour= beta, linetype=type))+
  geom_line(size=1.5)+
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',
                                        colour = "grey"))
ggsave("/mnt/cargo/william.denault/graph_CFMR/power_eq_CFMR_TSMR.pdf")



temp <- out
temp$N <- as.factor(temp$N)
ggplot( temp, aes(beta, Beta_Est, colour= N))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  xlab("Theoretical beta")+
  ylab("Mean estimate of beta")+
  ggtitle("Mean estimate of CFMR for different \n value of beta and N  h2=20%")
ggsave("/mnt/cargo/william.denault/graph_CFMR/precision_CFMR_1K_10K.pdf")


ggplot( temp, aes(sd_theo, sd_emp, colour= N))+
  geom_point()+
  xlab("Theoretical standard deviation")+
  ylab("Empricial standard deviation ")+
  geom_abline(slope = 1, intercept = 0)+
  ggtitle("Empricial standard deviation of CFMR against theoretical standard deviation \n of CFMR for different value of beta and N, h2=20%")
ggsave("/mnt/cargo/william.denault/graph_CFMR/precision_CFMR_1K_10K.pdf")



res$power0.05 <- ifelse(res$Pr...t.. >0.05,0,1)
res$power0.01 <- ifelse(res$Pr...t.. >0.01,0,1)
res$power0.001 <- ifelse(res$Pr...t.. >0.001,0,1)

library(dplyr)
outb <- res %>% 
  group_by(N,beta)%>%
  summarise( power0.05=mean(power0.05),power0.01=mean(power0.01),power0.001=mean(power0.001),n=n() )
data.frame(outb)


typeI_error <- outb[which(outb$beta==0),]
typeI_error
xtable(typeI_error, digits=3)
