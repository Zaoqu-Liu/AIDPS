rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(tibble)
library(dplyr)
library(CoxBoost)
library(survivalsvm)
library(survival)

load('PACA_AU_Array.rda')
load('TCGA_PAAD.rda')

est_data <- PACA_AU_Array
val_data <- TCGA_PAAD
pre_var <- colnames(est_data)[-c(1:3)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd <- val_data[,c('OS.time','OS',pre_var)]

###CoxBoost
seed <- 39
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)

rid <- names(coef(fit)[which(coef(fit)!=0)])

####Survivalsvm
est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd2 <- val_data[,c('OS.time','OS',rid)]

x1 <- as.matrix(est_dd2[,rid])
x2 <- as.matrix(Surv(est_dd2$OS.time,est_dd2$OS))

fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)

PACA_AU_Array <- cbind(est_dd2[,1:2],RS=as.numeric(predict(fit, est_dd2)$predicted))
TCGA_PAAD <- cbind(val_dd2[,1:2],RS=as.numeric(predict(fit, val_dd2)$predicted))

save(rid,fit,file = "AIDPS.rda")









