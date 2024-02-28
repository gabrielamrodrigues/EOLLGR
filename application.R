rm(list=ls(all=TRUE))

##'Packages

library(ggplot2)
library(readxl)
library(gamlss)
library(gamlss.cens)
library(caret)
library(mlbench)
library(randomForestSRC)
library(randomForest)
library(ggRandomForests)
library(knitr)
library(dplyr)
library(survival)
library(survminer)
library(ggRandomForests)
library(pec)


data1 <- read_excel("new_campinas_2020julho.xlsx")

source('EOLLGR_gamlss.R')


##' Kaplan-Meier


km1 <- survfit(Surv(times, cens) ~sex, data = data1)

figure5a <- ggsurvplot(km1, data = data1, 
                       palette = c("#CC0000",1),    
                       xlab='Time (days)',
                       pval = TRUE, pval.coord = c(0, 0.1),  
                       legend.title="Sex",
                       legend.labs = c("1", "0"),
                       legend.position = c(0.3,0.90),
                       risk.table = F)

figure5a$plot+theme_gray()+
  theme(legend.background = element_rect(fill = "transparent"),
        axis.title = element_text( size = (12)),
        axis.text = element_text(size=11),
        legend.position = c(0.8,0.85))


km2 <- survfit(Surv(times, cens) ~heart, data = data1)

figure5b <- ggsurvplot(km2, data = data1, 
                       palette = c("#CC0000",1),    
                       xlab='Time (days)',
                       pval = TRUE, pval.coord = c(0, 0.1),  
                       legend.title="Chronic \n cardiovascular",
                       legend.labs = c("1", "0"),
                       risk.table = F)


figure5b$plot+theme_gray()+
  theme(legend.background = element_rect(fill = "transparent"),
        axis.title = element_text( size = (12)),
        axis.text = element_text(size=11),
        legend.position = c(0.8,0.85))

km3 <- survfit(Surv(times, cens) ~diab, data = data1)

figure5c <- ggsurvplot(km3, data = data1, 
                       palette = c("#CC0000",1),    
                       xlab='Time (days)',
                       pval = TRUE, pval.coord = c(0, 0.1),  
                       legend.title="Diabetes",
                       legend.labs = c("1", "0"),
                       risk.table = F)

figure5c$plot+theme_gray()+
  theme(legend.background = element_rect(fill = "transparent"),
        axis.title = element_text( size = (12)),
        axis.text = element_text(size=11),
        legend.position = c(0.8,0.85))


km4 <- survfit(Surv(times, cens) ~neuro, data = data1)

figure5d <- ggsurvplot(km4, data = data1, 
                       palette = c("#CC0000",1),    
                       xlab='Time (days)',
                       pval = TRUE, pval.coord = c(0, 0.1),  
                       legend.title="Chronic \n neurological",
                       legend.labs = c("1", "0"),
                       risk.table = F)

figure5d$plot+theme_gray()+
  theme(legend.background = element_rect(fill = "transparent"),
        axis.title = element_text( size = (12)),
        axis.text = element_text(size=11),
        legend.position = c(0.8,0.85))


km7 <- survfit(Surv(times, cens) ~renal, data = data1)

figure5e <- ggsurvplot(km7, data = data1, 
                       palette = c("#CC0000",1),    
                       xlab='Time (days)',
                       pval = TRUE, pval.coord = c(120, 0.4),  
                       legend.title="Chronic \n Kidney",
                       legend.labs = c("1", "0"),
                       risk.table = F)


figure5e$plot+theme_gray()+
  theme(legend.background = element_rect(fill = "transparent"),
        axis.title = element_text( size = (12)),
        axis.text = element_text(size=11),
        legend.position = c(0.8,0.85))


km6 <- survfit(Surv(times, cens) ~asthma, data = data1)

figure5f <- ggsurvplot(km6, data = data1, 
                       palette = c("#CC0000",1),    
                       xlab='Time (days)',
                       pval = TRUE, pval.coord = c(0, 0.1),  
                       legend.title="Asthma",
                       legend.labs = c("1", "0"),
                       risk.table = F)

figure5f$plot+theme_gray()+
  theme(legend.background = element_rect(fill = "transparent"),
        axis.title = element_text( size = (12)),
        axis.text = element_text(size=11),
        legend.position = c(0.8,0.85))

# Null models

fit.gr.m0  <-gamlss(Surv(times,cens)~ 1,family=cens(EOLLGR)(),data=data1, c.crit=0.01,sigma.start = 1, n.cyc=2000,tau.start = 1,tau.fix = T, nu.start = 1,nu.fix = T)


fit.ollgr.m0  <-gamlss(Surv(times,cens)~ 1,family=cens(EOLLGR)(),data=data1, c.crit=0.01,sigma.start = fit.gr.m0$sigma.fv, n.cyc=2000, mu.start = fit.gr.m0$mu.fv, tau.start = 1,tau.fix = T)

fit.egr.m0  <-gamlss(Surv(times,cens)~ 1,family=cens(EOLLGR)(),data=data1, c.crit=0.01,sigma.start =  fit.gr.m0$sigma.fv, n.cyc=2000, mu.start =  fit.gr.m0$mu.fv, nu.start = 1,nu.fix = T)

fit.eollgr.m0  <- gamlss(formula = Surv(times, cens) ~ 1, family = cens(EOLLGR)(), data = data1, mu.start = fit.egr.m0$mu.fv,sigma.start =fit.egr.m0$sigma.fv, c.crit = 0.01, n.cyc = 2000)


# Covariate selection

fit.gr.m1 <- stepGAICAll.A(fit.gr.m0, scope=list(lower=~1,upper=~sex+age+heart+diab+neuro+renal+asthma), sigma.try = F,nu.try = FALSE, tau.try = FALSE)

fit.egr.m1 <- stepGAICAll.A(fit.egr.m0, scope=list(lower=~1,upper=~sex+age+heart+diab+neuro+renal+asthma), nu.try = FALSE, tau.try = FALSE, sigma.try = F)

fit.ollgr.m1 <- stepGAICAll.A(fit.ollgr.m0, scope=list(lower=~1,upper=~sex+age+heart+diab+neuro+renal+asthma), nu.try = FALSE, tau.try = FALSE, sigma.try = F)

fit.eollgr.m1 <- stepGAICAll.A(fit.eollgr.m0, scope=list(lower=~1,upper=~sex+age+heart+diab+neuro+renal+asthma), nu.try = FALSE, tau.try = FALSE, sigma.try = F)

##Results

table3 <- rbind(c(AIC(fit.eollgr.m1),AIC(fit.ollgr.m1),AIC(fit.egr.m1),AIC(fit.gr.m1)),
                      
                      c(BIC(fit.eollgr.m1),BIC(fit.ollgr.m1),BIC(fit.egr.m1),BIC(fit.gr.m1)),
                      
                      c(fit.eollgr.m1$G.deviance,fit.ollgr.m1$G.deviance,fit.egr.m1$G.deviance,fit.gr.m1$G.deviance)
)

table4 <- summary(fit.eollgr.m1, type='qr')


# Residuals

figure6a <- wp(fit.eollgr.m1)
figure6b <- wp(fit.egr.m1)
figure6c <- wp(fit.ollgr.m1)
figure6d <- wp(fit.gr.m1)



# Envelope

fit.eollgr.m1 <- fit.eollgr.m1
n <- fit.eollgr.m1$N
mu <-     fit.eollgr.m1$mu.fv;
sigma <-   fit.eollgr.m1$sigma.fv; 
nu <- fit.eollgr.m1$nu.fv
tau <- fit.eollgr.m1$tau.fv

set.seed(1317)
res <- matrix(0,n,100)
erd <- matrix(0,n,100)
i <- 0
while(i <= 100){
  resp <- c()
  for(j in 1:n){
    resp[j] <- rEOLLGR(n=1,mu =mu[j],sigma = sigma[1],nu=nu[1], tau=tau[1])}
  fit <- try(gamlss(formula = Surv(resp, cens) ~age + renal + neuro, family = cens(EOLLGR)(), data = data1, mu.start = mu,sigma.start =sigma, tau.start = tau,nu.start=nu, c.crit = 0.01, n.cyc = 2000))
  if((class(fit)[1] != "try-error")==T){
    res[,i] <- sort(fit$residuals)
    i <- i+1
  }else{i <- i}
}


res.q <- fit.eollgr.m1$residuals

e1<-numeric(n)
e2<-numeric(n)


for(i in 1:n){
  eo<-sort(res[i,])
  e1[i]<-eo[1]
  e2[i]<-eo[100]
}


theoretical.quant <-  qnorm(1:n/(n+1))
sample.quant <-  sort(res.q)
db <-  data.frame(theoretical.quant, sample.quant)
df <- db %>%
  mutate(
    points1 = case_when(
      sample.quant < e1 ~ "out1",
      sample.quant > e2 ~ "out1",
      T ~ "within1"
    )
  )

db1 <-  data.frame(e1=sort(e1), theoretical.quant,points1=df$points1)
db2 <-  data.frame(e2=sort(e2), theoretical.quant,points1=df$points1)
media=colMeans(rbind(e1,e2))
db3 <- data.frame(media=sort(media),theoretical.quant)
out1 <- df[df$points1 == 'out1',]
within1 <- df[df$points1 == 'within1',]
obsout1 <- nrow(out1)
porc <- paste(obsout1," (", round((100*obsout1)/n,2),"%", ")",sep="")
percentage <- paste("Points out of envelope: ", porc, sep='')
totalp <- paste("Total points: ", 195, sep='')


g2 <-   ggplot(df , aes(x=theoretical.quant, y=sample.quant)) + geom_point()+xlab('Theoretical quantile')+ylab('Sample quantile')
figure7a <-  g2 + geom_line(aes(y=e1),db1)+geom_line(aes(y=e2),db2)+geom_line(aes(y=media),db3,linetype = "dashed",color='black',size=0.7)+annotate("text", x=0.5, y=-3.5,size = 4.5, label= percentage)


index <- seq(1:length(res.q))
d1r <- data.frame(index,res.q)
figure7b <- ggplot(d1r, aes(x=index, y=res.q)) + ylim(-4,4)+
  geom_point()+ylab('Quantile residual')+
  xlab('Index')+
  geom_hline(yintercept = 3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = -3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = 0,color='darkgray', size=0.8)




## marginally

fit.eollgr.renal  <- gamlss(formula = Surv(times, cens) ~renal , family = cens(EOLLGR)(), data = data1, mu.start = fit.gr.m0$mu.fv,sigma.start =fit.gr.m0$sigma.fv, tau.start = fit.gr.m0$tau.fv, c.crit = 0.01, n.cyc = 2000)

fit.eollgr.neuro  <- gamlss(formula = Surv(times, cens) ~ neuro, family = cens(EOLLGR)(), data = data1, mu.start = fit.gr.m0$mu.fv,sigma.start =fit.gr.m0$sigma.fv, tau.start = fit.gr.m0$tau.fv, c.crit = 0.001, n.cyc = 2000)

# cross validation

k <- 5
set.seed(1317)

folds <- createFolds(data1$times, k = k)


## regression models

model.eollgr.m1 <- list()
model.ollgr.m1 <- list()
model.egr.m1 <- list()
model.gr.m1 <- list()

set.seed(1311)
for (i in 1:k) {
  train <- data1[-folds[[i]],]
  test <- data1[folds[[i]],]
  
  fit.gr.m0  <-gamlss(Surv(times,cens)~ 1,family=cens(EOLLGR)(),data=train, c.crit=0.01, n.cyc=2000, tau.start = 1,tau.fix = T,nu.start = 1,nu.fix = T)
  
  fit.gr.m1  <-gamlss(Surv(times,cens)~ age + renal + neuro,family=cens(EOLLGR)(),data=train, c.crit=0.01, n.cyc=2000, tau.start = 1,tau.fix = T,nu.start = 1,nu.fix = T, mu.start = fit.gr.m0$mu.fv, sigma.start = fit.gr.m0$sigma.fv)
  
  fit.ollgr.m1  <-gamlss(Surv(times,cens)~ age + renal + neuro,family=cens(EOLLGR)(),data=train, c.crit=0.01,sigma.start = fit.gr.m1$sigma.fv, n.cyc=2000, mu.start = fit.gr.m1$mu.fv,tau.start = 1,tau.fix = T)
  
  fit.egr.m1  <-gamlss(Surv(times,cens)~ age + renal + neuro,family=cens(EOLLGR)(),data=train, c.crit=0.01, n.cyc=2000, mu.start = fit.gr.m1$mu.fv, sigma.start = fit.gr.m1$sigma.fv, nu.start = 1,nu.fix = T)
  
  fit.eollgr.m1  <- gamlss(formula = Surv(times, cens) ~age + renal + neuro, family = cens(EOLLGR)(), data = train, mu.start = fit.egr.m1$mu.fv,sigma.start =fit.egr.m1$sigma.fv, c.crit = 0.01, n.cyc = 2000)
  
  
  model.eollgr.m1.k05[[i]] <- fit.eollgr.m1
  model.ollgr.m1.k05[[i]] <- fit.ollgr.m1
  model.egr.m1.k05[[i]] <- fit.egr.m1
  model.gr.m1.k05[[i]] <- fit.gr.m1
  
  
}

## Results k=5


k <- 5

set.seed(1317)

folds <- createFolds(data1$times, k = k)
folds


aic.gr.k05 <- c()
bic.gr.k05 <- c()
gd.gr.k05 <- c()
c.gr.k05 <- c()

for(i in 1:k){
  train <- data1[-folds[[i]],]
  test <- data1[folds[[i]],]
  m1 <-  model.gr.m1.k05[[i]]
  
  aic.gr.k05[i] <- AIC(m1)
  bic.gr.k05[i] <- BIC(m1)
  gd.gr.k05[i] <- m1$G.deviance
  temp <- concordance(m1,newdata=test)
  c.gr.k05[i]  <- temp$concordance
}

aic.eollgr.k05 <- c()
bic.eollgr.k05 <- c()
gd.eollgr.k05 <- c()
c.eollgr.k05 <- c()

for(i in 1:k){
  train <- data1[-folds[[i]],]
  test <- data1[folds[[i]],]
  m1 <-  model.eollgr.m1.k05[[i]]
  
  aic.eollgr.k05[i] <- AIC(m1)
  bic.eollgr.k05[i] <- BIC(m1)
  gd.eollgr.k05[i] <- m1$G.deviance
  temp <- concordance(m1,newdata=test)
  c.eollgr.k05[i]  <- temp$concordance
}

aic.ollgr.k05 <- c()
bic.ollgr.k05 <- c()
gd.ollgr.k05 <- c()
c.ollgr.k05 <- c()

for(i in 1:k){
  train <- data1[-folds[[i]],]
  test <- data1[folds[[i]],]
  m1 <-  model.ollgr.m1.k05[[i]]
  
  aic.ollgr.k05[i] <- AIC(m1)
  bic.ollgr.k05[i] <- BIC(m1)
  gd.ollgr.k05[i] <- m1$G.deviance
  temp <- concordance(m1,newdata=test)
  c.ollgr.k05[i]  <- temp$concordance
}

aic.egr.k05 <- c()
bic.egr.k05 <- c()
gd.egr.k05 <- c()
c.egr.k05 <- c()

for(i in 1:k){
  train <- data1[-folds[[i]],]
  test <- data1[folds[[i]],]
  m1 <-  model.egr.m1.k05[[i]]
  
  aic.egr.k05[i] <- AIC(m1)
  bic.egr.k05[i] <- BIC(m1)
  gd.egr.k05[i] <- m1$G.deviance
  temp <- concordance(m1,newdata=test)
  c.egr.k05[i]  <- temp$concordance
}



reg.aic.k05 <- c(mean(aic.eollgr.k05),
                 mean(aic.ollgr.k05),
                 mean(aic.egr.k05),
                 mean(aic.gr.k05));reg.aic.k05

reg.bic.k05 <- c(mean(bic.eollgr.k05),
                 mean(bic.ollgr.k05),
                 mean(bic.egr.k05),
                 mean(bic.gr.k05));reg.bic.k05

reg.gd.k05 <- c(mean(gd.eollgr.k05),
                mean(gd.ollgr.k05),
                mean(gd.egr.k05),
                mean(gd.gr.k05));reg.gd.k05


reg.k05 <- c(mean(c.eollgr.k05),
             mean(c.ollgr.k05),
             mean(c.egr.k05),
             mean(c.gr.k05)
)


k <- 10
set.seed(1317)
folds <- createFolds(data1$times, k = k)

## regression models

model.eollgr.m1 <- list()
model.ollgr.m1 <- list()
model.egr.m1 <- list()
model.gr.m1 <- list()

set.seed(1311)
for (i in 1:k) {
  train <- data1[-folds[[i]],]
  test <- data1[folds[[i]],]
  
  
  fit.gr.m0  <-gamlss(Surv(times,cens)~ 1,family=cens(EOLLGR)(),data=train, c.crit=0.01, n.cyc=2000, tau.start = 1,tau.fix = T,nu.start = 1,nu.fix = T)
  
  fit.gr.m1  <-gamlss(Surv(times,cens)~ age + renal + neuro,family=cens(EOLLGR)(),data=train, c.crit=0.01, n.cyc=2000, tau.start = 1,tau.fix = T,nu.start = 1,nu.fix = T, mu.start = fit.gr.m0$mu.fv, sigma.start = fit.gr.m0$sigma.fv)
  
  fit.ollgr.m1  <-gamlss(Surv(times,cens)~ age + renal + neuro,family=cens(EOLLGR)(),data=train, c.crit=0.01,sigma.start = fit.gr.m1$sigma.fv, n.cyc=2000, mu.start = fit.gr.m1$mu.fv,tau.start = 1,tau.fix = T)
  
  fit.egr.m1  <-gamlss(Surv(times,cens)~ age + renal + neuro,family=cens(EOLLGR)(),data=train, c.crit=0.01, n.cyc=2000, mu.start = fit.gr.m1$mu.fv, sigma.start = fit.gr.m1$sigma.fv, nu.start = 1,nu.fix = T)
  
  fit.eollgr.m1  <- gamlss(formula = Surv(times, cens) ~age + renal + neuro, family = cens(EOLLGR)(), data = train, mu.start = fit.egr.m1$mu.fv,sigma.start =fit.egr.m1$sigma.fv, c.crit = 0.01, n.cyc = 2000)
  
  model.eollgr.m1.k10[[i]] <- fit.eollgr.m1
  model.ollgr.m1.k10[[i]] <- fit.ollgr.m1
  model.egr.m1.k10[[i]] <- fit.egr.m1
  model.gr.m1.k10[[i]] <- fit.gr.m1
  
  
}


## results k=10


k <- 10

set.seed(1317)

folds <- createFolds(data1$times, k = k)
folds

aic.gr.k10 <- c()
bic.gr.k10 <- c()
gd.gr.k10 <- c()
c.gr.k10 <- c()

for(i in 1:k){
  train <- data1[-folds[[i]],]
  test <- data1[folds[[i]],]
  m1 <-  model.gr.m1.k10[[i]]
  
  aic.gr.k10[i] <- AIC(m1)
  bic.gr.k10[i] <- BIC(m1)
  gd.gr.k10[i] <- m1$G.deviance
  temp <- concordance(m1,newdata=test)
  c.gr.k10[i]  <- temp$concordance
}

aic.eollgr.k10 <- c()
bic.eollgr.k10 <- c()
gd.eollgr.k10 <- c()
c.eollgr.k10 <- c()

for(i in 1:k){
  train <- data1[-folds[[i]],]
  test <- data1[folds[[i]],]
  m1 <-  model.eollgr.m1.k10[[i]]
  
  aic.eollgr.k10[i] <- AIC(m1)
  bic.eollgr.k10[i] <- BIC(m1)
  gd.eollgr.k10[i] <- m1$G.deviance
  temp <- concordance(m1,newdata=test)
  c.eollgr.k10[i]  <- temp$concordance
}

aic.ollgr.k10 <- c()
bic.ollgr.k10 <- c()
gd.ollgr.k10 <- c()
c.ollgr.k10 <- c()

for(i in 1:k){
  train <- data1[-folds[[i]],]
  test <- data1[folds[[i]],]
  m1 <-  model.ollgr.m1.k10[[i]]
  
  aic.ollgr.k10[i] <- AIC(m1)
  bic.ollgr.k10[i] <- BIC(m1)
  gd.ollgr.k10[i] <- m1$G.deviance
  temp <- concordance(m1,newdata=test)
  c.ollgr.k10[i]  <- temp$concordance
}

aic.egr.k10 <- c()
bic.egr.k10 <- c()
gd.egr.k10 <- c()
c.egr.k10 <- c()

for(i in 1:k){
  train <- data1[-folds[[i]],]
  test <- data1[folds[[i]],]
  m1 <-  model.egr.m1.k10[[i]]
  
  aic.egr.k10[i] <- AIC(m1)
  bic.egr.k10[i] <- BIC(m1)
  gd.egr.k10[i] <- m1$G.deviance
  temp <- concordance(m1,newdata=test)
  c.egr.k10[i]  <- temp$concordance
}


reg.aic.k10 <- c(mean(aic.eollgr.k10),
                 mean(aic.ollgr.k10),
                 mean(aic.egr.k10),
                 mean(aic.gr.k10));reg.aic.k10

reg.bic.k10 <- c(mean(bic.eollgr.k10),
                 mean(bic.ollgr.k10),
                 mean(bic.egr.k10),
                 mean(bic.gr.k10));reg.bic.k10

reg.gd.k10 <- c(mean(gd.eollgr.k10),
                mean(gd.ollgr.k10),
                mean(gd.egr.k10),
                mean(gd.gr.k10));reg.gd.k10


reg.k10 <- c(mean(c.eollgr.k10),
             mean(c.ollgr.k10),
             mean(c.egr.k10),
             mean(c.gr.k10)
)

##' Forest cross validation

NTREE = 1000 
NODESIZE = 10 
NSPLIT = 10 
IMPORTANCE = 'permute'

k <- 5

set.seed(1317)

folds <- createFolds(data1$times, k = k)

model.forest3 <- list()
model.forest2 <- list()

set.seed(1317)
for (i in 1:k) {
  train <- data1[-folds[[i]],]
  test <- data1[folds[[i]],]
  
  
  fit.forest2 <- rfsrc(Surv(times, cens)~sex+age+heart+diab+neuro+renal+asthma, train,  ntree = NTREE, 
                       nodesize = NODESIZE, 
                       mtry = 2, 
                       nsplit = NSPLIT, 
                       importance = IMPORTANCE, seed=1311) 
  
  model.forest2.k05[[i]] <- fit.forest2
  
}

#results k=5

f2.k05<- c()
imp2.k05 <- list()
for (i in 1:k) {
  model <- model.forest2.k05[[i]]
  test <- data1[folds[[i]],]
  
  pred <- predict(model, test,  na.action = "na.impute")
  f2.k05[i] <- 1-get.cindex(time = pred$yvar[,1], 
                         censoring = pred$yvar[,2], 
                         predicted = pred$predicted)
  imp2.k05[[i]] <- data.frame(gg_vimp(model))
}


md.k05 <- list()
for(i in 1:5){
  varsel_2 <- var.select(model.forest2.k05[[i]])
  gg_md <- gg_minimal_depth(varsel_2)
  md.k05[[i]] <- gg_md
}



forest1 <- model.forest2.k05[[1]]
forest2 <- model.forest2.k05[[2]]
forest3 <- model.forest2.k05[[3]]
forest4 <- model.forest2.k05[[4]]
forest5 <- model.forest2.k05[[5]]


data.error <- rbind(gg_error(forest1),gg_error(forest2),gg_error(forest3),gg_error(forest4),gg_error(forest5))

data.error.05 <- data.frame(k=rep(5,nrow(data.error)),forest=rep(c(1:5),each=100),data.error)


## k=10

k <- 10
set.seed(1317)
folds <- createFolds(data1$times, k = k)


model.forest2 <-  list()

set.seed(1311)
for (i in 1:k) {
  train <- data1[-folds[[i]],]
  test <- data1[folds[[i]],]
  

  fit.forest2 <- rfsrc(Surv(times, cens)~sex+age+heart+diab+neuro+renal+asthma, train,  ntree = NTREE, 
                       nodesize = NODESIZE, 
                       mtry = 2, 
                       nsplit = NSPLIT, 
                       importance = IMPORTANCE, seed=1311) 
  
  

  model.forest2.k10[[i]] <- fit.forest2
  
}

##results

f2.k10.2 <- c()
imp2.k10.2 <- list()
for (i in 1:k){
  model <- model.forest2.k10.2[[i]]
  test <- data1[folds[[i]],]
  
  pred <- predict(model, test,  na.action = "na.impute")
  f2.k10.2[i] <- 1-get.cindex(time = pred$yvar[,1], 
                         censoring = pred$yvar[,2], 
                         predicted = pred$predicted)
  imp2.k10.2[[i]] <- data.frame(gg_vimp(model))
}


md.k10 <- list()
for(i in 1:k){
  varsel_2 <- var.select(model.forest2.k10.2[[i]])
  gg_md <- gg_minimal_depth(varsel_2)
  md.k10[[i]] <- gg_md
}

forest1 <- model.forest2.k10[[1]]
forest2 <- model.forest2.k10[[2]]
forest3 <- model.forest2.k10[[3]]
forest4 <- model.forest2.k10[[4]]
forest5 <- model.forest2.k10[[5]]
forest6 <- model.forest2.k10[[6]]
forest7 <- model.forest2.k10[[7]]
forest8 <- model.forest2.k10[[8]]
forest9 <- model.forest2.k10[[9]]
forest10 <- model.forest2.k10[[10]]

data.error <- rbind(gg_error(forest1),gg_error(forest2),gg_error(forest3),gg_error(forest4),gg_error(forest5),gg_error(forest6),gg_error(forest7),gg_error(forest8),gg_error(forest9),gg_error(forest10))

data.error.10 <- data.frame(k=rep(10,nrow(data.error)),forest=rep(c(1:10),each=100),data.error)

##Tables

table5 <- rbind(reg1,for2)
rownames(table5) <- c('EOLLGR','OLLGR','EGR','GR','Forest')
colnames(table5) <- c('K=5','K=10')


table6 <- rbind(rbind(reg.aic.k05, reg.bic.k05,reg.gd.k05),
                 rbind(reg.aic.k10, reg.bic.k10,reg.gd.k10))

colnames(table6) <- c('EOLLGR','OLLGR','EGR','GR')


#Variable Importance

data.imp.k05 <- rbind(imp2.k05[[1]],imp2.k05[[2]],imp2.k05[[3]],imp2.k05[[4]],imp2.k05[[5]])


names <- as.character(data.imp.k05$vars)

names[names == 'age'] <- '1 age'
names[names == 'renal'] <- '2 renal'
names[names == 'neuro'] <- '4 neuro'
names[names == 'heart'] <- '3 heart'
names[names == 'diab'] <- '6 diab'
names[names == 'sex'] <- '7 sex'
names[names == 'asthma'] <- '5 asthma'
names

data.imp.k05$vars <- as.factor(names)


figure9a <- ggplot(data.imp.k05, aes(x=vars ,y=vimp)) +
  geom_boxplot(fill='#A4A4A4', color="black",
               outlier.size=2)+
  stat_summary(fun.y=mean, geom="point", size=2, color="red", fill="black")+
  labs(title="",x="", y = "Variable Importance")

figure9a+theme(legend.background = element_rect(fill = "transparent"),
              axis.title = element_text( size = (12)),
              axis.text = element_text(size=11),
              legend.position = c(0.8,0.75)) +
  theme(legend.key.width = unit(1.2,"cm"))+
  scale_x_discrete(labels = c('Idade', 'Renal','Cardio','Neuro','Asma','Diabete', 'Sexo'))+geom_hline(yintercept=0, linetype="dashed", color = "red")


data.imp.k10 <- rbind(imp2.k10[[1]],imp2.k10[[2]],imp2.k10[[3]],imp2.k10[[4]],imp2.k10[[5]],imp2.k10[[6]],imp2.k10[[7]],imp2.k10[[8]],imp2.k10[[9]],imp2.k10[[10]])


names <- as.character(data.imp.k10$vars)

names[names == 'age'] <- '1 age'
names[names == 'renal'] <- '2 renal'
names[names == 'neuro'] <- '3 neuro'
names[names == 'heart'] <- '4 heart'
names[names == 'diab'] <- '6 diab'
names[names == 'sex'] <- '7 sex'
names[names == 'asthma'] <- '5 asthma'
names

data.imp.k10$vars <- as.factor(names)


figure9b <- ggplot(data.imp.k10, aes(x=vars ,y=vimp)) +
  geom_boxplot(fill='#A4A4A4', color="black",
               outlier.size=2)+
  stat_summary(fun.y=mean, geom="point", size=2, color="red", fill="black")+
  labs(title="",x="", y = "Variable Importance")

figure9b+theme(legend.background = element_rect(fill = "transparent"),
              axis.title = element_text( size = (12)),
              axis.text = element_text(size=11),
              legend.position = c(0.8,0.75)) +
  theme(legend.key.width = unit(1.2,"cm"))+
  scale_x_discrete(labels = c('Idade', 'Renal', 'Neuro','Cardio', 'Asma','Diabetes', 'Sexo'))+geom_hline(yintercept=0, linetype="dashed", color = "red")



#Minimal Depth 

data.md.k05 <- rbind(md.k05[[1]]$varselect,md.k05[[2]]$varselect,md.k05[[3]]$varselect,md.k05[[4]]$varselect,md.k05[[5]]$varselect)


data.md.k05 <- data.frame(k=as.factor(rep(seq(1:5), each=7)),data.md.k05)
colnames(data.md.k05) <- c('k','depth','vimp', 'vars')


names <- as.character(data.md.k05$vars)

names[names == 'age'] <- '1 age'
names[names == 'renal'] <- '2 renal'
names[names == 'neuro'] <- '3 neuro'
names[names == 'heart'] <- '4 heart'
names[names == 'diab'] <- '5 diab'
names[names == 'sex'] <- '6 sex'
names[names == 'asthma'] <- '7 asthma'
names

data.md.k05$vars <- as.factor(names)


figure10a <- ggplot(data.md.k05, aes(x=vars ,y=depth)) +
  geom_boxplot(fill='#A4A4A4', color="black",
               outlier.size=2)+
  stat_summary(fun.y=mean, geom="point", size=2, color="black", fill="black")+
  labs(title="",x="", y = "Minimal Depth of a Variable")

figure10a+theme(legend.background = element_rect(fill = "transparent"),
                axis.title = element_text( size = (12)),
                axis.text = element_text(size=11),
                legend.position = c(0.8,0.75)) +
  theme(legend.key.width = unit(1.2,"cm"))+
  scale_x_discrete(labels = c('Age', 'Kidney', 'Neuro', 'Cardio', 'Diabetes', 'Sex', 'Asthma'))




data.md.k10 <- rbind(md.k10[[1]]$varselect,md.k10[[2]]$varselect,md.k10[[3]]$varselect,md.k10[[4]]$varselect,md.k10[[5]]$varselect,md.k10[[6]]$varselect,md.k10[[7]]$varselect,md.k10[[8]]$varselect,md.k10[[9]]$varselect,md.k10[[10]]$varselect)


data.md.k10 <- data.frame(k=as.factor(rep(seq(1:10), each=7)),data.md.k10)

colnames(data.md.k10) <- c('k','depth','vimp', 'vars')

names <- as.character(data.md.k10$vars)

names[names == 'age'] <- '1 age'
names[names == 'renal'] <- '2 renal'
names[names == 'neuro'] <- '3 neuro'
names[names == 'heart'] <- '4 heart'
names[names == 'diab'] <- '5 diab'
names[names == 'sex'] <- '6 sex'
names[names == 'asthma'] <- '7 asthma'
names

data.md.k10$vars <- as.factor(names)

figure10b <- ggplot(data.md.k10, aes(x=vars ,y=depth)) +
  geom_boxplot(fill='#A4A4A4', color="black",
               outlier.size=2)+
  stat_summary(fun.y=mean, geom="point", size=2, color="black", fill="black")+
  labs(title="",x="", y = "Minimal Depth of a Variable")

figure10b+theme(legend.background = element_rect(fill = "transparent"),
                axis.title = element_text( size = (12)),
                axis.text = element_text(size=11),
                legend.position = c(0.8,0.75)) +
  theme(legend.key.width = unit(1.2,"cm"))+
  scale_x_discrete(labels = c('Age', 'Kidney', 'Neuro', 'Cardio', 'Diabetes', 'Sex', 'Asthma'))

#Plot errors

p05e <- ggplot(data.error.05, aes(x=ntree, y=error, colour=as.factor(forest))) +
  geom_line(size=0.95)+labs(colour="k", y='Taxa de erro OOB', x="Número de árvores")

figure11a <- p05e+theme(
  legend.background = element_rect(fill = "transparent")
)+ylim(0.32,.40)+
  theme(axis.title = element_text( size = (12)),
        axis.text = element_text(size=11),
        legend.position = c(0.4,0.90), 
        legend.direction = "horizontal")


p10e <- ggplot(data.error.10, aes(x=ntree, y=error, colour=as.factor(forest))) +
  geom_line(size=0.95)+labs(colour="k", y='Taxa de erro OOB', x="Número de árvores")

figure11b <- p10e+theme(
  legend.background = element_rect(fill = "transparent")
)+ylim(0.32,.40)+
  theme(axis.title = element_text( size = (12)),
        axis.text = element_text(size=11),
        legend.position = c(0.5,0.90), 
        legend.direction = "horizontal")


library("pec")


k <- 5
NTREE = 1000 
NODESIZE = 10 
NSPLIT = 10 
IMPORTANCE = 'permute'

mod1 <- rfsrc(Surv(times, cens)~sex+age+heart+diab+neuro+renal+asthma, data1,  ntree = NTREE, 
              nodesize = NODESIZE, 
              mtry = 2, 
              nsplit = NSPLIT, 
              importance = IMPORTANCE, seed=1311) 


newdata1 <- rbind(c(age=60,renal=1, neuro=0),
                  c(age=60,renal=0, neuro=1),
                  c(age=60,renal=1, neuro=1),
                  c(age=60,renal=0, neuro=0))

newdata1 <- data.frame(newdata1)
newdata1$sex <- rep(NA,4)
newdata1$heart <- rep(NA,4)
newdata1$asthma  <- rep(NA,4)
newdata1$diab <- rep(NA,4)


sEOLLGR2 <- function(x,age, renal, neuro, sigma = 0.4, nu=1, tau=0.5){
  fit.eollgr.m1$mu.coefficients
  mu <- exp(6.47527669-0.01773443*age-0.48796167*renal-0.33235557 *neuro)
  s <- 1-pEOLLGR(x,mu,sigma,nu,tau)
  s
}


sigma <- fit.eollgr.m1$sigma.fv[1]
nu <- fit.eollgr.m1$nu.fv[1]
tau <- fit.eollgr.m1$tau.fv[1]

pp1.15 <- predictSurvProb(mod1, newdata = newdata1, times = 15,na.action='na.impute')


gr15 <- c(sEOLLGR2(15, 60,1,0,sigma = sigma, nu=nu, tau=tau),
          sEOLLGR2(15, 60,0,1,sigma = sigma, nu=nu, tau=tau),
          sEOLLGR2(15, 60,1,1,sigma = sigma, nu=nu, tau=tau),
          sEOLLGR2(15, 60,0,0,sigma = sigma, nu=nu, tau=tau))


pp1.30 <- predictSurvProb(mod1, newdata = newdata1, times = 30,na.action='na.impute')


gr30 <- c(sEOLLGR2(30, 60,1,0,sigma = sigma, nu=nu, tau=tau),
          sEOLLGR2(30, 60,0,1,sigma = sigma, nu=nu, tau=tau),
          sEOLLGR2(30, 60,1,1,sigma = sigma, nu=nu, tau=tau),
          sEOLLGR2(30, 60,0,0,sigma = sigma, nu=nu, tau=tau))


pp1.45 <- predictSurvProb(mod1, newdata = newdata1, times = 45,na.action='na.impute')


gr45 <- c(sEOLLGR2(45, 60,1,0,sigma = sigma, nu=nu, tau=tau),
          sEOLLGR2(45, 60,0,1,sigma = sigma, nu=nu, tau=tau),
          sEOLLGR2(45, 60,1,1,sigma = sigma, nu=nu, tau=tau),
          sEOLLGR2(45, 60,0,0,sigma = sigma, nu=nu, tau=tau))


table8 <- data.frame(pac=c(1,2,3,4),e=c(1,4) ,gr15,pp1.15,e=c(1,4) ,gr30,pp1.30,e=c(1,4) , gr45,pp1.45)



figure12a <- ggplot()+
  stat_function(fun = sEOLLGR2,args = list(60,1,0,sigma = sigma, nu=nu, tau=tau),size = 1, aes(color = "A"))+
  stat_function(fun = sEOLLGR2,args = list(60,0,1,sigma = sigma, nu=nu, tau=tau),size = 1, aes(color = "B"))+
  stat_function(fun = sEOLLGR2,args = list(60,1,1,sigma = sigma, nu=nu, tau=tau),size = 1, aes(color = "C"))+
  stat_function(fun = sEOLLGR2,args = list(60,0,0,sigma = sigma, nu=nu, tau=tau),size = 1, aes(color = "D"))+
  xlim(0,150)

figure12a+theme(legend.background = element_rect(fill = "transparent"),
                axis.title = element_text( size = (12)),
                axis.text = element_text(size=11),
                legend.position = c(0.8,0.8)) +
  labs(colour = "Patient")+ylab('Survival probability')+xlab('Time (days)')

pp1 <- predict.rfsrc(mod1, newdata = newdata1[1,],  na.action = "na.impute")
dados.pp1 <- data.frame(times=pp1$time.interest,times.n= as.numeric(pp1$survival))
pp2 <- predict.rfsrc(mod1, newdata = newdata1[2,],  na.action = "na.impute")
dados.pp2 <- data.frame(times=pp2$time.interest,times.n= as.numeric(pp2$survival))
pp3 <- predict.rfsrc(mod1, newdata = newdata1[3,],  na.action = "na.impute")
dados.pp3 <- data.frame(times=pp3$time.interest,times.n= as.numeric(pp3$survival))
pp4 <- predict.rfsrc(mod1, newdata = newdata1[4,],  na.action = "na.impute")
dados.pp4 <- data.frame(times=pp4$time.interest,times.n= as.numeric(pp4$survival))

figure12b <- ggplot()+
  geom_step(data=dados.pp1, aes(x=times, y = times.n,color = "A"),size=1)+
  geom_step(data=dados.pp2, aes(x=times, y = times.n,color = "B"),size=1)+
  geom_step(data=dados.pp3, aes(x=times, y = times.n,color = "C"),size=1)+
  geom_step(data=dados.pp4, aes(x=times, y = times.n,color = "D"),size=1)+
  xlim(0,150)

figure12b+theme(legend.background = element_rect(fill = "transparent"),
                axis.title = element_text( size = (12)),
                axis.text = element_text(size=11),
                legend.position = c(0.8,0.8)) +
  labs(colour = "Patient")+ylab('Survival probability')+xlab('Time (days)')






