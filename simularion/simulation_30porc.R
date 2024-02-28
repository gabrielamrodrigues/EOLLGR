
rm(list = ls())
source('EOLLGR_gamlss.R')

library(gamlss.cens)
library(randomForestSRC)
library(randomForest)
library(ggRandomForests)

# Function to generate regression model with 2 covariates
reg_eollgr <- function(n,beta10,beta11,beta12,beta13,beta14,beta15,beta20,beta3,beta40){
  x1 <- rbinom(n,1,0.45)
  x2 <- rbinom(n,1,0.50)
  x3 <- rbinom(n,1,0.55)
  x4 <- rnorm(n,0,0.5)
  x5 <- rnorm(n,0,1)
  mu <- exp(beta10+beta11*x1+beta12*x2+beta13*x3+beta14*x4+beta15*x5)
  sigma <- exp(beta20)
  nu <- exp(beta30)
  tau <- exp(beta40)
  y <- rEOLLGR(n=n,mu=mu,sigma=sigma,nu=nu, tau=tau)
  return(list(y=y,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,mu=mu,sigma=sigma,nu=nu, tau=tau))
}

beta10 <- 1.5; beta11 <- -0.5; beta12 <- -0.4; beta13 <- -0.3;
beta14 <- -0.6; beta15 <- -0.2; 
beta20 <- -2
beta30 <- 0.8; beta40 <- 0.9


initial <- c(beta10,beta11,beta12,beta13,beta14,beta15,beta20,beta30,beta40)


# Generate censored data
estimativas_cens <- function(n,beta10,beta11,beta12,beta13,beta14,beta15,beta20,beta30,beta40,porc.cens,tau){
  reg <- reg_eollgr(n,beta10,beta11,beta12,beta13,beta14,beta15,beta20,beta30,beta40)
  tempos <- reg$y
  x1 <- reg$x1
  x2 <- reg$x2
  x3 <- reg$x3
  x4 <- reg$x4
  x5 <- reg$x5
  mu0 <- reg$mu
  sigma0 <- reg$sigma
  nu0 <- reg$nu
  tau0 <- reg$tau
  tempo.cens <- runif(n,0,tau)
  delta <- c()
  tempo <- c()
  k <- 0
  for(j in 1:n){
    if(tempo.cens[j] <= tempos[j] && k < round(n*porc.cens/100))
    {delta[j] <- 0
    tempo[j] <- tempo.cens[j]
    k <- k+1}
    else
    {delta[j] <- 1
    tempo[j] <- tempos[j]}
  }
  dados <- data.frame(times=tempo,cens=delta,x1,x2,x3,x4,x5,mu0,sigma0,nu0,tau0)
  return(dados)
}


###

#30% n=80

porc.cens <- 30
tau <- 1
npar <- 9
n <- 80
r <- 1000

theta <- matrix(0,r,npar)
res.rq = matrix(0,r,n)

md_forest <- list()
imp.comp <- list()
data.error.comp2 <- list()
data.error.comp3 <- list()
data.error.comp4 <- list()

se <- list()
se2 <- list()
i <- 1
iter <- 0
cens.vector <- c()
c.eollgr.test  <- c()
c.eollgr.train  <- c()
c.eollgr.complete <- c()


c.forest.train2  <- c()
c.forest.train3  <- c()
c.forest.train4 <- c()


set.seed(1311)
while(i<=r){
  
  
  data1 <- estimativas_cens(n,beta10,beta11,beta12,beta13,beta14,beta15,beta20,beta30,beta40,porc.cens,tau)
  
  mu1 <- exp(beta10+beta11*data1$x1+beta12*data1$x2+beta13*data1$x3+beta14*data1$x4+beta15*data1$x5)
  
  fit1=try(gamlss(Surv(times,cens)~x1+x2+x3+x4+x5,family = cens(EOLLGR),
                  c.crit=0.01,n.cyc=200,data=data1,sigma.start = data1$sigma0,nu.start= data1$nu0,tau.start =  data1$tau0,
                  mu.start = data1$mu0,trace=F))
  
  n1 <- sample(1:2,
               size = nrow(data1),
               replace = TRUE,
               prob=c(0.7, 0.3))
  
  train <- data1[n1==1,]
  test <- data1[n1==2,]
  
  mu0 <- data1$mu[n1==1]
  sigma0 <- data1$sigma[n1==1]
  nu0 <- data1$nu[n1==1]
  tau0 <- data1$tau[n1==1]
  
  fit2=try(gamlss(Surv(times,cens)~x1+x2+x3+x4+x5,family = cens(EOLLGR),
                  n.cyc=200,data=train,
                  c.crit=0.01,mu.start = train$mu0,sigma.start = train$sigma0,nu.start= train$nu0,tau.start =  train$tau0,
                  trace=F))
  
  fit.forest2 <- rfsrc(Surv(times, cens)~x1+x2+x3+x4+x5, train, ntree = 5000, 
                       mtry = 2,  
                       importance = 'permute', seed=1311) 
  
  fit.forest3 <- rfsrc(Surv(times, cens)~x1+x2+x3+x4+x5, train, ntree = 5000, 
                       mtry = 3,  
                       importance = 'permute', seed=1311)
  
  fit.forest4 <- rfsrc(Surv(times, cens)~x1+x2+x3+x4+x5, train, ntree = 5000, 
                       mtry = 4,  
                       importance = 'permute', seed=1311)
  
  if((class(fit1)[1] != "try-error" & class(fit2)[1] != "try-error")==T){
    teste <- fit1$converged
    teste2 <- fit2$converged
    
    if(teste == TRUE & teste2 == TRUE){
      betas <- c(mu=fit1$mu.coefficients ,sigma=fit1$sigma.coefficients,nu=fit1$nu.coefficients,tau=fit1$tau.coefficients)
      theta[i,] <-  betas
      res.rq[i,] = fit1$residuals
      se[[i]] <- summary(fit1,type='qr')
      se2[[i]] <- summary(fit2,type='qr')
      cens.vector[i] <- summary(as.factor(data1$cens))[1]
      
      
      temp <- concordance(fit2,newdata=train)
      c.eollgr.train[i]  <- temp$concordance
      
      
      pred <- predict(fit.forest2, train,  na.action = "na.impute")
      c.forest.train2[i] <- 1-get.cindex(time = pred$yvar[,1], 
                                         censoring = pred$yvar[,2], 
                                         predicted = pred$predicted)
      
      pred <- predict(fit.forest3, train,  na.action = "na.impute")
      c.forest.train3[i] <- 1-get.cindex(time = pred$yvar[,1], 
                                         censoring = pred$yvar[,2], 
                                         predicted = pred$predicted)
      
      pred <- predict(fit.forest4, train,  na.action = "na.impute")
      c.forest.train4[i] <- 1-get.cindex(time = pred$yvar[,1], 
                                         censoring = pred$yvar[,2], 
                                         predicted = pred$predicted)
      
      data.error.comp2[[i]] <- (gg_error(fit.forest2))
      data.error.comp3[[i]] <- (gg_error(fit.forest3))
      data.error.comp4[[i]] <- (gg_error(fit.forest4))
      
      i <- i+1
    }
    else{i <- i}
  }
  else{i <- i}
  cat("iteration = ", iter <- iter + 1, i,"\n")
}

theta80_cens30 <- theta

save(list = "theta80_cens30",
     file = "results/results_30/theta80_cens30.Rda",
     compress = "xz")

se80_cens30 <- se
p80_cens30 <- se2

save(list = "se80_cens30",
     file = "results/results_30/se80_cens30.Rda",
     compress = "xz")

save(list = "p80_cens30",
     file = "results/results_30/p80_cens30.Rda",
     compress = "xz")

cens.vector80.30 <- cens.vector

save(list = "cens.vector80.30",
     file = "results/results_30/cens.vector80.30.Rda",
     compress = "xz")

c.eollgr80.30 <- c.eollgr.train
c.forest80.30 <- data.frame(c.forest.train2,c.forest.train3,c.forest.train4)

save(list = "c.forest80.30",
     file = "results/results_30/c.forest80.30.Rda",
     compress = "xz")

save(list = "c.eollgr80.30",
     file = "results/results_30/c.eollgr80.30.Rda",
     compress = "xz")


error80_cens302 <- data.error.comp2
error80_cens303 <- data.error.comp3
error80_cens304 <- data.error.comp4

save(list = "error80_cens302",
     file = "results/results_30/error80_cens302.Rda",
     compress = "xz")

save(list = "error80_cens303",
     file = "results/results_30/error80_cens303.Rda",
     compress = "xz")

save(list = "error80_cens304",
     file = "results/results_30/error80_cens304.Rda",
     compress = "xz")


#30% n=150

porc.cens <- 30
tau <- 1
npar <- 9
n <- 150
r <- 1000

theta <- matrix(0,r,npar)
res.rq = matrix(0,r,n)

md_forest <- list()
imp.comp <- list()
data.error.comp2 <- list()
data.error.comp3 <- list()
data.error.comp4 <- list()

se <- list()
se2 <- list()
i <- 1
iter <- 0
cens.vector <- c()
c.eollgr.test  <- c()
c.eollgr.train  <- c()
c.eollgr.complete <- c()


c.forest.train2  <- c()
c.forest.train3  <- c()
c.forest.train4 <- c()


set.seed(1311)
while(i<=r){
  
  
  data1 <- estimativas_cens(n,beta10,beta11,beta12,beta13,beta14,beta15,beta20,beta30,beta40,porc.cens,tau)
  
  mu1 <- exp(beta10+beta11*data1$x1+beta12*data1$x2+beta13*data1$x3+beta14*data1$x4+beta15*data1$x5)
  
  fit1=try(gamlss(Surv(times,cens)~x1+x2+x3+x4+x5,family = cens(EOLLGR),
                  c.crit=0.01,n.cyc=200,data=data1,sigma.start = data1$sigma0,nu.start= data1$nu0,tau.start =  data1$tau0,
                  mu.start = data1$mu0,trace=F))
  
  n1 <- sample(1:2,
               size = nrow(data1),
               replace = TRUE,
               prob=c(0.7, 0.3))
  
  train <- data1[n1==1,]
  test <- data1[n1==2,]
  
  mu0 <- data1$mu[n1==1]
  sigma0 <- data1$sigma[n1==1]
  nu0 <- data1$nu[n1==1]
  tau0 <- data1$tau[n1==1]
  
  fit2=try(gamlss(Surv(times,cens)~x1+x2+x3+x4+x5,family = cens(EOLLGR),
                  n.cyc=200,data=train,
                  c.crit=0.01,mu.start = train$mu0,sigma.start = train$sigma0,nu.start= train$nu0,tau.start =  train$tau0,
                  trace=F))
  
  fit.forest2 <- rfsrc(Surv(times, cens)~x1+x2+x3+x4+x5, train, ntree = 2000, 
                       mtry = 2,  
                       importance = 'permute', seed=1311) 
  
  fit.forest3 <- rfsrc(Surv(times, cens)~x1+x2+x3+x4+x5, train, ntree = 2000, 
                       mtry = 3,  
                       importance = 'permute', seed=1311)
  
  fit.forest4 <- rfsrc(Surv(times, cens)~x1+x2+x3+x4+x5, train, ntree = 2000, 
                       mtry = 4,  
                       importance = 'permute', seed=1311)
  
  if((class(fit1)[1] != "try-error" & class(fit2)[1] != "try-error")==T){
    teste <- fit1$converged
    teste2 <- fit2$converged
    
    if(teste == TRUE & teste2 == TRUE){
      betas <- c(mu=fit1$mu.coefficients ,sigma=fit1$sigma.coefficients,nu=fit1$nu.coefficients,tau=fit1$tau.coefficients)
      theta[i,] <-  betas
      res.rq[i,] = fit1$residuals
      se[[i]] <- summary(fit1,type='qr')
      se2[[i]] <- summary(fit2,type='qr')
      cens.vector[i] <- summary(as.factor(data1$cens))[1]
      
      
      temp <- concordance(fit2,newdata=train)
      c.eollgr.train[i]  <- temp$concordance
      
      
      pred <- predict(fit.forest2, train,  na.action = "na.impute")
      c.forest.train2[i] <- 1-get.cindex(time = pred$yvar[,1], 
                                         censoring = pred$yvar[,2], 
                                         predicted = pred$predicted)
      
      pred <- predict(fit.forest3, train,  na.action = "na.impute")
      c.forest.train3[i] <- 1-get.cindex(time = pred$yvar[,1], 
                                         censoring = pred$yvar[,2], 
                                         predicted = pred$predicted)
      
      pred <- predict(fit.forest4, train,  na.action = "na.impute")
      c.forest.train4[i] <- 1-get.cindex(time = pred$yvar[,1], 
                                         censoring = pred$yvar[,2], 
                                         predicted = pred$predicted)
      
      data.error.comp2[[i]] <- (gg_error(fit.forest2))
      data.error.comp3[[i]] <- (gg_error(fit.forest3))
      data.error.comp4[[i]] <- (gg_error(fit.forest4))
      
      i <- i+1
    }
    else{i <- i}
  }
  else{i <- i}
  cat("iteration = ", iter <- iter + 1, i,"\n")
}

theta150_cens30 <- theta

save(list = "theta150_cens30",
     file = "results/results_30/theta150_cens30.Rda",
     compress = "xz")

se150_cens30 <- se
p150_cens30 <- se2

save(list = "se150_cens30",
     file = "results/results_30/se150_cens30.Rda",
     compress = "xz")

save(list = "p150_cens30",
     file = "results/results_30/p150_cens30.Rda",
     compress = "xz")

cens.vector150.30 <- cens.vector

save(list = "cens.vector150.30",
     file = "results/results_30/cens.vector150.30.Rda",
     compress = "xz")

c.eollgr150.30 <- c.eollgr.train
c.forest150.30 <- data.frame(c.forest.train2,c.forest.train3,c.forest.train4)

save(list = "c.forest150.30",
     file = "results/results_30/c.forest150.30.Rda",
     compress = "xz")

save(list = "c.eollgr150.30",
     file = "results/results_30/c.eollgr150.30.Rda",
     compress = "xz")


error150_cens302 <- data.error.comp2
error150_cens303 <- data.error.comp3
error150_cens304 <- data.error.comp4

save(list = "error150_cens302",
     file = "results/results_30/error150_cens302.Rda",
     compress = "xz")

save(list = "error150_cens303",
     file = "results/results_30/error150_cens303.Rda",
     compress = "xz")

save(list = "error150_cens304",
     file = "results/results_30/error150_cens304.Rda",
     compress = "xz")

#30% n=450

porc.cens <- 30
tau <- 1
npar <- 9
n <- 450
r <- 1000

theta <- matrix(0,r,npar)
res.rq = matrix(0,r,n)

md_forest <- list()
imp.comp <- list()
data.error.comp2 <- list()
data.error.comp3 <- list()
data.error.comp4 <- list()

se <- list()
se2 <- list()
i <- 1
iter <- 0
cens.vector <- c()
c.eollgr.test  <- c()
c.eollgr.train  <- c()
c.eollgr.complete <- c()


c.forest.train2  <- c()
c.forest.train3  <- c()
c.forest.train4 <- c()


set.seed(1311)
while(i<=r){
  
  
  data1 <- estimativas_cens(n,beta10,beta11,beta12,beta13,beta14,beta15,beta20,beta30,beta40,porc.cens,tau)
  
  mu1 <- exp(beta10+beta11*data1$x1+beta12*data1$x2+beta13*data1$x3+beta14*data1$x4+beta15*data1$x5)
  
  fit1=try(gamlss(Surv(times,cens)~x1+x2+x3+x4+x5,family = cens(EOLLGR),
                  c.crit=0.01,n.cyc=200,data=data1,sigma.start = data1$sigma0,nu.start= data1$nu0,tau.start =  data1$tau0,
                  mu.start = data1$mu0,trace=F))
  
  n1 <- sample(1:2,
               size = nrow(data1),
               replace = TRUE,
               prob=c(0.7, 0.3))
  
  train <- data1[n1==1,]
  test <- data1[n1==2,]
  
  mu0 <- data1$mu[n1==1]
  sigma0 <- data1$sigma[n1==1]
  nu0 <- data1$nu[n1==1]
  tau0 <- data1$tau[n1==1]
  
  fit2=try(gamlss(Surv(times,cens)~x1+x2+x3+x4+x5,family = cens(EOLLGR),
                  n.cyc=200,data=train,
                  c.crit=0.01,mu.start = train$mu0,sigma.start = train$sigma0,nu.start= train$nu0,tau.start =  train$tau0,
                  trace=F))
  
  fit.forest2 <- rfsrc(Surv(times, cens)~x1+x2+x3+x4+x5, train, ntree = 1000, 
                       mtry = 2,  
                       importance = 'permute', seed=1311) 
  
  fit.forest3 <- rfsrc(Surv(times, cens)~x1+x2+x3+x4+x5, train, ntree = 1000, 
                       mtry = 3,  
                       importance = 'permute', seed=1311)
  
  fit.forest4 <- rfsrc(Surv(times, cens)~x1+x2+x3+x4+x5, train, ntree = 1000, 
                       mtry = 4,  
                       importance = 'permute', seed=1311)
  
  if((class(fit1)[1] != "try-error" & class(fit2)[1] != "try-error")==T){
    teste <- fit1$converged
    teste2 <- fit2$converged
    
    if(teste == TRUE & teste2 == TRUE){
      betas <- c(mu=fit1$mu.coefficients ,sigma=fit1$sigma.coefficients,nu=fit1$nu.coefficients,tau=fit1$tau.coefficients)
      theta[i,] <-  betas
      res.rq[i,] = fit1$residuals
      se[[i]] <- summary(fit1,type='qr')
      se2[[i]] <- summary(fit2,type='qr')
      cens.vector[i] <- summary(as.factor(data1$cens))[1]
      
      
      temp <- concordance(fit2,newdata=train)
      c.eollgr.train[i]  <- temp$concordance
      
      
      pred <- predict(fit.forest2, train,  na.action = "na.impute")
      c.forest.train2[i] <- 1-get.cindex(time = pred$yvar[,1], 
                                         censoring = pred$yvar[,2], 
                                         predicted = pred$predicted)
      
      pred <- predict(fit.forest3, train,  na.action = "na.impute")
      c.forest.train3[i] <- 1-get.cindex(time = pred$yvar[,1], 
                                         censoring = pred$yvar[,2], 
                                         predicted = pred$predicted)
      
      pred <- predict(fit.forest4, train,  na.action = "na.impute")
      c.forest.train4[i] <- 1-get.cindex(time = pred$yvar[,1], 
                                         censoring = pred$yvar[,2], 
                                         predicted = pred$predicted)
      
      data.error.comp2[[i]] <- (gg_error(fit.forest2))
      data.error.comp3[[i]] <- (gg_error(fit.forest3))
      data.error.comp4[[i]] <- (gg_error(fit.forest4))
      
      i <- i+1
    }
    else{i <- i}
  }
  else{i <- i}
  cat("iteration = ", iter <- iter + 1, i,"\n")
}

theta450_cens30 <- theta

save(list = "theta450_cens30",
     file = "results/results_30/theta450_cens30.Rda",
     compress = "xz")

se450_cens30 <- se
p450_cens30 <- se2

save(list = "se450_cens30",
     file = "results/results_30/se450_cens30.Rda",
     compress = "xz")

save(list = "p450_cens30",
     file = "results/results_30/p450_cens30.Rda",
     compress = "xz")

cens.vector450.30 <- cens.vector

save(list = "cens.vector450.30",
     file = "results/results_30/cens.vector450.30.Rda",
     compress = "xz")

c.eollgr450.30 <- c.eollgr.train
c.forest450.30 <- data.frame(c.forest.train2,c.forest.train3,c.forest.train4)

save(list = "c.forest450.30",
     file = "results/results_30/c.forest450.30.Rda",
     compress = "xz")

save(list = "c.eollgr450.30",
     file = "results/results_30/c.eollgr450.30.Rda",
     compress = "xz")


error450_cens302 <- data.error.comp2
error450_cens303 <- data.error.comp3
error450_cens304 <- data.error.comp4

save(list = "error450_cens302",
     file = "results/results_30/error450_cens302.Rda",
     compress = "xz")

save(list = "error450_cens303",
     file = "results/results_30/error450_cens303.Rda",
     compress = "xz")

save(list = "error450_cens304",
     file = "results/results_30/error450_cens304.Rda",
     compress = "xz")
