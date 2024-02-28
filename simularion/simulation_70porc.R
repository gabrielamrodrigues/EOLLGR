
rm(list = ls())
source('~/Esalq doutorado/Artigos/regressao_distribuicao/EOLLGR/EOLLGR_gamlss.R')

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
  # x6 <- rnorm(n,0,1)
  # x7 <- rnorm(n,0,0.5)
  
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

# Calculate results
resultado_reg <- function(estimativas,n,valores){
  media <- apply(estimativas, 2, mean) 
  sd <- apply(estimativas, 2, sd)
  errop <- sd/sqrt(n) 
  var <- apply(estimativas, 2,var) 
  eqm <- c( var[1]+(media[1]-valores[1])^2,
            var[2]+(media[2]-valores[2])^2,
            var[3]+(media[3]-valores[3])^2,
            var[4]+(media[4]-valores[4])^2,
            var[5]+(media[5]-valores[5])^2,
            var[6]+(media[6]-valores[6])^2,
            var[7]+(media[7]-valores[7])^2,
            var[8]+(media[8]-valores[8])^2,
            var[9]+(media[9]-valores[9])^2)
  vies <- c( (media[1]-valores[1]),
             (media[2]-valores[2]),
             (media[3]-valores[3]),
             (media[4]-valores[4]),
             (media[5]-valores[5]),
             (media[6]-valores[6]),
             (media[7]-valores[7]),
             (media[8]-valores[8]),
             (media[9]-valores[9]))
  parametro <- rep(c('beta10', 'beta11','beta12', 'beta13','beta14','beta15','beta20', 'beta30','beta40'))
  tamanho <- rep(n,9)
  resultado <- cbind(tamanho,valores,parametro, data.frame(media, vies,eqm))
  colnames(resultado) <- c('n','Valor real','Parameters', 'AEs','Biases', 'MSEs')
  return(resultado)
}


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

#70% n=80

porc.cens <- 70
tau <- 0.4
npar <- 9
n <- 1107
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
  
  fit.forest2 <- rfsrc(Surv(times, cens)~x1+x2+x3+x4+x5, train, ntree =
                         1000, nodesize = 10, mtry = 2, nsplit = 10, importance = "permute", splitrule = "logrank", seed=1311) 
  
  
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
      
      
      
      
      i <- i+1
    }
    else{i <- i}
  }
  else{i <- i}
  cat("iteration = ", iter <- iter + 1, i,"\n")
}


theta1107_cens70 <- theta

save(list = "theta1107_cens70",
     file = "results/results_70/theta1107_cens70.Rda",
     compress = "xz")

se1107_cens70 <- se
p1107_cens70 <- se2

save(list = "se1107_cens70",
     file = "results/results_70/se1107_cens70.Rda",
     compress = "xz")

save(list = "p1107_cens70",
     file = "results/results_70/p1107_cens70.Rda",
     compress = "xz")

cens.vector1107.70 <- cens.vector

save(list = "cens.vector1107.70",
     file = "results/results_70/cens.vector1107.70.Rda",
     compress = "xz")

c.eollgr1107.70 <- c.eollgr.train
c.forest1107.70 <- c.forest.train2

save(list = "c.forest1107.70",
     file = "results/results_70/c.forest1107.70.Rda",
     compress = "xz")

save(list = "c.eollgr1107.70",
     file = "results/results_70/c.eollgr1107.70.Rda",
     compress = "xz")


error1107_cens702 <- data.error.comp2
error1107_cens703 <- data.error.comp3
error1107_cens704 <- data.error.comp4

save(list = "error1107_cens702",
     file = "results/results_70/error1107_cens702.Rda",
     compress = "xz")

save(list = "error1107_cens703",
     file = "results/results_70/error1107_cens703.Rda",
     compress = "xz")

save(list = "error1107_cens704",
     file = "results/results_70/error1107_cens704.Rda",
     compress = "xz")

