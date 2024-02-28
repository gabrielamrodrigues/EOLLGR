rm(list = ls())

source('EOLLGR_gamlss.R')

library(gamlss.cens)
library(randomForestSRC)
library(randomForest)
library(ggplot2)
library(xtable)

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
  resultado <- cbind(tamanho,parametro,valores, data.frame(media, vies,eqm))
  colnames(resultado) <- c('n','Parameters','Valor real', 'AEs','Biases', 'MSEs')
  return(resultado)
}

#Table 0% 10% and 30%

##0%

load("results/results_0/cens.vector80.0.Rda")
load("results/results_0/cens.vector150.0.Rda")
load("results/results_0/cens.vector450.0.Rda")
load("results/results_0/p80_cens0.Rda")
load("results/results_0/p150_cens0.Rda")
load("results/results_0/p450_cens0.Rda")
load("results/results_0/se80_cens0.Rda")
load("results/results_0/se150_cens0.Rda")
load("results/results_0/se450_cens0.Rda")
load("results/results_0/theta80_cens0.Rda")
load("results/results_0/theta150_cens0.Rda")
load("results/results_0/theta450_cens0.Rda")


t1_0 <- resultado_reg(theta80_cens0,80,initial)
t2_0 <- resultado_reg(theta150_cens0,150,initial)
t3_0 <- resultado_reg(theta450_cens0,450,initial)

tabela_0 <- cbind(t1_0[,1:3],rep(1,9), t1_0[,4:6],rep(1,9),t2_0[,4:6],rep(1,9),t3_0[,4:6])

xtable(tabela_0,digits = 2)


#10%

load("results/results_10/cens.vector80.10.Rda")
load("results/results_10/cens.vector150.10.Rda")
load("results/results_10/cens.vector450.10.Rda")
load("results/results_10/p80_cens10.Rda")
load("results/results_10/p150_cens10.Rda")
load("results/results_10/p450_cens10.Rda")
load("results/results_10/se80_cens10.Rda")
load("results/results_10/se150_cens10.Rda")
load("results/results_10/se450_cens10.Rda")
load("results/results_10/theta80_cens10.Rda")
load("results/results_10/theta150_cens10.Rda")
load("results/results_10/theta450_cens10.Rda")




t1_10 <- resultado_reg(theta80_cens10,80,initial)
t2_10 <- resultado_reg(theta150_cens10,150,initial)
t3_10 <- resultado_reg(theta450_cens10,450,initial)

tabela_10 <- cbind(t1_10[,1:3],rep(1,9), t1_10[,4:6],rep(1,9),t2_10[,4:6],rep(1,9),t3_10[,4:6])

xtable(tabela_10,digits = 2)


##30%

load("results/results_30/cens.vector80.30.Rda")
load("results/results_30/cens.vector150.30.Rda")
load("results/results_30/cens.vector450.30.Rda")
load("results/results_30/se80_cens30.Rda")
load("results/results_30/se150_cens30.Rda")
load("results/results_30/se450_cens30.Rda")
load("results/results_30/theta80_cens30.Rda")
load("results/results_30/theta150_cens30.Rda")
load("results/results_30/theta450_cens30.Rda")



t1_30 <- resultado_reg(theta80_cens30,80,initial)
t2_30 <- resultado_reg(theta150_cens30,150,initial)
t3_30 <- resultado_reg(theta450_cens30,450,initial)

tabela_30 <- cbind(t1_30[,1:3],rep(1,9), t1_30[,4:6],rep(1,9),t2_30[,4:6],rep(1,9),t3_30[,4:6])

xtable(tabela_30,digits = 2)


## cindex 0% 10% and 30%


load("results/results_0/c.eollgr80.0.Rda")
load("results/results_0/c.eollgr150.0.Rda")
load("results/results_0/c.eollgr450.0.Rda")
load("results/results_0/c.forest80.0.Rda")
load("results/results_0/c.forest150.0.Rda")
load("results/results_0/c.forest450.0.Rda")

load("results/results_10/c.eollgr80.10.Rda")
load("results/results_10/c.eollgr150.10.Rda")
load("results/results_10/c.eollgr450.10.Rda")
load("results/results_10/c.forest80.10.Rda")
load("results/results_10/c.forest150.10.Rda")
load("results/results_10/c.forest450.10.Rda")

load("results/results_30/c.eollgr80.30.Rda")
load("results/results_30/c.eollgr150.30.Rda")
load("results/results_30/c.eollgr450.30.Rda")
load("results/results_30/c.forest80.30.Rda")
load("results/results_30/c.forest150.30.Rda")
load("results/results_30/c.forest450.30.Rda")




model=rep('eollgr',3000)
tab0.eollgr <- data.frame(n=as.factor(c(c(rep(80,1000),rep(150,1000),rep(450,1000)))),
                          model=rep('EOLLGR',3000),
                          cindex=c(c.eollgr80.0,c.eollgr150.0,c.eollgr450.0))   


tab0.forest <- data.frame(n=as.factor(c(c(rep(80,1000),rep(150,1000),rep(450,1000)))),
                          model=rep('Forest',3000),
                          cindex=c(c.forest80.0[,2],c.forest150.0[,2],c.forest450.0[,2]))

tab0 <- rbind(tab0.eollgr, tab0.forest)
summary(tab0$n)

p0 <- ggplot(tab0, aes(x=as.factor(n) ,y=cindex, colour=as.factor(model))) +
  geom_boxplot(outlier.size=2)+
  labs(title="",x="Sample size", y = "C-Index")

p0+theme(legend.background = element_rect(fill = "transparent"),
         axis.title = element_text( size = (12)),
         axis.text = element_text(size=11),
         legend.position = c(0.83,0.15)) +
  theme(legend.key.width = unit(1.2,"cm"))+
  labs(colour = "Model")+ylim(0.45,1)


tab10.eollgr <- data.frame(n=as.factor(c(c(rep(80,1000),rep(150,1000),rep(450,1000)))),
                           model=rep('EOLLGR',3000),
                           cindex=c(c.eollgr80.10,c.eollgr150.10,c.eollgr450.10))   


tab10.forest <- data.frame(n=as.factor(c(c(rep(80,1000),rep(150,1000),rep(450,1000)))),
                           model=rep('Forest',3000),
                           cindex=c(c.forest80.10[,2],c.forest150.10[,2],c.forest450.10[,2]))

tab10 <- rbind(tab10.eollgr, tab10.forest)


p10 <- ggplot(tab10, aes(x=as.factor(n) ,y=cindex, colour=as.factor(model))) +
  geom_boxplot(outlier.size=2)+
  labs(title="",x="Tamanho amostral", y = "C-Index")

p10+theme(legend.background = element_rect(fill = "transparent"),
          axis.title = element_text( size = (12)),
          axis.text = element_text(size=11),
          legend.position = c(0.83,0.15)) +
  theme(legend.key.width = unit(1.2,"cm"))+
  labs(colour = "Modelo")+ylim(0.45,1)


tab30.eollgr <- data.frame(n=as.factor(c(c(rep(80,1000),rep(150,1000),rep(450,1000)))),
                           model=rep('EOLLGR',3000),
                           cindex=c(c.eollgr80.30,c.eollgr150.30,c.eollgr450.30))   


tab30.forest <- data.frame(n=as.factor(c(c(rep(80,1000),rep(150,1000),rep(450,1000)))),
                           model=rep('Forest',3000),
                           cindex=c(c.forest80.30[,3],c.forest150.30[,3],c.forest450.30[,3]))

tab30 <- rbind(tab30.eollgr, tab30.forest)


p30 <- ggplot(tab30, aes(x=as.factor(n) ,y=cindex, colour=as.factor(model))) +
  geom_boxplot(outlier.size=2)+
  labs(title="",x="Tamanho amostral", y = "C-Index")

p30+theme(legend.background = element_rect(fill = "transparent"),
          axis.title = element_text( size = (12)),
          axis.text = element_text(size=11),
          legend.position = c(0.83,0.15)) +
  theme(legend.key.width = unit(1.2,"cm"))+
  labs(colour = "Modelo")+ylim(0.45,1)


## Results 70%

## Table


##70%

load("results/results_70/cens.vector1107.70.Rda")
load("results/results_70/se1107_cens70.Rda")
load("results/results_70/theta1107_cens70.Rda")


t1_70 <- resultado_reg(theta1107_cens70,1107,initial)
tabela_70 <- t1_70
xtable(tabela_70,digits = 2)


load("results/results_70/c.eollgr1107.70.Rda")
load("results/results_70/c.forest1107.70.Rda")


data.frame(rep(1107,1000),c.eollgr1107.70)



model=rep('eollgr',1000)
tab0.eollgr <- data.frame(n=as.factor(c(c(rep(1107,1000)))),
                          model=rep('EOLLGR',1000),
                          cindex=c(c.eollgr1107.70))   


tab0.forest <- data.frame(n=as.factor(c(c(rep(1107,1000)))),
                          model=rep('Forest',1000),
                          cindex=c(c.forest1107.70))

tab0 <- rbind(tab0.eollgr, tab0.forest)

p0 <- ggplot(tab0, aes(x=as.factor(n) ,y=cindex, colour=as.factor(model))) +
  geom_boxplot(outlier.size=2)+
  labs(title="",x="Sample size", y = "C-Index")

p0+theme(legend.background = element_rect(fill = "transparent"),
         axis.title = element_text( size = (12)),
         axis.text = element_text(size=11),
         legend.position = c(0.83,0.15)) +
  theme(legend.key.width = unit(1.2,"cm"))+
  labs(colour = "Model")+ylim(0.45,1)
