###############################################################################
## Example 1: optimized and non-optimized renewable resource management
##            under convexity
## 3.1.2. Non-GBM stochasticity
## Figure 5 and Figure 6
###############################################################################
## R-packages used in the examples
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(rootSolve)
library(tidyverse)

rm(list=ls())
## load pre-defined functions
source(url('https://raw.githubusercontent.com/ysd2004/stochasticcapnJAERE/refs/heads/main/data_and_fncs/functions.R'))

## parameters (see Table C.1)
b <- 1
eta <- 1/2
c <- 5
gamma <- 2
delta <- 0.05


R <- 1
M <- 0.53
C <- 0.03
r <- R*(M-C)
K1 <- 114.3
L <- 1/K1
K <- 100
G <- ((M-C)/K)-M*L
mu <- 1

## net growth function
gfun <- function(s){
  netgrowth <- r*s*(1-(s/K))
  return(netgrowth)
}

## catch function
qfun <- function(s,Vs){
  qout <- b*((Vs+c/(s^gamma))^(-eta))
  return(qout)
}

## drift function
mufun <- function(s,q){
  muout <- growthfun(s) - q
  return(muout)
}


###############################################################################
## Define domains and nodes

## approximation nodes: 30 nodes in (202,114.3)
nnode <- 30
lower <- 20
upper <- K1
s <- chebnodegen(nnode,lower,upper)
Aspace <- aproxdef(nnode,lower,upper,delta)

## prediction nodes: nnode*10+1 following Sims, Horan, and Meadows (2018)
s1 <- seq(lower,upper,length.out=(nnode*10+1))
s1 <- s1[s1<K]

###############################################################################
## Deterministic
###############################################################################
cvDET <- vaproxsc(Aspace,s,gfun)
vDET <- vsim(cvDET,s1)

qoptDET <- qfun(s,vsim(cvDET,s)$shadowp)
###############################################################################

###############################################################################
## Figure 5 in the manuscript 
## GBM with theta0 = 0.1
###############################################################################
theta0 <- 0.1
sigsGBM <- as.matrix((theta0*s)^2,col=1)

cvGBM <- vaproxsc(Aspace,s,gfun,sigsGBM)
vGBM <- vsim(cvGBM,s1)

## GBM steady-state
ssGBM <- 57.94068

## DEM+ENV
thetaM <- 0.1481263

sigsDEM <- as.matrix(r*s*(((M+C)/(M-C))+(1-2*mu)*(s/K))+((thetaM*R*s*(1-s/K1))^2),col=1)

cvDEM <- vaproxsc(Aspace,s,gfun,sigsDEM)
vDEM <- vsim(cvDEM,s1)

## Figure data setup
set1 <- data.frame(stock=vDET$stock,vfun=vDET$vfun*100,shadowp=vDET$shadowp*100,
                   key='Deterministic')
set2 <- data.frame(stock=vGBM$stock,vfun=vGBM$vfun*100,shadowp=vGBM$shadowp*100,
                   key='GBM')
set3 <- data.frame(stock=vDEM$stock,vfun=vDEM$vfun*100,shadowp=vDEM$shadowp*100,
                   key='DEM+ENV')
setall <- rbind(set1,set2,set3)
colnames(setall) <- c('stock','vfun','shadowp','key')
setall$key2 <- factor(setall$key, 
                      levels = c('Deterministic','GBM','DEM+ENV'))

g1 <- ggplot(data=setall,aes(x=stock,y=vfun,col=key2,linetype=key2)) + 
  coord_cartesian(xlim=c(20,100), ylim=c(-300,-190)) +
  geom_line(lwd=0.8) + 
  scale_linetype_manual(name='', values=c(1,2,3),
                        labels = c('Deterministic',
                                   expression( paste("GBM (",theta[0]," = ","0.1)")),
                                   expression( paste("DEM+ENV (",theta[M]," = ","0.1481 at ssGBM)")))) +
  scale_color_manual(name='', values=c('#000000','#3366CC','#9900FF'),
                     labels = c('Deterministic',
                                expression( paste("GBM (",theta[0]," = ","0.1)")),
                                expression( paste("DEM+ENV (",theta[M]," = ","0.1481 at ssGBM)")))) +
  geom_vline(xintercept=ssGBM,linetype=2,colour='#3366CC') +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.box = 'vertical',
        legend.box.just = 'left',
        legend.spacing.y = unit(-0.2, "cm"),
        legend.text = element_text(size=14)) +
  annotate(geom="text", x=ssGBM, y=-300, label="GBM Steady State (ssGBM = 57.9407)",color="black") +
  labs(x='Stock (% K)',y='Value Function (V)')

g2 <- ggplot(data=setall,aes(x=stock,y=shadowp,col=key2,linetype=key2)) + 
  coord_cartesian(xlim=c(20,100), ylim=c(0,6.5)) +
  geom_line(lwd=0.8) + 
  scale_linetype_manual(name='', values=c(1,2,3),
                        labels = c('Deterministic',
                                   expression( paste("GBM (",theta[0]," = ","0.1)")),
                                   expression( paste("DEM+ENV (",theta[M]," = ","0.1481 at ssGBM)")))) +
  scale_color_manual(name='', values=c('#000000','#3366CC','#9900FF'),
                     labels = c('Deterministic',
                                expression( paste("GBM (",theta[0]," = ","0.1)")),
                                expression( paste("DEM+ENV (",theta[M]," = ","0.1481 at ssGBM)")))) +
  geom_vline(xintercept=ssGBM,linetype=2,colour='#3366CC') +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.box = 'vertical',
        legend.box.just = 'left',
        legend.spacing.y = unit(-0.2, "cm"),
        legend.text = element_text(size=14)) +
  annotate(geom="text", x=ssGBM, y=6.5, label="GBM Steady State (ssGBM = 57.9407)",color="black") +
  labs(x='Stock (% K)',y='Accounting Price (p)')

Figure5 <- ggarrange(g1, g2, nrow=1,ncol=2, common.legend = TRUE, legend="bottom")

Figure5
###############################################################################

###############################################################################
## Figure 6 in the manuscript 
## GBM with theta0 = 0.15
###############################################################################
theta0 <- 0.15
sigsGBM <- as.matrix((theta0*s)^2,col=1)

cvGBM <- vaproxsc(Aspace,s,gfun,sigsGBM)
vGBM <- vsim(cvGBM,s1)

## GBM steady-state
ssGBM <- 58.57256

## DEM+ENV
thetaM <- 0.2747174

sigsDEM <- as.matrix(r*s*(((M+C)/(M-C))+(1-2*mu)*(s/K))+((thetaM*R*s*(1-s/K1))^2),col=1)

cvDEM <- vaproxsc(Aspace,s,gfun,sigsDEM)
vDEM <- vsim(cvDEM,s1)

## Figure data setup
set1 <- data.frame(stock=vDET$stock,vfun=vDET$vfun*100,shadowp=vDET$shadowp*100,
                   key='Deterministic')
set2 <- data.frame(stock=vGBM$stock,vfun=vGBM$vfun*100,shadowp=vGBM$shadowp*100,
                   key='GBM')
set3 <- data.frame(stock=vDEM$stock,vfun=vDEM$vfun*100,shadowp=vDEM$shadowp*100,
                   key='DEM+ENV')
setall <- rbind(set1,set2,set3)
colnames(setall) <- c('stock','vfun','shadowp','key')
setall$key2 <- factor(setall$key, 
                      levels = c('Deterministic','GBM','DEM+ENV'))

g1 <- ggplot(data=setall,aes(x=stock,y=vfun,col=key2,linetype=key2)) + 
  coord_cartesian(xlim=c(20,100), ylim=c(-300,-190)) +
  geom_line(lwd=0.8) + 
  scale_linetype_manual(name='', values=c(1,2,3),
                        labels = c('Deterministic',
                                   expression( paste("GBM (",theta[0]," = ","0.15)")),
                                   expression( paste("DEM+ENV (",theta[M]," = ","0.2747 at ssGBM)")))) +
  scale_color_manual(name='', values=c('#000000','#3366CC','#9900FF'),
                     labels = c('Deterministic',
                                expression( paste("GBM (",theta[0]," = ","0.15)")),
                                expression( paste("DEM+ENV (",theta[M]," = ","0.2747 at ssGBM)")))) +
  geom_vline(xintercept=ssGBM,linetype=2,colour='#3366CC') +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.box = 'vertical',
        legend.box.just = 'left',
        legend.spacing.y = unit(-0.2, "cm"),
        legend.text = element_text(size=14)) +
  annotate(geom="text", x=ssGBM, y=-300, label="GBM Steady State (ssGBM = 58.5726)",color="black") +
  labs(x='Stock (% K)',y='Value Function (V)')

g2 <- ggplot(data=setall,aes(x=stock,y=shadowp,col=key2,linetype=key2)) + 
  coord_cartesian(xlim=c(20,100), ylim=c(0,6.5)) +
  geom_line(lwd=0.8) + 
  scale_linetype_manual(name='', values=c(1,2,3),
                        labels = c('Deterministic',
                                   expression( paste("GBM (",theta[0]," = ","0.15)")),
                                   expression( paste("DEM+ENV (",theta[M]," = ","0.2747 at ssGBM)")))) +
  scale_color_manual(name='', values=c('#000000','#3366CC','#9900FF'),
                     labels = c('Deterministic',
                                expression( paste("GBM (",theta[0]," = ","0.15)")),
                                expression( paste("DEM+ENV (",theta[M]," = ","0.2747 at ssGBM)")))) +
  geom_vline(xintercept=ssGBM,linetype=2,colour='#3366CC') +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.box = 'vertical',
        legend.box.just = 'left',
        legend.spacing.y = unit(-0.2, "cm"),
        legend.text = element_text(size=14)) +
  annotate(geom="text", x=ssGBM, y=6.5, label="GBM Steady State (ssGBM = 58.5726)",color="black") +
  labs(x='Stock (% K)',y='Accounting Price (p)')

Figure6 <- ggarrange(g1, g2, nrow=1,ncol=2, common.legend = TRUE, legend="bottom")

Figure6
###############################################################################

