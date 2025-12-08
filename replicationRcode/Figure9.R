###############################################################################
## Example 3: A stock with convex drift and stochastic non-convexity
## Figure 9
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
  muout <- r*s*(1-(s/K)) - q
  return(muout)
}

## profit function
wfun <- function(s,q){
  if (eta == 1){
    wout <- b*log(q)-(c/(s^gamma))*q
  } else {
    wout <- (b^(1/eta))/(1-(1/eta))*(q^(1-(1/eta)))-(c/(s^gamma))*q
  }
  return(wout)
}

###############################################################################
## Define domains and nodes

## approximation nodes: 30 nodes in (20,114.3)
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
## Possion Jump Process
###############################################################################
## the minimum stock level
slow <- 10
vlow <- vsim(cvDET,as.matrix(slow))
qlow <- qfun(slow,vlow$shadowp)
z <- c(wfun(slow,qlow)/delta)
mus <- mufun(s,qoptDET)

###############################################################################
## alpha = 0.5
alpha <- 1/2

hs <- alpha/(alpha+s)

w <- wfun(s,qoptDET) + hs*z

cvPSShalf <- vaproxpss(Aspace,s,mus,w,hs)

vPSShalf <- vsim(cvPSShalf,s1)

###############################################################################
## alpha = 2
alpha <- 2 

hs <- alpha/(alpha+s)

w <- wfun(s,qoptDET) + hs*z

cvPSS2 <- vaproxpss(Aspace,s,mus,w,hs)

vPSS2 <- vsim(cvPSS2,s1)

###############################################################################
## alpha = 5
alpha <- 5

hs <- alpha/(alpha+s)

w <- wfun(s,qoptDET) + hs*z

cvPSS5 <- vaproxpss(Aspace,s,mus,w,hs)

vPSS5 <- vsim(cvPSS5,s1)

####################################################
## Figure 9
####################################################
set1 <- data.frame(stock=vDET$stock,vfun=vDET$vfun*100,shadowp=vDET$shadowp*100,
                   key='Deterministic')
set2 <- data.frame(stock=vPSShalf$stock,vfun=vPSShalf$vfun*100,shadowp=vPSShalf$shadowp*100,
                   key='Poisson Jump 1/2')
set3 <- data.frame(stock=vPSS2$stock,vfun=vPSS2$vfun*100,shadowp=vPSS2$shadowp*100,
                   key='Poisson Jump 2')
set4 <- data.frame(stock=vPSS5$stock,vfun=vPSS5$vfun*100,shadowp=vPSS5$shadowp*100,
                   key='Poisson Jump 5')

setall <- rbind(set1,set2,set3,set4)
colnames(setall) <- c('stock','vfun','shadowp','key')
setall$key2 <- factor(setall$key, 
                      levels = c('Deterministic','Poisson Jump 1/2','Poisson Jump 2',
                                 'Poisson Jump 5'))

g1 <- ggplot(data=setall,aes(x=stock,y=vfun,col=key2,linetype=key2)) + 
  coord_cartesian(xlim=c(20,100), ylim=c(-1250,-200)) +
  geom_line(lwd=0.8) + 
  scale_linetype_manual(name='', values=c(1,2,3,4),
                        labels = c('Deterministic',
                                   expression( paste("Poisson Jump (",alpha," = ","1/2)")),
                                   expression( paste("Poisson Jump (",alpha," = ","2)")),
                                   expression( paste("Poisson Jump (",alpha," = ","5)")))) +
  scale_color_manual(name='', values=c('#000000','#3366CC','#0000FF','#9900FF'),
                     labels = c('Deterministic',
                                expression( paste("Poisson Jump (",alpha," = ","1/2)")),
                                expression( paste("Poisson Jump (",alpha," = ","2)")),
                                expression( paste("Poisson Jump (",alpha," = ","5)")))) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.box = 'vertical',
        legend.box.just = 'left',
        legend.spacing.y = unit(-0.2, "cm"),
        legend.text = element_text(size=14)) +
  labs(x='Stock',y='Value Function (V)')

g2 <- ggplot(data=setall,aes(x=stock,y=shadowp,col=key2,linetype=key2)) + 
  coord_cartesian(xlim=c(20,100), ylim=c(0,9)) +
  geom_line(lwd=0.8) + 
  scale_linetype_manual(name='', values=c(1,2,3,4,5),
                        labels = c('Deterministic',
                                   expression( paste("Poisson Jump (",alpha," = ","1/2)")),
                                   expression( paste("Poisson Jump (",alpha," = ","2)")),
                                   expression( paste("Poisson Jump (",alpha," = ","5)")),
                                   expression( paste("Poisson Jump (",alpha," = ","30)")))) +
  scale_color_manual(name='', values=c('#000000','#3366CC','#0000FF','#9900FF'),
                     labels = c('Deterministic',
                                expression( paste("Poisson Jump (",alpha," = ","1/2)")),
                                expression( paste("Poisson Jump (",alpha," = ","2)")),
                                expression( paste("Poisson Jump (",alpha," = ","5)")),
                                expression( paste("Poisson Jump (",alpha," = ","30)")))) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.box = 'vertical',
        legend.box.just = 'left',
        legend.spacing.y = unit(-0.2, "cm"),
        legend.text = element_text(size=14)) +
  labs(x='Stock',y='Accounting Price (p)')

Figure9 <- ggarrange(g1, g2, nrow=1,ncol=2, common.legend = TRUE, legend="bottom")

Figure9
