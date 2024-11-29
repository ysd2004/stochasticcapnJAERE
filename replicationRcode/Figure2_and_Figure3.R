###############################################################################
## Example 1: optimized and non-optimized renewable resource management
##            under convexity
## Figure 2 and Figure 3
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

###############################################################################
## Define parameters
r <- 0.5      # instrinsic growth rate
K <- 1        # carrying capacity
b <- 1
eta <- 0.5
c <- 5
gamma <- 2
delta <- 0.05

###############################################################################
## Steady state: deterministic
sig <- 0
phi <- (2*(b^2)+2*b*(((b^2)+c*((r+delta-(sig^2))^2)))^(1/2))/(((r+delta-(sig^2))^2))

ssfun <- function(s){
  ss <- r*s*(1-(s/K)) - b*((phi +c)^(-1/2))*s
  return(ss)
}

ssd <- uniroot(ssfun,c(0.2,0.6))$root

## Steady state: stochastic when sig = 0.1
sig <- 0.1
phi <- (2*(b^2)+2*b*(((b^2)+c*((r+delta-(sig^2))^2)))^(1/2))/(((r+delta-(sig^2))^2))

ssfun <- function(s){
  ss <- r*s*(1-(s/K)) - b*((phi +c)^(-1/2))*s
  return(ss)
}

ss <- uniroot(ssfun,c(0.3,0.8))$root

#########################################################################
## Define domains and nodes

## approximation nodes: 35 nodes in (0.2,1.4)
nnode <- 35
lower <- 0.2
upper <- 1.4
stockapx <- chebnodegen(35,0.2,1.4)

## simulation: 35 nodes in (0.2,1.1)
stock0 <- chebnodegen(35,0.2,1.1)

## Fig 1: 50 nodes in (0.01,1)
stockfig <- chebnodegen(50,0.01,1)

## Optimal
#########################################################################
sig <- 0.1
catch <- as.matrix(b*((phi +c)^(-1/2))*stockapx,ncol=1)
profit <- as.matrix(-((b^2)/catch)-(c/(stockapx^2))*catch,ncol=1)
mus <- as.matrix(r*stockapx*(1-(stockapx/K)) - catch,ncol=1)
sigs <- as.matrix((sig*stockapx)^2,ncol=1)

catch.opt <- catch
catch.opt.fig <- as.matrix(b*((phi +c)^(-1/2))*stockfig,ncol=1)

## stochastic sig=0.1
Aspace <- aproxdef(nnode,lower,upper,delta)
vC.opt.s <- vaprox(Aspace,stockapx,mus,profit,sigs)
opt.s <- vsim(vC.opt.s,as.matrix(stock0,ncol=1),profit)

## deterministic sig=0.1
vC.opt.d <- vaprox(Aspace,stockapx,mus,profit)
opt.d <- vsim(vC.opt.d,as.matrix(stock0,ncol=1),profit)


### 1/2 catch
#########################################################################
catch <- as.matrix(0.5*b*((phi +c)^(-1/2))*stockapx,ncol=1)
profit <- as.matrix(-((b^2)/catch)-(c/(stockapx^2))*catch,ncol=1)
mus <- as.matrix(r*stockapx*(1-(stockapx/K)) - catch,ncol=1)
sigs <- as.matrix((sig*stockapx)^2,ncol=1)

catch.half <- catch
catch.half.fig <- as.matrix(0.5*b*((phi +c)^(-1/2))*stockfig,ncol=1)

## stochastic sig=0.1
vC.opt.half.s <- vaprox(Aspace,stockapx,mus,profit,sigs)
opt.half.s <- vsim(vC.opt.half.s,as.matrix(stock0,ncol=1),profit)

## deterministic sig=0.1
vC.opt.half.d <- vaprox(Aspace,stockapx,mus,profit)
opt.half.d <- vsim(vC.opt.half.d,as.matrix(stock0,ncol=1),profit)


### 1.5 catch
#########################################################################
catch <- as.matrix(1.5*b*((phi +c)^(-1/2))*stockapx,ncol=1)
profit <- as.matrix(-((b^2)/catch)-(c/(stockapx^2))*catch,ncol=1)
mus <- as.matrix(r*stockapx*(1-(stockapx/K)) - catch,ncol=1)
sigs <- as.matrix((sig*stockapx)^2,ncol=1)

catch.2x <- catch
catch.2x.fig <- as.matrix(1.5*b*((phi +c)^(-1/2))*stockfig,ncol=1)

## stochastic sig=0.1
vC.opt.2x.s <- vaprox(Aspace,stockapx,mus,profit,sigs)
opt.2x.s <- vsim(vC.opt.2x.s,as.matrix(stock0,ncol=1),profit)

## deterministic sig=0.1
vC.opt.2x.d <- vaprox(Aspace,stockapx,mus,profit)
opt.2x.d <- vsim(vC.opt.2x.d,as.matrix(stock0,ncol=1),profit)

### precautious: adaptive
#########################################################################
adjust <- 1/(1+1.1*(ss-stockapx))

catch <- as.matrix(adjust*b*((phi +c)^(-1/2))*stockapx,ncol=1)
profit <- as.matrix(-((b^2)/catch)-(c/(stockapx^2))*catch,ncol=1)
mus <- as.matrix(r*stockapx*(1-(stockapx/K)) - catch,ncol=1)
sigs <- as.matrix((sig*stockapx)^2,ncol=1)

catch.adp <- catch

adjust.fig <- 1/(1+1.1*(ss-stockfig))

catch.fig <- as.matrix(adjust.fig*b*((phi +c)^(-1/2))*stockfig,ncol=1)
catch.adp.fig <- as.matrix(adjust.fig*b*((phi +c)^(-1/2))*stockfig,ncol=1)

## stochastic sig=0.1
vC.adaptive.s <- vaprox(Aspace,stockapx,mus,profit,sigs)
adaptive.s <- vsim(vC.adaptive.s,as.matrix(stock0,ncol=1),profit)

## deterministic sig=0.1
vC.adaptive.d <- vaprox(Aspace,stockapx,mus,profit)
adaptive.d <- vsim(vC.adaptive.d,as.matrix(stock0,ncol=1),profit)

## Optimal: deterministic sig=0
vC.opt.d0 <- vaprox(Aspace,stockapx,mus,profit)
opt.d0 <- vsim(vC.opt.d0,as.matrix(stock0,ncol=1),profit)

## (1/2) catch
catch <- as.matrix(0.5*b*((phi +c)^(-1/2))*stockapx,ncol=1)
profit <- as.matrix(-((b^2)/catch)-(c/(stockapx^2))*catch,ncol=1)
mus <- as.matrix(r*stockapx*(1-(stockapx/K)) - catch,ncol=1)

catch.half.d0 <- catch

vC.half.d0 <- vaprox(Aspace,stockapx,mus,profit)
opt.half.d0 <- vsim(vC.half.d0,as.matrix(stock0,ncol=1),profit)

## (1.5) catch
catch <- as.matrix(1.5*b*((phi +c)^(-1/2))*stockapx,ncol=1)
profit <- as.matrix(-((b^2)/catch)-(c/(stockapx^2))*catch,ncol=1)
mus <- as.matrix(r*stockapx*(1-(stockapx/K)) - catch,ncol=1)

catch.2x.d0 <- catch

vC.2x.d0 <- vaprox(Aspace,stockapx,mus,profit)
opt.2x.d0 <- vsim(vC.2x.d0,as.matrix(stock0,ncol=1),profit)

## Adaptive
adjust <- 1/(1+1.1*(ss-stockapx))

catch <- as.matrix(adjust*b*((phi +c)^(-1/2))*stockapx,ncol=1)
profit <- as.matrix(-((b^2)/catch)-(c/(stockapx^2))*catch,ncol=1)
mus <- as.matrix(r*stockapx*(1-(stockapx/K)) - catch,ncol=1)
sigs <- as.matrix((sig*stockapx)^2,ncol=1)

catch.adp0 <- catch

## deterministic
vC.adaptive.d0 <- vaprox(Aspace,stockapx,mus,profit)
adaptive.d0 <- vsim(vC.adaptive.d0,as.matrix(stock0,ncol=1),profit)

#########################################################################
## Figure 2 in the manuscript
set1 <- data.frame(stock=opt.s$stock,vfun=opt.s$vfun,shadowp=opt.s$shadowp,
                   key1='Optimal',key2='Stochastic')
set2 <- data.frame(stock=opt.d0$stock,vfun=opt.d0$vfun,shadowp=opt.d0$shadowp,
                   key1='Optimal',key2='Deterministic')

set3 <- data.frame(stock=opt.half.s$stock,vfun=opt.half.s$vfun,shadowp=opt.half.s$shadowp,
                   key1='(1/2) x Optimal',key2='Stochastic')
set4 <- data.frame(stock=opt.half.d$stock,vfun=opt.half.d$vfun,shadowp=opt.half.d0$shadowp,
                   key1='(1/2) x Optimal',key2='Deterministic')
set5 <- data.frame(stock=opt.2x.s$stock,vfun=opt.2x.s$vfun,shadowp=opt.2x.s$shadowp,
                   key1='(1.5) x Optimal',key2='Stochastic')
set6 <- data.frame(stock=opt.2x.d$stock,vfun=opt.2x.d$vfun,shadowp=opt.2x.d0$shadowp,
                   key1='(1.5) x Optimal',key2='Deterministic')

setall <- rbind(set1,set2,set3,set4,set5,set6)

colnames(setall) <- c('stock','vfun','shadowp','key1','key2')
setall$key1 <- factor(setall$key1, levels = c('Optimal','(1/2) x Optimal','(1.5) x Optimal'))
setall$key2 <- factor(setall$key2, levels = c('Deterministic',
                                              'Stochastic'))

g1 <- ggplot(data=setall,aes(x=stock,y=vfun,col=key1,linetype=key2)) + 
  coord_cartesian(xlim=c(0.2,1.1), ylim=c(-350,-180)) +
  geom_line(lwd=0.8) + 
  scale_linetype_manual(name='Stochasticity Level', values=c(1,2)) +
  scale_color_manual(name='Economic Program', values=c("#000000","#3366CC","#FF6633"),
                     labels = c(expression( paste( "Optimal (",sigma," varying)")),
                                expression( paste( "(0.5) x Optimal (",sigma," = ","0.1)" )),
                                expression( paste("1.5 x Optimal (",sigma," = ","0.1)" )))) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.box = 'vertical',
        legend.box.just = 'left',
        legend.spacing.y = unit(-0.2, "cm"),
        legend.text = element_text(size=10)) +
  labs(x='Stock',y='Intertemporal Welfare (V)')

g2 <- ggplot(data=setall,aes(x=stock,y=shadowp,col=key1,linetype=key2)) + 
  coord_cartesian(xlim=c(0.2,1.1), ylim=c(0,600)) +
  geom_line(lwd=0.8) + 
  scale_linetype_manual(name='Stochasticity Level', values=c(1,2)) +
  scale_color_manual(name='Economic Program', values=c("#000000","#3366CC","#FF6633"),
                     labels = c(expression( paste( "Optimal (",sigma," varying)")),
                                expression( paste( "(0.5) x Optimal (",sigma," = ","0.1)" )),
                                expression( paste("1.5 x Optimal (",sigma," = ","0.1)" )))) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.box = 'vertical',
        legend.box.just = 'left',
        legend.spacing.y = unit(-0.2, "cm"),
        legend.text = element_text(size=10)) +
  labs(x='Stock',y='Accounting Price (p)')

Figure2 <- ggarrange(g1, g2, nrow=1,ncol=2, common.legend = TRUE, legend="bottom")

Figure2

#########################################################################
## Figure 3 in the manuscript
set1 <- data.frame(stock=opt.s$stock,vfun=opt.s$vfun,shadowp=opt.s$shadowp,
                   key1='Optimal',key2='Stochastic')
set2 <- data.frame(stock=opt.d$stock,vfun=opt.d$vfun,shadowp=opt.d$shadowp,
                   key1='Optimal',key2='Deterministic')
set3 <- data.frame(stock=adaptive.s$stock,vfun=adaptive.s$vfun,shadowp=adaptive.s$shadowp,
                   key1='Adaptive',key2='Stochastic')
set4 <- data.frame(stock=adaptive.d$stock,vfun=adaptive.d$vfun,shadowp=adaptive.d$shadowp,
                   key1='Adaptive',key2='Deterministic')

setall <- rbind(set1,set2,set3,set4)

colnames(setall) <- c('stock','vfun','shadowp','key1','key2')
setall$key1 <- factor(setall$key1, levels = c('Optimal','Adaptive'))
setall$key2 <- factor(setall$key2, levels = c('Stochastic',
                                              'Deterministic'))

g1 <- ggplot(data=setall,aes(x=stock,y=vfun,col=key1,linetype=key2)) + 
  coord_cartesian(xlim=c(0.2,1.1), ylim=c(-280,-180)) +
  geom_line(lwd=0.8) + 
  scale_linetype_manual(name='Stochasticity Level', values=c(2,1)) +
  scale_color_manual(name='Economic Program', values=c('#000000','#CC0000'),
                     labels = c('Optimal','Adaptive')) +
  geom_vline(xintercept=ss,linetype=2) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.box = 'vertical',
        legend.box.just = 'left',
        legend.spacing.y = unit(-0.2, "cm"),
        legend.text = element_text(size=15)) +
  labs(x='Stock',y='Intertemporal Welfare (V)') +
  annotate(geom="text", x=ss, y=-279, label="Steady State",color="black") 

g2 <- ggplot(data=setall,aes(x=stock,y=shadowp,col=key1,linetype=key2)) + 
  coord_cartesian(xlim=c(0.2,1.1), ylim=c(0,500)) +
  geom_line(lwd=0.8) + 
  scale_linetype_manual(name='Stochasticity Level', values=c(2,1)) +
  scale_color_manual(name='Economic Program', values=c('#000000','#CC0000'),
                     labels = c('Optimal','Adaptive')) +
  geom_vline(xintercept=ss,linetype=2) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.box = 'vertical',
        legend.box.just = 'left',
        legend.spacing.y = unit(-0.2, "cm"),
        legend.text = element_text(size=15)) +
  labs(x='Stock',y='Accounting Price (p)') +
  annotate(geom="text", x=ss, y=10, label="Steady State",color="black") 

Figure3 <- ggarrange(g1, g2, nrow=1,ncol=2, common.legend = TRUE, legend="bottom")

Figure3