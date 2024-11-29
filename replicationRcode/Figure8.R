###############################################################################
## Example 2: A stock with non-convex drift and smooth stochasticity
## Figure 8
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

#########################################################################
## load simulation data
load(url('https://raw.githubusercontent.com/ysd2004/stochasticcapnJAERE/refs/heads/main/data_and_fncs/simresult.rda'))

#########################################################################
## Getting accounting price curves with spline and Chebyshev polynomials

## deterministic case
determp1 <- splinefun(sto00$stock[1:160],sto00$vfun[1:160])
determp2 <- splinefun(sto00$stock[161:600],sto00$vfun[161:600])
x <- sto00$stock[1:160]
nodes <- chebnodegen(100,x[1]*1.05,x[160]*0.98)
Aspace <- aproxdef(40,x[1]*1.05,x[160]*0.98,0.02)
vc0p1 <- chebc(Aspace,nodes,determp1(nodes))
nodesp1 <- chebnodegen(100,x[1]*1.1,x[160]*0.98)
vsim0p1 <- vsim(vc0p1,nodesp1)
x2 <- sto00$stock[161:600]
nodes2 <- chebnodegen(100,x2[1]*1.005,x2[440]*0.95)
Aspace2 <- aproxdef(50,x2[1]*1.005,x2[440]*0.95,0.02)
vc0p2 <- chebc(Aspace2,nodes2,determp2(nodes2))
nodesp2 <- chebnodegen(100,x2[1],x2[440]*0.93)
vsim0p2 <- vsim(vc0p2,nodes2)

## stochastic cases
sig01 <- splinefun(sto01$stock,sto01$vfun)
sig02 <- splinefun(sto02$stock,sto02$vfun)
sig03 <- splinefun(sto03$stock,sto03$vfun)

## Chevyshev nodes
x <- sto01$stock
nodes <- chebnodegen(300,x[1]*1.05,x[600]*0.99)
nodes2 <- nodes[5:277]
nodes3 <- c(vsim0p1$stock,vsim0p2$stock)

Aspace <- aproxdef(62,x[1]*1.05,x[600]*0.99,0.02)
vc1 <- chebc(Aspace,nodes,sig01(nodes))
vsim1 <- vsim(vc1,nodes2)
vsim12 <- vsim(vc1,nodes3)

Aspace <- aproxdef(27,x[1]*1.05,x[600]*0.99,0.02)
vc2 <- chebc(Aspace,nodes,sig02(nodes))
vsim2 <- vsim(vc2,nodes2)
vsim22 <- vsim(vc2,nodes3)

Aspace <- aproxdef(22,x[1]*1.05,x[600]*0.99,0.02)
vc3 <- chebc(Aspace,nodes,sig03(nodes))
vsim3 <- vsim(vc3,nodes2)
vsim32 <- vsim(vc3,nodes3)

## Unit adjustment: $1,000/lb to $/lb
vsim0p1$shadowp <- vsim0p1$shadowp*1000
vsim0p2$shadowp <- vsim0p2$shadowp*1000
vsim1$shadowp <- vsim1$shadowp*1000
vsim2$shadowp <- vsim2$shadowp*1000
vsim3$shadowp <- vsim3$shadowp*1000
vsim12$shadowp <- vsim12$shadowp*1000
vsim22$shadowp <- vsim22$shadowp*1000
vsim32$shadowp <- vsim32$shadowp*1000


####################################################
## Figure 7(a)
####################################################
set1 <- data.frame(stock=sto00$stock[1:160],vfun=sto00$vfun[1:160],key='determ1')
set2 <- data.frame(stock=sto00$stock[161:600],vfun=sto00$vfun[161:600],key='determ2')
set3 <- data.frame(stock=sto01$stock,vfun=sto01$vfun,key='sto01')
set4 <- data.frame(stock=sto02$stock,vfun=sto02$vfun,key='sto02')
set5 <- data.frame(stock=sto03$stock,vfun=sto03$vfun,key='sto03')

setall <- rbind(set2,set3,set4,set5)

colnames(setall) <- c('stock','vfun','key')
setall$key <- factor(setall$key, 
                     levels = c('determ2','sto01',
                                'sto02','sto03'))

ggplot(data=setall,aes(x=stock,y=vfun)) + 
  coord_cartesian(xlim=c(0,350), ylim=c(0,12)) +
  geom_line(lwd=0.8,aes(col=key,linetype=key)) + 
  scale_linetype_manual(name='',values=c(1,2,2,2),
                        labels = c('Deterministic',
                                   expression( paste("Stochastic (",sigma," = ","0.1)")),
                                   expression( paste("Stochastic (",sigma," = ","0.2)")),
                                   expression( paste("Stochastic (",sigma," = ","0.3)")))) +
  #  scale_color_manual(name='',values=c('black','grey40','grey60','grey80'),
  scale_color_manual(name='', values = c('#000000','#3366CC','#0000FF','#9900FF'),
                     labels = c('Deterministic',
                                expression( paste("Stochastic (",sigma," = ","0.1)")),
                                expression( paste("Stochastic (",sigma," = ","0.2)")),
                                expression( paste("Stochastic (",sigma," = ","0.3)")))) +
  geom_line(data=set1,lwd=0.8,col='#000000') +
  geom_vline(xintercept=ss1,linetype=2) + 
  geom_vline(xintercept=ss2,linetype=2) + 
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.box = 'vertical',
        legend.position = c(0.6,0.2),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        legend.text = element_text(size=12),
        legend.text.align = 0) +
  annotate(geom="text", x=ss1, y=12.2, label="Steady State 1",color="black") +
  annotate(geom="text", x=ss2, y=12.2, label="Steady State 2",color="black") +
  labs(x='Stock (Million lb)',y='Intertemproal Welafare (V in $ Billions)')

####################################################
## Figure 7(b)
####################################################
set1 <- data.frame(stock=vsim0p1$stock,shadowp=vsim0p1$shadowp,key='determ1')
set2 <- data.frame(stock=vsim0p2$stock,shadowp=vsim0p2$shadowp,key='determ2')
set3 <- data.frame(stock=vsim1$stock,shadowp=vsim1$shadowp,key='sto01')
set4 <- data.frame(stock=vsim2$stock,shadowp=vsim2$shadowp,key='sto02')
set5 <- data.frame(stock=vsim3$stock,shadowp=vsim3$shadowp,key='sto03')

setall <- rbind(set2,set3,set4,set5)

colnames(setall) <- c('stock','shadowp','key')
colnames(set1) <- c('stock','shadowp','key')

setall$shadowp <- setall$shadowp

colnames(setall) <- c('stock','shadowp','key')
colnames(set1) <- c('stock','shadowp','key')
setall$key <- factor(setall$key, 
                     levels = c('determ2','sto01',
                                'sto02','sto03'))


ggplot(data=setall,aes(x=stock,y=shadowp)) + 
  #  coord_cartesian(xlim=c(0,350), ylim=c(0,15)) +
  geom_line(lwd=0.8,aes(col=key,linetype=key)) + 
  scale_linetype_manual(name='',values=c(1,2,2,2),
                        labels = c('Deterministic',
                                   expression( paste("Stochastic (",sigma," = ","0.1)")),
                                   expression( paste("Stochastic (",sigma," = ","0.2)")),
                                   expression( paste("Stochastic (",sigma," = ","0.3)")))) +
  #  scale_color_manual(name='',values=c('black','grey40','grey60','grey80'),
  scale_color_manual(name='', values = c('#000000','#3366CC','#0000FF','#9900FF'),
                     labels = c('Deterministic',
                                expression( paste("Stochastic (",sigma," = ","0.1)")),
                                expression( paste("Stochastic (",sigma," = ","0.2)")),
                                expression( paste("Stochastic (",sigma," = ","0.3)")))) +
  geom_line(data=set1,lwd=0.8,col='#000000') +
  geom_vline(xintercept=ss1,linetype=2) + 
  geom_vline(xintercept=ss2,linetype=2) + 
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.box = 'vertical',
        legend.position = c(0.6,0.3),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        legend.text = element_text(size=12),
        legend.text.align = 0) +
  annotate(geom="text", x=ss1, y=-30, label="Steady State 1",color="black") +
  annotate(geom="text", x=ss2, y=-30, label="Steady State 2",color="black") +
  labs(x='Stock (Million lb)',y='Accounting Price (p in $/lb)')

####################################################
## Figure 7(c)
####################################################
set1 <- data.frame(stock=vsim0p1$stock,shadowp=vsim0p1$shadowp,key='determ1')
set2 <- data.frame(stock=vsim0p2$stock,shadowp=vsim0p2$shadowp,key='determ2')
set3 <- data.frame(stock=vsim12$stock[1:100],shadowp=vsim12$shadowp[1:100],key='sto11')
set4 <- data.frame(stock=vsim12$stock[101:200],shadowp=vsim12$shadowp[101:200],key='sto12')
set5 <- data.frame(stock=vsim22$stock[1:100],shadowp=vsim22$shadowp[1:100],key='sto21')
set6 <- data.frame(stock=vsim22$stock[101:200],shadowp=vsim22$shadowp[101:200],key='sto22')
set7 <- data.frame(stock=vsim32$stock[1:100],shadowp=vsim32$shadowp[1:100],key='sto31')
set8 <- data.frame(stock=vsim32$stock[101:200],shadowp=vsim32$shadowp[101:200],key='sto32')

set3$diff <- set3$shadowp - set1$acc.price1
set4$diff <- set4$shadowp - set2$acc.price1
set5$diff <- set5$shadowp - set1$acc.price1
set6$diff <- set6$shadowp - set2$acc.price1
set7$diff <- set7$shadowp - set1$acc.price1
set8$diff <- set8$shadowp - set2$acc.price1


setall <- rbind(set3,set5,set7)
colnames(setall) <- c('stock','shadowp','key','diff')
setall$key <- factor(setall$key, levels = c('sto11','sto21','sto31'))


ggplot(data=setall,aes(x=stock,y=diff)) + 
  coord_cartesian(xlim=c(0,350), ylim=c(-500,400)) +
  geom_line(lwd=0.8,aes(col=key,linetype=key)) + 
  scale_linetype_manual(name='',values=c(2,2,2),
                        labels = c(expression( paste("Stochastic (",sigma," = ","0.1)")),
                                   expression( paste("Stochastic (",sigma," = ","0.2)")),
                                   expression( paste("Stochastic (",sigma," = ","0.3)")))) +
  scale_color_manual(name='',values=c('#3366CC','#0000FF','#9900FF'),
                     labels = c(expression( paste("Stochastic (",sigma," = ","0.1)")),
                                expression( paste("Stochastic (",sigma," = ","0.2)")),
                                expression( paste("Stochastic (",sigma," = ","0.3)")))) +
  geom_line(data=set4,lty=2,lwd=0.8,col='#3366CC') +
  geom_line(data=set6,lty=2,lwd=0.8,col='#0000FF') +
  geom_line(data=set8,lty=2,lwd=0.8,col='#9900FF') +
  geom_vline(xintercept=ss1,linetype=2) + 
  geom_vline(xintercept=ss2,linetype=2) + 
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.box = 'vertical',
        legend.position = c(0.6,0.8),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        legend.text = element_text(size=12),
        legend.text.align = 0) +
  annotate(geom="text", x=ss1, y=400, label="Steady State 1",color="black") +
  annotate(geom="text", x=ss2, y=400, label="Steady State 2",color="black") +
  labs(x='Stock (Million lb)',y='Acc. Price Diff. from Deterministic ($/lb)')