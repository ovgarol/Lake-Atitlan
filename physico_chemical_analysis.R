#!/usr/bin/env Rscript
# SPDX-FileCopyrightText: 2025 Ovidio Garcia-Oliva
# SPDX-License-Identifier: CC-BY-4.0
# SPDX-FileContributor Ovidio Garcia-Oliva <ovidio.garcia@hereon.de
#
# This script analyzes physicochemical variables of lake. 

library(cmocean)
library(trend)
library(latex2exp)
library(lubridate)
library(mblm)
library(scales)
library(rLakeAnalyzer)

date.lim = as.Date(c('2010-01-01','2025-01-01'))

plot = function(...) graphics::plot(...,     xaxs='i',yaxs='i')

title = function(...) graphics::title(...,font.main=1)

mean = function(...)base::mean(...,na.rm=T)

diff0 = function(x,...){
  y=diff(x,...)
  return(c(y,y[length(y)]))
}

column.filler = function(db){
  col.z = NULL
  col.depth = NULL
  col.date = NULL
  for(i in unique(db$date)){
    column = db[db$date==i,]
    if(length(column$depth)<2)next
    #x = min(column$depth):max(column$depth)
    x = seq(round(min(column$depth)),round(max(column$depth)+0.1),1)
    y = approx(column$depth,column$z,x)$y
    col.depth = c(col.depth,x)
    col.z = c(col.z,y)
    col.date = c(col.date,rep(i,length(x)))
  }
  db.r = as.data.frame(list(date=as.Date(col.date,'1970-01-01'),depth=col.depth,z=col.z))
  return(db.r)
}

get.variable = function(var.name){
  CWG.xx = subset(CWG.all,select = c('date','depth',var.name))
  CWG.xx$z = CWG.xx[[3]]
  CWG.xx = na.omit(CWG.xx)
  CWG.xx = column.filler(CWG.xx)
  return(CWG.xx)
}

date.axis = function(...){
  abline(v=as.Date(paste0(2011:2024,'-01-01')),col=c(rep('gray95',4),'gray95',rep('gray95',4),'gray95',rep('gray95',3)))
  axis(side=1, at =as.Date(paste0(c(2010,2012,2014,2016,2018,2020,2022,2024),'-01-01')),labels=c(2010,2012,2014,2016,2018,2020,2022,2024))
}

## lake bathymetry after Reyes Morales et al., 2017

bthD = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330)
bthA = c(125780000,120911738,117933973,114988032,112438621,109757028,105889298,103898886,101762812,99483844,97035391,94336082,91454624,88473435,85322780,82529547,80012657,77529869,74974090,72289184,69502890,66784155,63997799,61396806,58672558,55857048,52520835,49523589,46370480,42618706,38108775,31682062,14877818,0)

bottomd = 100
surfaced = 30

###

weather = read.csv("./data/weather.csv",comment.char = '#')
weather$date = as.Date(weather$date)

ENSO.orig = read.csv("./data/ENSO.csv",comment.char = '#')
ENSO.orig$date = ymd(paste(ENSO.orig$year,ENSO.orig$month,'15'))

###

CWG.all = read.csv("./data/Lake_Atitlan_water_quality_CWG_2010-2024.csv",comment.char = '#')
CWG.all$date = ymd(CWG.all$date)

CWG.temp = subset(CWG.all,select = c(date,depth,temperature))
CWG.temp = na.omit(CWG.temp)
CWG.temp$z = CWG.temp$temperature
CWG.temp.0 = CWG.temp # used for uncertainty calculation
CWG.temp = column.filler(CWG.temp)

CWG.do = subset(CWG.all,select = c(date,depth,dissolved_oxygen))
CWG.do = na.omit(CWG.do)
CWG.do$z = CWG.do$dissolved_oxygen
CWG.do = column.filler(CWG.do)

CWG.summary = NULL
ldepth = 10
ldepth2 = 50

for(i in sort(unique(CWG.temp$date))){
  wtr = CWG.temp[CWG.temp$date==i,]
  wtr.0 = CWG.temp.0[CWG.temp.0$date==i,]
  tclina = thermo.depth(wtr$z[wtr$depth>ldepth],wtr$depth[wtr$depth>ldepth],seasonal = T,Smin=.1)
  SS = schmidt.stability(wtr$z,wtr$depth,bthA,bthD)
  SS.unc = 2*schmidt.stability((0.02+abs(diff0(wtr.0$z)))/diff0(wtr.0$depth)*min(mean(diff(wtr.0$depth)),50),wtr.0$depth,bthA,bthD) 
  max.depth = max(wtr$depth)
  ie = internal.energy(wtr$z,wtr$depth, bthA,bthD)
  lake.temperature = whole.lake.temperature(wtr$z,wtr$depth,bthA,bthD)
  max.lake.temperature = max(wtr$z,na.rm=T)
  lake.temp.unc = 2*min(mean(diff(wtr.0$depth)),50)*abs(whole.lake.temperature((0.02+abs(diff0(wtr.0$z)))/diff0(wtr.0$depth),wtr.0$depth,bthA,bthD))+0.02
  metal = meta.depths(wtr$z,wtr$depth,slope=0.01)
  min.temperature = min(wtr$z)
  oxy = CWG.do[CWG.do$date==i,]
  oxycline = thermo.depth(oxy$z[oxy$depth>ldepth],oxy$depth[oxy$depth>ldepth],seasonal=T,Smin=.1)

  CWG.summary = rbind(CWG.summary,c(i ,tclina,SS,max.depth,ie,lake.temperature,lake.temp.unc,SS.unc,metal,min.temperature,oxycline,max.lake.temperature))
}

CWG = NULL
CWG$date = as.Date(CWG.summary[,1],origin = '1970-01-01')
CWG$tclina = CWG.summary[,2]
CWG$SS = CWG.summary[,3]/1e3
CWG$max.depth = CWG.summary[,4]
CWG$ie = CWG.summary[,5]/1e9
CWG$lake.temperature = CWG.summary[,6]
CWG$lake.temperature.unc = CWG.summary[,7]
CWG$SS.unc = CWG.summary[,8]/1e3
CWG$meta1 =  CWG.summary[,9]
CWG$meta2 =  CWG.summary[,10]
CWG$min.temp =  CWG.summary[,11]
CWG$oxycline =  CWG.summary[,12]
CWG$max.temp =  CWG.summary[,13]

CWG = as.data.frame(CWG)
CWG$Group.1 = paste(year(CWG$date),month(CWG$date))
CWG$tclina = ifelse(CWG$tclina>ldepth&CWG$max.depth>=ldepth2,CWG$tclina,NA)
CWG$ie = ifelse(CWG$max.depth>ldepth2,CWG$ie,NA)
CWG$SS = ifelse(CWG$max.depth>ldepth2,CWG$SS,NA)
CWG$lake.temperature = ifelse(CWG$max.depth>ldepth2,CWG$lake.temperature,NA)
CWG$min.temp = ifelse(CWG$max.depth>ldepth2,CWG$min.temp,NA)

CWG$meta1 = ifelse(CWG$meta1>ldepth,CWG$meta1,NA)
CWG$meta2 = ifelse(CWG$meta2>ldepth,CWG$meta2,NA)

################################################################################
pdf('./fig/Fig_1.pdf',w=8.5,h=11)
par(mfrow=c(6,1),mai=0.1+c(0.4,0.5,0.1,0.1),oma=c(2.5,.5,0,0),lend=1)
par(bty='o',las=1,family='Helvetica')

plot(ENSO.orig$date,ENSO.orig$ONI,
     col='black',
     lwd=1,
     type='n',
     xlim=date.lim,
     ylim=c(-3,3),
     ylab='°C',
     xlab='',xaxt='n'
     )
title(main='a. Oceanic Niño Index, ONI',adj=0.)
abline(h=0.,lty=2)
lines(ENSO.orig$date,ENSO.orig$ONI,lwd=1,type='l',col='black')
date.axis()



plot(weather$date,weather$ec.t.mean,
     type='n',xlim=date.lim,
     ylim=c(-1.5,1.5),pch=19,col='gray95',
     ylab='°C',
     xlab='',xaxt='n'
     )
title(main='b. air temperature anomaly (monthly average)',adj=0.)
abline(h=0,lty=2)
date.axis()

#lines(c(clima$date,sa.mean$date),sam-sa.mean0,col=alpha('darkgray',0.5),type='p',pch=15,lwd=1,cex=0.5)
#lines(c(clima$date,ec.mean$date),ecm-ec.mean0,col=alpha('black',0.5),type='p',pch=15,lwd=1,cex=0.5)

lines(supsmu(weather$date,weather$sa.t.mean-mean(weather$sa.t.mean),span=0.005),col='black',type='l',pch=15,lwd=1)
lines(supsmu(weather$date,weather$ec.t.mean-mean(weather$ec.t.mean),span=0.005),col='darkgray',type='l',pch=15,lwd=1)
legend('topright',
       bg=alpha('white',0.66),
       bty='n',
       horiz=F,
       cex=1,
       legend=c('Santiago','El Capitan'),
       pch=15,
       pt.cex = 2,
       col=alpha(c('darkgray','black'),1.5)
)




plot(weather$date,weather$ec.wind.vel*0.2778,
     type='n',
     xlim=date.lim,
     ylim=c(0,4),
     pch=19,col='black',lwd=2,
     ylab=latex2exp::TeX('m s$^{-1}$'),
     xlab='',xaxt='n'
)
title(main='c. wind velocity (monthly average)',adj=0.)
date.axis()
lines(supsmu(weather$date,weather$sa.wind.vel*0.2778,span = 0.05),col='darkgray')
lines(supsmu(weather$date,weather$ec.wind.vel*0.2778,span = 0.05),col='black')
legend('topright',
       bg=alpha('white',0.66),
       bty='n',
       horiz=F,
       cex=1,
       legend=c('Santiago','El Capitan'),
       pt.cex = 2,
       pch=15,
       col=alpha(c('darkgray','black'),1.5)
)



plot(lake.temperature~date,data=CWG,type='p',xlim=date.lim,col='black',pch=19,ylim=c(20,23),
     ylab='°C',lwd=2,cex=0.25,
     xlab='',xaxt='n',yaxt='n')
title(main='d. whole lake temperature',adj=0.)
date.axis()
axis(side=2,20:23)

for(i in 1:length(CWG$date))lines(CWG$date[i]+c(0,0),CWG$lake.temperature.unc[i]*c(-1,1)+CWG$lake.temperature[i],col='black',lwd=0.5)
for(i in 1:length(CWG$date))points(CWG$date[i]+c(0,0),CWG$lake.temperature.unc[i]*c(-1,1)+CWG$lake.temperature[i],col='black',pch='_')
summary(lm(lake.temperature~date,data=CWG[year(CWG$date)<=2015,]))

plot(tclina~date,data=CWG,
     type='p',
     xlim=date.lim,
     pch=19,
     ylab='m',
     lwd=2,
     cex=0.25,
     xlab='',
     xaxt='n',
     yaxt='l',
     ylim=c(80,0)
)
title(main='e. thermocline depth',adj=0.)
date.axis()

plot(SS~date,data=CWG,type='n',xlim=date.lim,pch=19,ylim=c(-10,50),
     ylab=latex2exp::TeX('kJ m$^{-2}$'),lwd=2,cex=0.25,
     xlab='',xaxt='n',yaxt='l')
title(main='f. Schmidt stability',adj=0.)
abline(h=1,lty=2)
date.axis()

for(i in 1:length(CWG$date))lines(CWG$date[i]+c(0,0),2*CWG$SS.unc[i]*c(-1,1)+CWG$SS[i],col='black',lwd=0.5)
for(i in 1:length(CWG$date))points(CWG$date[i]+c(0,0),2*CWG$SS.unc[i]*c(-1,1)+CWG$SS[i],col='black',pch='_')
points(SS~date,data=CWG,type='p',xlim=date.lim,cex=0.5+.5*as.numeric(CWG$SS<1),pch=19,col=c('black','tomato')[1+as.numeric(CWG$SS<1)])

dev.off()

################################################################################
pdf('./fig/Fig_2.pdf',w=8.5,h=11)
par(mfrow=c(6,1),mai=0.1+c(0.4,0.5,0.1,0.1),oma=c(2.5,.5,0,0),lend=1)
par(bty='o',las=1,family='Helvetica')

CWG.do.bottom = aggregate(CWG.do[CWG.do$depth>bottomd,],by=list(CWG.do$date[CWG.do$depth>bottomd]),mean)

plot(CWG.do.bottom$date,CWG.do.bottom$z,
     xlim=date.lim,type='n',pch=19,ylim=c(0,5),
     ylab=latex2exp::TeX('mg L$^{-1}$'),lwd=2,cex=0.25,
     xlab='',xaxt='n',yaxt='l')
title(main='a. bottom dissolved oxygen (depth > 100 m)',adj=0.)
date.axis()
lines(smooth.spline(CWG.do.bottom$date,CWG.do.bottom$z),type='l',pch=19,cex=1,lwd=2)
points(CWG.do.bottom$date,CWG.do.bottom$z,type='p',pch=19,cex=1,col=alpha('black',0.5))

CWG.nox = get.variable('nitrate_plus_nitrite')
CWG.nox$z = (CWG.nox$z*(CWG.nox$z>2)+1)
CWG.nox.bottom = aggregate(CWG.nox[CWG.nox$depth>bottomd,],by=list(CWG.nox$date[CWG.nox$depth>bottomd]),mean)
CWG.nox.surface = aggregate(CWG.nox[CWG.nox$depth<surfaced,],by=list(CWG.nox$date[CWG.nox$depth<surfaced]),mean)

CWG.nh4 = get.variable('ammonia')
CWG.nh4$z = (CWG.nh4$z*(CWG.nh4$z>2)+1)
CWG.nh4.bottom = aggregate(CWG.nh4[CWG.nh4$depth>bottomd,],by=list(CWG.nh4$date[CWG.nh4$depth>bottomd]),mean)
CWG.nh4.surface = aggregate(CWG.nh4[CWG.nh4$depth<surfaced,],by=list(CWG.nh4$date[CWG.nh4$depth<surfaced]),mean)

CWG.po4 = get.variable('ortho.phosphate')
CWG.po4$z = (CWG.po4$z*(CWG.po4$z>2)+1)
CWG.po4.bottom = aggregate(CWG.po4[CWG.po4$depth>bottomd,],by=list(CWG.po4$date[CWG.po4$depth>bottomd]),mean)
CWG.po4.surface = aggregate(CWG.po4[CWG.po4$depth<surfaced,],by=list(CWG.po4$date[CWG.po4$depth<surfaced]),mean)

CWG.tp = get.variable('total_phosphorus')
CWG.tp$z = (CWG.tp$z*(CWG.tp$z>2)+1)
CWG.tp.bottom = aggregate(CWG.tp[CWG.tp$depth>bottomd,],by=list(CWG.tp$date[CWG.tp$depth>bottomd]),mean)
CWG.tp.surface = aggregate(CWG.tp[CWG.tp$depth<surfaced,],by=list(CWG.tp$date[CWG.tp$depth<surfaced]),mean)

CWG.din.bottom = merge(CWG.nh4.bottom,CWG.nox.bottom,by="date",suffixes=c('.nh4','.nox'))
CWG.din.bottom$z = CWG.din.bottom$z.nh4+CWG.din.bottom$z.nox

CWG.din.surface = merge(CWG.nh4.surface,CWG.nox.surface,by="date",suffixes=c('.nh4','.nox'))
CWG.din.surface$z = CWG.din.surface$z.nh4+CWG.din.surface$z.nox


plot(CWG.po4.bottom$date,CWG.po4.bottom$z,xlim=date.lim,type='n',ylim=c(10,200.),
     ylab=latex2exp::TeX('$\\mu$g-P L$^{-1}$'),lwd=2,cex=0.25,log='y',
     xlab='',xaxt='n',yaxt='l')
title(main='b. bottom phosphorus (depth > 100 m)',adj=0.)
date.axis()

points(CWG.po4.bottom$date,CWG.po4.bottom$z,col=alpha('tomato',0.5),pch=19,type='p',cex=1.)
lines(smooth.spline(CWG.po4.bottom$date,(CWG.po4.bottom$z)),col='tomato',lwd=2)
points(CWG.tp.bottom$date,(CWG.tp.bottom$z),col=alpha('darkred',0.5),pch=15,type='p',cex=1.)
lines(smooth.spline(CWG.tp.bottom$date,(CWG.tp.bottom$z)),col='darkred',lwd=2)

legend('bottomleft',
       bg=alpha('white',0.66),
       bty='n',
       horiz=T,
       legend=c('SRP','TP'),
       pch=c(19,15),
       col=alpha(c('tomato','darkred'),0.5)
)



plot(CWG.nox.bottom$date,(CWG.nox.bottom$z),xlim=date.lim,type='n',ylim=c(1,500),
     ylab=latex2exp::TeX('$\\mu$g-N L$^{-1}$'),lwd=2,cex=0.25,log='y',
     xlab='',xaxt='n',yaxt='l')
title(main='c. bottom inorganic nitrogen species (depth > 100 m)',adj=0.)
date.axis()

points(CWG.nox.bottom$date,(CWG.nox.bottom$z),col=alpha('darkcyan',0.5),pch=19,type='p',cex=1.)
lines(smooth.spline(CWG.nox.bottom$date,(CWG.nox.bottom$z)),col='darkcyan',lwd=2)
points(CWG.nh4.bottom$date,(CWG.nh4.bottom$z),col=alpha('aquamarine3',0.4),pch=20,type='p',cex=1.)
lines(smooth.spline(CWG.nh4.bottom$date,(CWG.nh4.bottom$z)),col='aquamarine3',lwd=2)

legend('bottomleft',
       bg=alpha('white',0.66),
       bty='n',
       horiz=T,
       legend=c(latex2exp::TeX('NO$_2$+NO$_3$'),latex2exp::TeX('NH$_4$')),
       pch=c(19,20),
       col=alpha(c('darkcyan','aquamarine3'),0.5)
)



plot(CWG.din.bottom$date,(CWG.din.bottom$z.nh4/CWG.din.bottom$z),
     xlim=date.lim,
     type='n',
     ylim=c(0,1),
     pch=19,
     ylab='',lwd=2,cex=0.25,
     xlab='',xaxt='n',yaxt='l'
     )
title(main=latex2exp::TeX("d. bottom NH$_4$-to-DIN ratio (depth > 100 m)"),adj=0.)
date.axis()

points(CWG.din.bottom$date,(CWG.din.bottom$z.nh4/CWG.din.bottom$z),type='p',pch=19,cex=1,col=alpha('black',0.5))
lines(smooth.spline(CWG.din.bottom$date,(CWG.din.bottom$z.nh4/CWG.din.bottom$z)),type='l',pch=19,cex=1,lwd=2)




plot(CWG.nox.surface$date,CWG.nox.surface$z,xlim=date.lim,type='n',ylim=c(1,100),
     ylab=latex2exp::TeX('$\\mu$g-N/-P L$^{-1}$'),lwd=2,cex=0.25,log='y',
     xlab='',xaxt='n',yaxt='l')
title(main='e. surface nutrient (depth < 30 m)',adj=0.)
date.axis()

points(CWG.po4.surface$date,(CWG.po4.surface$z),col=alpha('tomato',0.5),pch=19,type='p',cex=1)
lines(smooth.spline(CWG.po4.surface$date,(CWG.po4.surface$z)),col='tomato',lwd=2)

points(CWG.din.surface$date,(CWG.din.surface$z),col=alpha('darkcyan',0.5),pch=19,type='p',cex=1.)
lines(smooth.spline(CWG.din.surface$date,(CWG.din.surface$z)),col='darkcyan',lwd=2)

legend('bottomleft',
       bg=alpha('white',0.66),
       bty='n',
       horiz=T,
       legend=c('SRP','DIN'),
       pch=c(19,19),
       col=alpha(c('tomato','darkcyan'),0.5)
)

CWG.rnp=merge(CWG.po4.surface,CWG.din.surface,by='date')
CWG.rnp$z=CWG.rnp$z.y/CWG.rnp$z.x


plot(CWG.rnp$date,log10(CWG.rnp$z),
     xlim=date.lim,type='n',ylim=c(0.03,30),pch=19,
     ylab='',lwd=2,cex=0.25,log='y',
     xlab='',xaxt='n',yaxt='l'
     )
title(main='f. surface DIN-to-SRP ratio (depth < 30 m)',adj=0.)
date.axis()

abline(h=1,lty=2)
abline(h=(7.22)*c(1),col='tomato',lty=1)
points(CWG.rnp$date,(CWG.rnp$z),type='p',pch=19,cex=0.75,col=alpha('black',0.5))
lines(smooth.spline(CWG.rnp$date,(CWG.rnp$z)),type='l',pch=19,cex=1,lwd=2)

dev.off()

###




write.csv(subset(CWG.nh4.bottom,select=c(date,z)),'./out/CWG_bottom_NH4.csv',row.names = F)
write.csv(subset(CWG.nox.bottom,select=c(date,z)),'./out/CWG_bottom_NOx.csv',row.names = F)
write.csv(subset(CWG.din.bottom,select=c(date,z)),'./out/CWG_bottom_DIN.csv',row.names = F)
write.csv(subset(CWG.po4.bottom,select=c(date,z)),'./out/CWG_bottom_SRP.csv',row.names = F)
write.csv(subset(CWG.tp.bottom,select=c(date,z)),'./out/CWG_bottom_TP.csv',row.names = F)
write.csv(subset(CWG.do.bottom,select=c(date,z)),'./out/CWG_bottom_DO.csv',row.names = F)

write.csv(subset(CWG.nh4.surface,select=c(date,z)),'./out/CWG_surface_NH4.csv',row.names = F)
write.csv(subset(CWG.nox.surface,select=c(date,z)),'./out/CWG_surface_NOx.csv',row.names = F)
write.csv(subset(CWG.din.surface,select=c(date,z)),'./out/CWG_surface_DIN.csv',row.names = F)
write.csv(subset(CWG.po4.surface,select=c(date,z)),'./out/CWG_surface_SRP.csv',row.names = F)
write.csv(subset(CWG.tp.surface,select=c(date,z)),'./out/CWG_surface_TP.csv',row.names = F)
write.csv(subset(CWG.rnp,select=c(date,z)),'./out/CWG_DIN2SRP.csv',row.names = F)

###################
ENSO.orig$Group.1 = paste(ENSO.orig$year,ENSO.orig$month)
ENSO.orig$date = ENSO.orig$date.y
ENSO = merge(ENSO.orig,CWG,by='Group.1')
ENSO = ENSO[order(as.POSIXct(ENSO$date)),]

pdf('./fig/Fig_S2.pdf',w=7,h=4)
par(mfrow=c(2,3),mai=0.75*c(0.8, 0.8, 0.2, 0.2),las=1,oma=c(1,0,0,0))
layout(matrix(c(1,2,3,4,5,6),nrow=2))

ENSO = ENSO[!is.na(ENSO$ie) & ENSO$max.depth>70,]
ENSO = ENSO[ENSO$year>2010,]

s = lm(max.temp~date+ONI+as.factor(month),data=ENSO,weights = 1/ENSO$lake.temperature.unc)
summary(s)
s0 = glm(max.temp~date+as.factor(month),data=ENSO,weights = 1/ENSO$lake.temperature.unc)
summary(s0)

anova(s,s0)
AIC(s)
AIC(s0)

rr = summary(s)
rr$coefficients[2,1]*365.25*10
rr$coefficients[2,2]*365.25*10
rr$coefficients[2,4]

rr$coefficients[3,1]
rr$coefficients[3,2]
rr$coefficients[3,4]

p = predict(s,newdata = ENSO)
p0 = predict(s0,newdata = ENSO)

plot(ENSO$date,ENSO$max.temp,
     pch=19,
     col='gray',
     xlab='',
     ylab='temperature °C',
     xaxt='n'
)
date.axis()
title(main='a. surface lake temperature',adj=0)
lines(ENSO$date,p,pch=19,col='black',cex=0.2)
#lines(ENSO$date,p0,pch=19,col='green')

plot(ENSO$max.temp,p,ylim=c(20,25),
     xlim=c(20,25),
     pch=19,
     col='gray',
     xlab='',
     ylab='modeled temperature °C'
)
title(main='d.',adj=0)

#points(ENSO$max.temp,p0,col='green')
abline(a=0,b=1,lty=3)


###################
ENSO = ENSO[!is.na(ENSO$ie) & ENSO$max.depth>70,]
ENSO = ENSO[ENSO$year > 2018,]

s = lm(lake.temperature~date+ONI+as.factor(month),data=ENSO,weights = 1/ENSO$lake.temperature.unc)
summary(s)
s0 = lm(lake.temperature~date+as.factor(month),data=ENSO,weights = 1/ENSO$lake.temperature.unc)
summary(s0)

anova(s,s0)
AIC(s)
AIC(s0)

rr=summary(s)
rr$coefficients[2,1]*365.25*10
rr$coefficients[2,2]*365.25*10
rr$coefficients[2,4]

rr$coefficients[3,1]
rr$coefficients[3,2]
rr$coefficients[3,4]

p = predict(s,newdata = ENSO)
p0 = predict(s0,newdata = ENSO)

plot(ENSO$date,ENSO$lake.temperature,
     xlab='',ylab='',
     pch=19,
     col='gray',
     xaxt='n'
)
title(main='b. whole-lake temperature',adj=0)
date.axis()
lines(ENSO$date,p,pch=19,col='black')
#lines(ENSO$date,p0,pch=19,col='green')

plot(ENSO$lake.temperature,p,
     pch=19,col='gray',
     ylim=c(20.3,21.3),xlim=c(20.3,21.3),
     xlab='',ylab='')
title(main='e.',adj=0)
mtext(side=1, 'observed temperature °C',outer=T,cex=0.8,line=-2)
#points(ENSO$lake.temperature,p0,col='green')
abline(a=0,b=1,lty=3)

###################

s = glm(min.temp~date+ONI+as.factor(month),data=ENSO,weights = 1/ENSO$lake.temperature.unc)
summary(s)
s0 = glm(min.temp~date+as.factor(month),data=ENSO,weights = 1/ENSO$lake.temperature.unc)
summary(s0)

anova(s,s0)
AIC(s)
AIC(s0)

rr=summary(s)
rr$coefficients[2,1]*365.25*10
rr$coefficients[2,2]*365.25*10
rr$coefficients[2,4]

rr$coefficients[3,1]
rr$coefficients[3,2]
rr$coefficients[3,4]

p = predict(s,newdata = ENSO)
p0 = predict(s0,newdata = ENSO)

plot(ENSO$date,ENSO$min.temp,
     pch=19,col='gray',
     xlab='',ylab='',xaxt='n')
date.axis()

title(main='c. minimum lake temperature',adj=0)
lines(ENSO$date,p,pch=19,col='black')
#lines(ENSO$date,p0,pch=19,col='green')

plot(ENSO$min.temp,p,ylim=c(20.2,20.5),xlim=c(20.2,20.5),
     pch=19,col='gray',
     xlab='',ylab='')
#points(ENSO$min.temp,p0,col='green')
abline(a=0,b=1,lty=3)
title(main='f.',adj=0)


##########

ENSO = ENSO[!is.na(ENSO$ie) & ENSO$max.depth>70,]
ENSO = ENSO[ENSO$year>2018,]

s = lm(SS~date+ONI+as.factor(month),data=ENSO)#,weights = 1/abs(ENSO$SS.unc))
summary(s)
s0 = glm(SS~date+as.factor(month),data=ENSO)#,weights = 1/abs(ENSO$SS.unc))
summary(s0)

anova(s,s0)
AIC(s)
AIC(s0)

rr=summary(s)
rr$coefficients[2,1]*365.25*10
rr$coefficients[2,2]*365.25*10
rr$coefficients[2,4]

rr$coefficients[3,1]
rr$coefficients[3,2]
rr$coefficients[3,4]

dev.off()

################################################################################
sc

ENSO.orig$Group.1 = paste(ENSO.orig$year,ENSO.orig$month)
ENSO.orig$date = ENSO.orig$date.y
ENSO = merge(ENSO.orig,CWG,by='Group.1')

weather$Group.1 = paste(weather$year,weather$month)
all.data = merge(weather,ENSO,by='Group.1')

library(corrgram)
sub.all.data = subset(all.data,select = c(ec.t.mean,sa.t.mean,ec.wind.vel,sa.wind.vel,ONI))
colnames(sub.all.data) = c('air Temp. EC','air Temp. SA','wind vel. EC','wind vel. SA','ONI')

pdf('./fig/Fig_S3.pdf',w=6,h=6)
corrgram(sub.all.data,
         order=TRUE, 
         pch=20,
         lower.panel=panel.pts, 
         upper.panel=panel.cor,
         #diag.panel=panel.density
         )
title(main="Correlation of weather and climate variables",adj=0)
dev.off()

