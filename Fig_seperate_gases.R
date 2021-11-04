library("tidyverse")
library("lubridate")

#data import ####

#ABISS data
ox1 <- read_csv("~/Dropbox/Harvard/ABISS/Analysis/FinalData/optode1.csv")
ox2 <- read_csv("~/Dropbox/Harvard/ABISS/Analysis/FinalData/optode2.csv")
rga <- read_csv("~/Dropbox/Harvard/ABISS/Analysis/FinalData/rga.csv")

#calibrated ch4 data
conc_dat0 <- read_csv("~/Dropbox/Harvard/ABISS/Analysis/FinalData/gas_conc.csv")


# timeline ####
#event times in UTC, from Dan's handwritten notes
launch = ymd_hm(201809151300)
ROV_pu = ymd_hm(201809151631,201809151845) #lander picked up by ROV
ROV_pd = ymd_hm(201809151726,201809151857) #lander put down by ROV
ROV_clear = ymd_hm(201809152000) #ROV cleared the area 
lander_pos = ymd_hm(201809151910,201809151936)
wand_pos = ymd_hm(201809151942,201809151951)
power_blip = ymd_hm(201809151735)
bubble_fly = ymd_hm(201809151921) #flew through bubble stream
wand_in_pos = ymd_hm(201809151951) #start of real data

# inferred directly from log...
power_on = ox1$datetime[1]
rga_start = rga$datetime[1]
warmed_up = rga$datetime[129] #give it 30 mins

# indices of phases of deployment 
wp <- c(13332:13451,13870:13943) #weird peaks
warmup = which(rga$datetime < warmed_up)
background = which(rga$datetime >= power_blip+20*60 & rga$datetime < ROV_pu[2]) #lander in the mud while ROV searched for site, starting 20 mins after restart
active_onsite = which(rga$datetime >= wand_in_pos)
positioning = which((rga$datetime >= ROV_pu[1] & rga$datetime < ROV_pd[1]) | (rga$datetime >= ROV_pu[2] & rga$datetime < wand_pos[2]))
ebullitionA = which(rga$datetime >= wand_in_pos & rga$datetime < ymd_hm(201809160210)) #ymd_hm(201809152125)
ebullitionB = which(rga$datetime >= ymd_hm(201809160955) & rga$datetime < ymd_hm(201809162155))
ebullitionC = which(rga$datetime >= ymd_hm(201809170155) & rga$datetime < ymd_hm(201809170740))
ebullitionC[ebullitionC %in% wp] = NA
ebullitionC <- na.omit(ebullitionC)
ebullition = c(ebullitionA,ebullitionB,ebullitionC)
quiet1 <- 4506:4893
quiet2 <- 6110:6686
quiet3 <- 12101:12867
quiet <- c(quiet1,quiet2,quiet3)

# calculate moving average ####
dat = rga # 1 minute simple moving average
for(i in 2:644){
  dat[,i] <- stats::filter(dat[,i],rep(1/7,7),sides = 2) 
}


#o2 and h2o detrend ####
#cut out everything before accidental on/off period
use_dat <- tibble(datetime=as.numeric(rga$datetime-rga$datetime[1068]),
                  o2=rowMeans(rga[,312:314]),
                  h2o=rowMeans(rga[,172:174])) %>%
  slice(-(1068:1268)) %>%
  slice(-(1:1067)) 

#find initial guesses for 2nd order exponential decay model
plot(o2~datetime,data=use_dat,pch=".")
ig1 = c(yf=5e-9,y0=8e-9,y1=8e-9,log_alpha1=-9.5,log_alpha2=-10)
mod_o2 = ig1["yf"]+(ig1["y0"]-ig1["yf"])*exp(-exp(ig1["log_alpha1"])*use_dat$datetime)+(ig1["y1"]-ig1["yf"])*exp(-exp(ig1["log_alpha1"])*use_dat$datetime)
lines(use_dat$datetime,mod_o2,col="green")

plot(h2o~datetime,data=use_dat,pch=".")
ig2 = c(yf=3e-7,y0=7e-7,y1=7e-7,log_alpha1=-9.5,log_alpha2=-10)
mod_h2o = ig2[1]+(ig2[2]-ig2[1])*exp(-exp(ig2[4])*use_dat$datetime)+(ig2[3]-ig2[1])*exp(-exp(ig2[5])*use_dat$datetime)
lines(use_dat$datetime,mod_h2o,col="green")

dt_o2 <- nls(o2 ~ yf+(y0-yf)*exp(-exp(log_alpha1)*datetime)+(y1-yf)*exp(-exp(log_alpha2)*datetime),
             data=use_dat,
             start=as.list(ig1))

dt_h2o <- nls(h2o ~ yf+(y0-yf)*exp(-exp(log_alpha1)*datetime)+(y1-yf)*exp(-exp(log_alpha2)*datetime),
              data=use_dat,
              start=as.list(ig2)) 

pars1 <- dt_o2$m$getPars()
pars2 <- dt_h2o$m$getPars()

use_dat2 <- tibble(datetime=as.numeric(conc_dat0$date_time-conc_dat0$date_time[1068]),
                   o2=conc_dat0$oxygen_uM) %>%
  slice(-(1068:1268)) %>%
  slice(-(1:1067)) 

plot(o2~datetime,data=use_dat2,pch=".")
ig3 = c(yf=1,y0=100,y1=100,log_alpha1=-8.6,log_alpha2=-10.5)
mod_o2 = ig3["yf"]+(ig3["y0"]-ig3["yf"])*exp(-exp(ig3["log_alpha1"])*use_dat2$datetime)+(ig3["y1"]-ig3["yf"])*exp(-exp(ig3["log_alpha1"])*use_dat2$datetime)
lines(use_dat2$datetime,mod_o2,col="green")

dt_o2conc <- nls(o2 ~ yf+(y0-yf)*exp(-exp(log_alpha1)*datetime)+(y1-yf)*exp(-exp(log_alpha2)*datetime),
                 data=use_dat2,
                 start=as.list(ig3))

pars3 <- dt_o2conc$m$getPars()

#check fits
par(mfrow=c(3,1))
plot(o2~datetime,data=use_dat,pch=".",ylim=c(0,2e-8))
mod_o2 = pars1[1]+(pars1[2]-pars1[1])*exp(-exp(pars1[4])*use_dat$datetime)+(pars1[3]-pars1[1])*exp(-exp(pars1[5])*use_dat$datetime)
lines(use_dat$datetime,mod_o2,col="red")
plot(h2o~datetime,data=use_dat,pch=".",ylim=c(0,1e-6))
mod_h2o = pars2[1]+(pars2[2]-pars2[1])*exp(-exp(pars2[4])*use_dat$datetime)+(pars2[3]-pars2[1])*exp(-exp(pars2[5])*use_dat$datetime)
lines(use_dat$datetime,mod_h2o,col="red")
plot(o2~datetime,data=use_dat2,pch=".")
mod_o2 = pars3[1]+(pars3[2]-pars3[1])*exp(-exp(pars3[4])*use_dat2$datetime)+(pars3[3]-pars3[1])*exp(-exp(pars3[5])*use_dat2$datetime)
lines(use_dat2$datetime,mod_o2,col=2)

# objects with only the good parts ####
dat1 <- dat 
dat1[-active_onsite,-1] <- NA
dat1[background,] <- dat[background,]

conc_dat<-conc_dat0
conc_dat[is.na(dat1$total),2:11] <- NA

ox1[ox1$datetime < power_blip+600 | (ox1$datetime >= ROV_pu[2] & ox1$datetime<ROV_clear),-1] <- NA
ox2[ox2$datetime < power_blip+600 | (ox2$datetime >= ROV_pu[2] & ox2$datetime<ROV_clear),-1] <- NA

#take out weird peaks
dat2 <- dat1 
dat2[wp,-1] <- NA

#dat1_old=dat1

#dat1<-dat2
#dat1<-dat1_old

# pull out relevent vectors for gases  ####
datetime <- dat1$datetime#[-1:-100]
h2 <- rowMeans(dat1[,12:14])
h2o <- rowMeans(dat1[,172:174])
n2 <- rowMeans(dat1[,272:274])
o2 <- rowMeans(dat1[,312:314])
h2s <- rowMeans(dat1[,332:334])
ar <- rowMeans(dat1[,392:394])
co2 <- rowMeans(dat1[,432:434])
ch4 <- rowMeans(dat1[,142:144])
c2h6 <- rowMeans(dat1[,282:284])
noise <- dat1$MassSpectrum_641
noises <- rowMeans(dat1[,c(463:642)]) #c(23:103,342:382,443:642)
noises[noises<1e-10]=1e-10
total <-dat1$total

#create detrended vectors
xs <- as.numeric(dat1$datetime-dat1$datetime[1068])
o2_dt  <-  o2-(pars1["y0"]-pars1["yf"])*exp(-exp(pars1["log_alpha1"])*xs) - (pars1["y1"]-pars1["yf"])*exp(-exp(pars1["log_alpha2"])*xs)
h2o_dt <- h2o-(pars2["y0"]-pars2["yf"])*exp(-exp(pars2["log_alpha1"])*xs) - (pars2["y1"]-pars2["yf"])*exp(-exp(pars2["log_alpha2"])*xs)
o2_conc_dt <- conc_dat$oxygen_uM-(pars3["y0"]-pars3["yf"])*exp(-exp(pars3["log_alpha1"])*xs) - (pars3["y1"]-pars3["yf"])*exp(-exp(pars3["log_alpha2"])*xs)

#figure area
start = ymd_hm(201809151630)
end = ymd_hm(201809170800)
Time = ymd_hm(c(201809151200,201809151600,201809152000,201809160000,201809160400,201809160800,201809161200,201809161600,201809162000,201809170000,201809170400,201809170800))
time_names = c("12:00","16:00","20:00", "16-Sep","4:00","8:00","12:00","16:00","20:00","17-Sep","4:00","8:00")


#pub figure current selection! Use this! ####
norm = h2o_dt/mean(h2o_dt,na.rm = TRUE)*ar/mean(ar,na.rm = TRUE)
norm=h2o_dt/mean(h2o_dt,na.rm = TRUE)

setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/calibrated_gases.eps",height=9,width=7.5)
par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(2,.2,0),tcl=.25,las=1)

par(fig=c(.12,.9,.875,.99))
plot(datetime,conc_dat$methane_uM/1000/norm,type="l",col=2,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-10,120),
     frame.plot=TRUE,ann=FALSE,axes=FALSE)
axis(2,seq(-50,200,by=50))
mtext(expression(CH[4],(mM)),2,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
#axis(4,seq(-50,150,by=50)*1.49,seq(-50,150,by=50))
#mtext(expression(CH[4]~("%"~saturation)),4,1,las=0,cex=.8)

par(fig=c(.12,.9,.76,.875),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,10e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(expression(C[2]*H[6],(nTorr)),4,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(4,0:2*5e-9,0:2*.5)

par(fig=c(.12,.9,.645,.76),new=TRUE)
plot(datetime,conc_dat$hydrogen_uM/1000/norm,type="l",col=1,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-.10,.50),
     frame.plot=TRUE,ann=FALSE,axes=FALSE)
axis(2,seq(-.3,1.2,by=.3))
mtext(expression(H[2],(mM)),2,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
#axis(4,seq(-.5,1.5,by=.5)*.65,seq(-.5,1.5,by=.5))
#mtext(expression(H[2]~("%"~saturation)),4,1,las=0,cex=.8)

par(fig=c(.12,.9,.53,.645),new=TRUE)
plot(datetime,o2_dt/norm*1e9+5,type="l",col=6,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(5,16),
     frame.plot=TRUE,ann=FALSE,axes=FALSE)
axis(2,seq(-10,20,by=5))
mtext(expression(O[2],(µM)),2,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
lines(ox1$datetime,ox1$O2_uM_corr)
#lines(datetime,o2_conc_dt/norm/10+10,type="l")
mtext(expression(O[2],(nTorr)),4,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(4,0:3*5,0:3*5-5)

par(fig=c(.12,.9,.415,.53),new=TRUE)
plot(datetime,conc_dat$co2_uM/1000/norm,type="l",col=8,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(.5,4.9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(expression(CO[2],(mM)),2,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:6)

par(fig=c(.12,.9,.3,.415),new=TRUE)
plot(datetime,n2/norm,type="l",col=3,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(50e-9,450e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(expression(N[2],(nTorr)),4,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(4,(0:5)*100e-9,0:5)

par(fig=c(.12,.9,.185,.3),new=TRUE)
plot(datetime,ar/(h2o_dt/mean(h2o_dt,na.rm = TRUE)),type="l",col=5,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,22e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(expression(Ar,(nTorr)),4,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(4,0:2*10e-9,0:2*10)

par(fig=c(.12,.9,.07,.185),new=TRUE)
plot(datetime,h2o_dt/1,type="l",col=4,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(50e-9,1000e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(expression(H[2]*O,(nTorr)),4,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(4,0:3*200e-9,0:3*200)
axis(1,Time,time_names,cex.axis=.8)


dev.off()











# Pub Figure: seperate time series option 0: not normalized ####
norm = 1

setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/timeseries0.eps",height=9,width=7.5)
par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(2,.2,0),tcl=.25,las=1)
par(fig=c(.12,.99,.875,.99))
plot(datetime,ch4/norm,type="l",col=2,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(-.1e-6,3e-6)*3e6,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CH[4]),expression(uTorr)),2,3,padj=c(-.5,.7),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*1e-6,0:3)

par(fig=c(.12,.99,.76,.875),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,11e-9)*3e6,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(C[2]*H[6]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*5e-9,0:2*.5)

par(fig=c(.12,.99,.645,.76),new=TRUE)
plot(datetime,h2/norm,type="l",col=1,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,110e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(H[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*50e-9,0:2*50)

par(fig=c(.12,.99,.53,.645),new=TRUE)
plot(datetime,o2_dt/norm,type="l",col=6,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,1.1e-8),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(O[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*5e-9,0:3*5)

par(fig=c(.12,.99,.415,.53),new=TRUE)
plot(datetime,co2/norm,type="l",col=8,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,33e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CO[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*10e-9,0:3*10)

par(fig=c(.12,.99,.3,.415),new=TRUE)
plot(datetime,n2/norm,type="l",col=3,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,500e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(N[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*200e-9,0:2*200)

par(fig=c(.12,.99,.185,.3),new=TRUE)
plot(datetime,ar/norm,type="l",col=5,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,22e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(Ar),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*10e-9,0:2*10)

par(fig=c(.12,.99,.07,.185),new=TRUE)
plot(datetime,h2o_dt/norm,type="l",col=4,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,22e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(H[2]*O),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*10e-9,0:2*10)
axis(1,Time,time_names,cex.axis=.8)

dev.off()

# Pub Figure: seperate time series option 1: water normalized ####
norm=h2o_dt

setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/timeseries1.eps",height=9,width=7.5)
par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(2,.2,0),tcl=.25,las=1)
par(fig=c(.12,.99,.86,.99))
plot(datetime,ch4/norm,type="l",col=2,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-.1e-6,3e-6)*3e6,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CH[4]),expression(uTorr)),2,3,padj=c(-.5,.7),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*1e-6,0:3)

par(fig=c(.15,.39,.93,.98),new=TRUE)
plot(datetime,ch4/norm,type="l",col=2,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(-1e-6,20e-6),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
axis(4)#axis(4,0:2*10e-6,0:2*10)

par(fig=c(.12,.99,.73,.86),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,11e-9)*3e6,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(C[2]*H[6]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*5e-9,0:2*.5)

par(fig=c(.15,.39,.79,.84),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(-10*2e-9,220e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
axis(4)#axis(4,0:2*100e-9,0:2*100)

par(fig=c(.12,.99,.6,.73),new=TRUE)
plot(datetime,h2/norm,type="l",col=1,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,110e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(H[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*50e-9,0:2*50)

par(fig=c(.12,.99,.47,.6),new=TRUE)
plot(datetime,o2_dt/norm,type="l",col=6,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,1.1e-8),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(O[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*5e-9,0:3*5)

par(fig=c(.12,.99,.34,.47),new=TRUE)
plot(datetime,co2/norm,type="l",col=8,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,33e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CO[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*10e-9,0:3*10)

par(fig=c(.12,.99,.21,.34),new=TRUE)
plot(datetime,n2/norm,type="l",col=3,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,500e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(N[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*200e-9,0:2*200)

par(fig=c(.12,.99,.08,.21),new=TRUE)
plot(datetime,ar/norm,type="l",col=5,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,22e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(Ar),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*10e-9,0:2*10)
axis(1,Time,time_names,cex.axis=.8)

dev.off()

# Pub Figure: seperate time series option 2: water-anomaly normalized ####
norm=h2o_dt/mean(h2o_dt,na.rm = TRUE)
#norm=h2o/mean(h2o,na.rm = TRUE)

setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/timeseries2.eps",height=9,width=7.5)
par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(2,.2,0),tcl=.25,las=1)
par(fig=c(.12,.99,.86,.99))
plot(datetime,ch4/norm,type="l",col=2,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-.1e-6,3e-6),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CH[4]),expression(uTorr)),2,3,padj=c(-.5,.7),adj=.5,cex=c(1,.8))
axis(2,0:3*1e-6,0:3)

par(fig=c(.15,.39,.93,.98),new=TRUE)
plot(datetime,ch4/norm,type="l",col=2,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-1e-6,20e-6),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
axis(4,0:2*10e-6,0:2*10)

par(fig=c(.12,.99,.73,.86),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,11e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(C[2]*H[6]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:2*5e-9,0:2*.5)

par(fig=c(.15,.39,.79,.84),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-10*2e-9,220e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
axis(4,0:2*100e-9,0:2*100)

par(fig=c(.12,.99,.6,.73),new=TRUE)
plot(datetime,h2/norm,type="l",col=1,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,110e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(H[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:2*50e-9,0:2*50)

par(fig=c(.12,.99,.47,.6),new=TRUE)
plot(datetime,o2_dt/norm,type="l",col=6,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,1.1e-8),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(O[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:3*5e-9,0:3*5)

par(fig=c(.12,.99,.34,.47),new=TRUE)
plot(datetime,co2/norm,type="l",col=8,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,33e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CO[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:3*10e-9,0:3*10)

par(fig=c(.12,.99,.21,.34),new=TRUE)
plot(datetime,n2/norm,type="l",col=3,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,500e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(N[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:2*200e-9,0:2*200)

par(fig=c(.12,.99,.08,.21),new=TRUE)
plot(datetime,ar/norm,type="l",col=5,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,22e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(Ar),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:2*10e-9,0:2*10)

axis(1,Time,time_names,cex.axis=.8)

dev.off()


# Pub Figure: seperate time series option 3: argon anomaly normalized ####
norm = ar/mean(ar,na.rm = TRUE)

setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/timeseries3.eps",height=9,width=7.5)
par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(2,.2,0),tcl=.25,las=1)
par(fig=c(.12,.99,.86,.99))
plot(datetime,ch4/norm,type="l",col=2,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-.1e-6,3000e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CH[4]),expression(nTorr)),2,3,padj=c(-.5,.7),adj=.5,cex=c(1,.8))
axis(2,0:3*1000e-9,c(0,paste0(1:3,"k")))

par(fig=c(.15,.39,.93,.98),new=TRUE)
plot(datetime,ch4/norm,type="l",col=2,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-1e-6,20000e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
axis(4,0:2*10000e-9,c(0,paste0(1:2*10,"k")))

par(fig=c(.12,.99,.73,.86),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,11e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(C[2]*H[6]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:2*5e-9,0:2*.5)

par(fig=c(.15,.39,.8,.85),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-10*2e-9,220e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
axis(4,0:2*100e-9,0:2*100)

par(fig=c(.12,.99,.6,.73),new=TRUE)
plot(datetime,h2/norm,type="l",col=1,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,110e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(H[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:2*50e-9,0:2*50)

par(fig=c(.12,.99,.47,.6),new=TRUE)
plot(datetime,o2_dt/norm,type="l",col=6,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(1,9)*1e-9,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(O[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:4*2e-9,0:4*2)

par(fig=c(.12,.99,.34,.47),new=TRUE)
plot(datetime,h2o_dt/norm,type="l",col=4,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(50,700)*1e-9,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(H[2]*O),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:4*200e-9,0:4*200)

par(fig=c(.12,.99,.21,.34),new=TRUE)
plot(datetime,n2/norm,type="l",col=3,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(170,320)*1e-9,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(N[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:3*100e-9,0:3*100)

par(fig=c(.12,.99,.08,.21),new=TRUE)
plot(datetime,co2/norm,type="l",col=8,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(5,28)*1e-9,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CO[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:3*10e-9,0:3*10)

axis(1,Time,time_names,cex.axis=.8)

dev.off()

# Pub Figure: seperate time series option 4: water AND argon normalized ####
norm = h2o_dt/mean(h2o_dt,na.rm = TRUE)*ar/mean(ar,na.rm = TRUE)

setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/timeseries4.eps",height=9,width=7.5)
par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(2,.2,0),tcl=.25,las=1)
par(fig=c(.12,.99,.83,.99))
plot(datetime,ch4/norm,type="l",col=2,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-.1e-6,3000e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CH[4]),expression(nTorr)),2,3,padj=c(-.5,.7),adj=.5,cex=c(1,.8))
axis(2,0:3*1000e-9,c(0,paste0(1:3,"k")))

par(fig=c(.12,.99,.67,.83),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,6e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(C[2]*H[6]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:2*5e-9,0:2*.5)

par(fig=c(.12,.99,.51,.67),new=TRUE)
plot(datetime,h2/norm,type="l",col=1,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,110e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(H[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:2*50e-9,0:2*50)

par(fig=c(.12,.99,.35,.51),new=TRUE)
plot(datetime,o2_dt/norm,type="l",col=6,xaxs="i",#log="y",
     xlim=c(start,end),ylim=c(1,9)*1e-9,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(O[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:4*2e-9,0:4*2)

par(fig=c(.12,.99,.19,.35),new=TRUE)
plot(datetime,n2/norm,type="l",col=3,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(170,320)*1e-9,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(N[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:3*100e-9,0:3*100)

par(fig=c(.12,.99,.03,.19),new=TRUE)
plot(datetime,co2/norm,type="l",col=8,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(5,28)*1e-9,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CO[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:3*10e-9,0:3*10)

axis(1,Time,time_names,cex.axis=.8)

dev.off()

# Pub Figure: seperate time series option 5: total-P normalized ####
norm = total

setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/timeseries5.eps",height=9,width=7.5)
par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(2,.2,0),tcl=.25,las=1)
par(fig=c(.12,.99,.875,.99))
plot(datetime,ch4/norm,type="l",col=2,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(-.1e-6,3e-6)*3e6,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CH[4]),expression(uTorr)),2,3,padj=c(-.5,.7),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*1e-6,0:3)

par(fig=c(.12,.99,.76,.875),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,11e-9)*3e6,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(C[2]*H[6]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*5e-9,0:2*.5)

par(fig=c(.12,.99,.645,.76),new=TRUE)
plot(datetime,h2/norm,type="l",col=1,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,110e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(H[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*50e-9,0:2*50)

par(fig=c(.12,.99,.53,.645),new=TRUE)
plot(datetime,o2_dt/norm,type="l",col=6,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,1.1e-8),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(O[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*5e-9,0:3*5)

par(fig=c(.12,.99,.415,.53),new=TRUE)
plot(datetime,co2/norm,type="l",col=8,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,33e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CO[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*10e-9,0:3*10)

par(fig=c(.12,.99,.3,.415),new=TRUE)
plot(datetime,n2/norm,type="l",col=3,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,500e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(N[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*200e-9,0:2*200)

par(fig=c(.12,.99,.185,.3),new=TRUE)
plot(datetime,ar/norm,type="l",col=5,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,22e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(Ar),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*10e-9,0:2*10)

par(fig=c(.12,.99,.07,.185),new=TRUE)
plot(datetime,h2o_dt/norm,type="l",col=4,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,22e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(H[2]*O),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*10e-9,0:2*10)
axis(1,Time,time_names,cex.axis=.8)

dev.off()

# Pub Figure: seperate time series option 6: noise normalized ####
norm = (mean(noise,na.rm=TRUE)-noise)/mean(noise,na.rm=TRUE)
norm = noises

setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/timeseries6.eps",height=9,width=7.5)
par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(2,.2,0),tcl=.25,las=1)
par(fig=c(.12,.99,.875,.99))
plot(datetime,ch4/norm,type="l",col=2,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(-.1e-6,3e-6)*3e6,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CH[4]),expression(uTorr)),2,3,padj=c(-.5,.7),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*1e-6,0:3)

par(fig=c(.12,.99,.76,.875),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,11e-9)*3e6,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(C[2]*H[6]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*5e-9,0:2*.5)

par(fig=c(.12,.99,.645,.76),new=TRUE)
plot(datetime,h2/norm,type="l",col=1,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,110e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(H[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*50e-9,0:2*50)

par(fig=c(.12,.99,.53,.645),new=TRUE)
plot(datetime,o2_dt/norm,type="l",col=6,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,1.1e-8),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(O[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*5e-9,0:3*5)

par(fig=c(.12,.99,.415,.53),new=TRUE)
plot(datetime,co2/norm,type="l",col=8,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,33e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CO[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*10e-9,0:3*10)

par(fig=c(.12,.99,.3,.415),new=TRUE)
plot(datetime,n2/norm,type="l",col=3,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,500e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(N[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*200e-9,0:2*200)

par(fig=c(.12,.99,.185,.3),new=TRUE)
plot(datetime,ar/norm,type="l",col=5,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,22e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(Ar),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*10e-9,0:2*10)

par(fig=c(.12,.99,.07,.185),new=TRUE)
plot(datetime,h2o_dt/norm,type="l",col=4,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,22e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(H[2]*O),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*10e-9,0:2*10)
axis(1,Time,time_names,cex.axis=.8)

dev.off()

# Pub Figure: seperate time series option 7: noise and water normalized ####
norm = noises*h2o_dt

setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/timeseries7.eps",height=9,width=7.5)
par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(2,.2,0),tcl=.25,las=1)

par(fig=c(.12,.99,.86,.99))
plot(datetime,ch4/norm,type="l",col=2,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(-.1e-6,3e-6)*3e6,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CH[4]),expression(uTorr)),2,3,padj=c(-.5,.7),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*1e-6,0:3)

par(fig=c(.15,.39,.93,.98),new=TRUE)
plot(datetime,ch4/norm,type="l",col=2,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(-1e-6,20e-6),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
axis(4)#axis(4,0:2*10e-6,0:2*10)

par(fig=c(.12,.99,.73,.86),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,11e-9)*3e6,
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(C[2]*H[6]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*5e-9,0:2*.5)

par(fig=c(.15,.39,.79,.84),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(-10*2e-9,220e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
axis(4)#axis(4,0:2*100e-9,0:2*100)

par(fig=c(.12,.99,.6,.73),new=TRUE)
plot(datetime,h2/norm,type="l",col=1,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,110e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(H[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*50e-9,0:2*50)

par(fig=c(.12,.99,.47,.6),new=TRUE)
plot(datetime,o2_dt/norm,type="l",col=6,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,1.1e-8),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(O[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*5e-9,0:3*5)

par(fig=c(.12,.99,.34,.47),new=TRUE)
plot(datetime,co2/norm,type="l",col=8,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,33e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(CO[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:3*10e-9,0:3*10)

par(fig=c(.12,.99,.21,.34),new=TRUE)
plot(datetime,n2/norm,type="l",col=3,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,500e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(N[2]),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*200e-9,0:2*200)

par(fig=c(.12,.99,.08,.21),new=TRUE)
plot(datetime,ar/norm,type="l",col=5,xaxs="i",#log="y",
     xlim=c(start,end),#ylim=c(0,22e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(c(expression(Ar),expression(nTorr)),2,3,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2)#axis(2,0:2*10e-9,0:2*10)
axis(1,Time,time_names,cex.axis=.8)

dev.off()

# oxygen optode/rga compare ####

par(mfrow=c(3,1),mar=c(1,1,0,0),#c(0,0,0,0),
    mgp=c(2,.2,0),tcl=.25,las=3,cex.axis=.8)
#par(fig=c(.12,.88,.17,.99))
#plot(datetime,o2_conc_dt,ylim=c(0,50),type="l",col=6)
plot(datetime,o2_conc_dt,col=2,
     type="l",
     #xlim=c(ymd_hm(201809151900),ymd_hm(201809170430)),
     #ylim=c(5,20),
     #xaxs="i",yaxs="i",
     ann=FALSE,axes=FALSE,frame.plot = TRUE)
lines(ox1$datetime,ox1$O2_uM_corr)
lines(datetime,o2_dt/ar*mean(ar,na.rm = TRUE)*1e9+5,col=6,)
axis(2,c(0:8*2))
mtext(expression(DO~(uM)),2,1,padj=-.5)
axis(4,c(0:8*2),c(0:8*2)-5)
mtext(c("nTorr at 32 ± 0.1 m/z"),4,2,padj=c(-.5))
#axis(1,at=Time,labels=time_names)
legend("topright",c("Optode","Ar-normalized ISMS"),col=c(1,6),lwd=2,bty="n")

plot(datetime,o2_dt/h2o_dt*mean(h2o_dt,na.rm = TRUE)*1e9+5,col=6,type="l",
     #xlim=c(ymd_hm(201809151900),ymd_hm(201809170430)),
     ylim=c(5,20),
     #xaxs="i",yaxs="i",
     ann=FALSE,axes=FALSE,frame.plot = TRUE)
lines(ox1$datetime,ox1$O2_uM_corr)
legend("topright",c("Optode","water-normalized ISMS"),col=c(1,6),lwd=2,bty="n")
plot(datetime,o2_dt/h2o_dt*mean(h2o_dt,na.rm = TRUE)/ar*mean(ar,na.rm = TRUE)*1e9+5,col=6,type="l",
     #xlim=c(ymd_hm(201809151900),ymd_hm(201809170430)),
     ylim=c(5,20),
     #xaxs="i",yaxs="i",
     ann=FALSE,axes=FALSE,frame.plot = TRUE)
lines(ox1$datetime,ox1$O2_uM_corr)
legend("topright",c("Optode","Water&Ar-normalized ISMS"),col=c(1,6),lwd=2,bty="n")


# Pub Plot(s): with concentration ####
norm = h2o_dt/mean(h2o_dt,na.rm = TRUE)*ar/mean(ar,na.rm = TRUE)


setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/calibrated_gases.eps",height=9,width=7.5)
par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(2,.2,0),tcl=.25,las=1)

par(fig=c(.12,.9,.875,.99))
plot(datetime,conc_dat$methane_uM/1000/norm,type="l",col=2,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-10,120),
     frame.plot=TRUE,ann=FALSE,axes=FALSE)
axis(2,seq(-50,200,by=50))
mtext(expression(CH[4],(mM)),2,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
#axis(4,seq(-50,150,by=50)*1.49,seq(-50,150,by=50))
#mtext(expression(CH[4]~("%"~saturation)),4,1,las=0,cex=.8)

par(fig=c(.12,.9,.76,.875),new=TRUE)
plot(datetime,c2h6/norm,type="l",col="purple",
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,10e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(expression(C[2]*H[6],(nTorr)),4,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(4,0:2*5e-9,0:2*.5)

par(fig=c(.12,.9,.645,.76),new=TRUE)
plot(datetime,conc_dat$hydrogen_uM/1000/norm,type="l",col=1,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-.010,.800),
     frame.plot=TRUE,ann=FALSE,axes=FALSE)
axis(2,seq(-.3,1.2,by=.3))
mtext(expression(H[2],(mM)),2,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
#axis(4,seq(-.5,1.5,by=.5)*.65,seq(-.5,1.5,by=.5))
#mtext(expression(H[2]~("%"~saturation)),4,1,las=0,cex=.8)

par(fig=c(.12,.9,.53,.645),new=TRUE)
plot(datetime,o2_dt/norm*1e9+5,type="l",col=6,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(5,16),
     frame.plot=TRUE,ann=FALSE,axes=FALSE)
axis(2,seq(-10,20,by=5))
mtext(expression(O[2],(µM)),2,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
lines(ox1$datetime,ox1$O2_uM_corr)
#lines(datetime,o2_conc_dt/norm/10+10,type="l")
mtext(expression(O[2],(nTorr)),4,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(4,0:3*5,0:3*5-5)

par(fig=c(.12,.9,.415,.53),new=TRUE)
plot(datetime,conc_dat$co2_uM/1000/norm,type="l",col=8,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(.5,4.9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(expression(CO[2],(mM)),2,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(2,0:6)

par(fig=c(.12,.9,.3,.415),new=TRUE)
plot(datetime,n2/norm,type="l",col=3,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(100e-9,350e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(expression(N[2],(nTorr)),4,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(4,(0:4)*100e-9,0:4)

par(fig=c(.12,.9,.185,.3),new=TRUE)
plot(datetime,ar/(h2o_dt/mean(h2o_dt,na.rm = TRUE)),type="l",col=5,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(0,22e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(expression(Ar,(nTorr)),4,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(4,0:2*10e-9,0:2*10)

par(fig=c(.12,.9,.07,.185),new=TRUE)
plot(datetime,h2o_dt/(ar/mean(ar,na.rm = TRUE)),type="l",col=4,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(50e-9,590e-9),
     ann=FALSE,axes=FALSE,frame.plot=TRUE)
mtext(expression(H[2]*O,(nTorr)),4,2,padj=c(-.5,1.5),adj=.5,cex=c(1,.8))
axis(4,0:3*200e-9,0:3*200)
axis(1,Time,time_names,cex.axis=.8)


dev.off()



par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(2,.2,0),tcl=.25,las=1,cex.axis=.8)
par(fig=c(.12,.9,.68,.99))
plot(conc_dat$date_time,conc_dat$methane_uM/1000,type="l",col=2,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-10,150),
     frame.plot=TRUE,ann=FALSE,axes=FALSE)
axis(2,seq(-50,200,by=50))
mtext(expression(CH[4]~(mM)),2,1.5,las=0,cex=.8)
axis(4,seq(-50,150,by=50)*1.49,seq(-50,150,by=50))
mtext(expression(CH[4]~("%"~saturation)),4,1,las=0,cex=.8)

par(fig=c(.12,.9,.37,.68),new=TRUE)
plot(conc_dat$date_time,conc_dat$hydrogen_uM/1000,type="l",col=1,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-.010,.800),
     frame.plot=TRUE,ann=FALSE,axes=FALSE)
axis(2,seq(-.3,1.2,by=.3))
mtext(expression(H[2]~(mM)),2,1.5,las=0,cex=.8)
axis(4,seq(-.5,1.5,by=.5)*.65,seq(-.5,1.5,by=.5))
mtext(expression(H[2]~("%"~saturation)),4,1,las=0,cex=.8)

par(fig=c(.12,.9,.06,.37),new=TRUE)
plot(conc_dat$date_time,o2_conc_dt,type="l",col=6,
     xaxs="i",yaxs="i",#log="y",
     xlim=c(start,end),ylim=c(-20,100),
     frame.plot=TRUE,ann=FALSE,axes=FALSE)
axis(2,seq(-50,200,by=50))
mtext(expression(O[2]~(uM)),2,1.5,las=0,cex=.8)

axis(1,Time,time_names,cex.axis=.8)

dev.off()

