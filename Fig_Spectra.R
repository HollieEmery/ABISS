library("tidyverse")

#data import ####
rga <- read_csv("rga.csv")

# calculate moving average ####
dat = rga # 1 minute simple moving average
for(i in 2:644){
  dat[,i] <- stats::filter(dat[,i],rep(1/7,7),sides = 2) 
}

# indices of phases of deployment 
background = which(rga$datetime >= ymd_hm(201809151735)+20*60 & rga$datetime < ymd_hm(201809151845)) #lander in the mud while ROV searched for site, starting 10 mins after restart

ebullitionA = which(rga$datetime >= ymd_hm(201809151951) & rga$datetime < ymd_hm(201809160210)) #ymd_hm(201809152125)
ebullitionB = which(rga$datetime >= ymd_hm(201809160955) & rga$datetime < ymd_hm(201809162155))
ebullitionC = which(rga$datetime >= ymd_hm(201809170155) & rga$datetime < ymd_hm(201809170740))
ebullition = c(ebullitionA,ebullitionB,ebullitionC)

wp <- c(13332:13451,13870:13943)
quiet <- c(4506:4893,6110:6686,12101:12867)


#calculate means in phases
bck_means = colMeans(dat[background,2:642],na.rm=TRUE)
qui_means = colMeans(dat[quiet,2:642],na.rm=TRUE)
ebl_means = colMeans(dat[ebullition,2:642],na.rm=TRUE)
wp_means = colMeans(dat[wp,2:642],na.rm=TRUE)

# Pub Figure: Spectra in and out option 0: mulitple lines ####
setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/Spectras.eps",height=4,width=5.5)
par(mfrow=c(1,1),mai=c(0, 0, 0, 0),tcl=.25,mgp=c(3,.1,0),cex.axis=.8)
par(fig=c(.1,.99,.12,.98))
plot(1,1,type="l",ylim=c(1e-11,1e-4),log="y",xlim=c(0,55),xaxs="i",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
axis(2,c((1:9)*10^-12,(1:9)*10^-11,(1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-12}),rep("",17),expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
lines(2:642/10+.8,wp_means,col="red") 
lines(2:642/10+.8,ebl_means,col=3) 
lines(2:642/10+.8,qui_means,col=4) 
lines(2:642/10+.8,bck_means,col="grey60") 

axis(1)
mtext("m/z",1,1,font = 3)
legend("topright",c("Off Seep","Quiescence","Ebullition","Free Gas"),col=c("grey60",4,3,"red"),lty=1,bty="n",cex=.8)

dev.off()

# Pub Figure: Spectra in and out option 1: bold lines ####
setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/Spectra1.eps",height=4,width=5.5)

par(mfrow=c(1,1),mai=c(0, 0, 0, 0),tcl=.25,mgp=c(3,.1,0),cex.axis=.8)
par(fig=c(.1,.99,.54,.98))
plot(1,1,type="l",ylim=c(1e-10,1e-4),log="y",xlim=c(0,55),xaxs="i",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
axis(2,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
abline(h=c(1e-10,1e-8,1e-6,1e-4),col="grey95")
segments(c(14,28),1e-11,y1=bck_means[c(132,272)],col=3,lwd=4) #n2
segments(c(1,17,18),1e-11,y1=bck_means[c(2,162,172)],col=4,lwd=4) #h2o
segments(c(20,40),1e-11,y1=bck_means[c(192,392)],col=5,lwd=4) #ar
segments(c(16,32),1e-11,y1=bck_means[c(151,312)],col=6,lwd=4) #o2
segments(c(12,16,44),1e-11,y1=bck_means[c(112,153,432)],col=8,lty=c(1,2),lwd=4) #co2
legend("topright",expression(N[2],H[2]*O,Ar,O[2],CO[2]),lty=1,lwd=2,col=c(3:6,8),cex=.7,bty="n")
lines(2:642/10+.8,bck_means) 
mtext("a) Background Period",line=-1,adj=0.02,cex=.8)

par(fig=c(.1,.99,.1,.54),new=TRUE)
plot(1,1,type="l",ylim=c(1e-10,1e-4),log="y",xlim=c(0,55),xaxs="i",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
axis(2,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
abline(h=c(1e-10,1e-8,1e-6,1e-4),col="grey95")
segments(2,1e-11,y1=ebl_means[12],lwd=4) #h2
#segments(34,1e-11,y1=ebl_means[332],col=7,lwd=2) #h2s
segments(c(12:17),1e-11,y1=ebl_means[c(112,122,132,142,152,162)],col=2,lwd=4) #methane
segments(c(26:30),1e-11,y1=ebl_means[c(252,262,272,282,292)],col="purple",lwd=4) #ethane
legend("topright",expression(CH[4],C[2]*H[6],H[2]),col=c(2,"purple",1),lty=1,cex=.7,bty="n",lwd=2)
lines(2:642/10+.8,ebl_means)
mtext("b) Peak Ebullition",line=-1,adj=.02,cex=.8)

axis(1)
mtext("m/z",1,1,font = 3)

dev.off()

# Pub Figure: Spectra in and out option 2: grey spectra lines ####

setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/Spectra2.eps",height=4,width=5.5)

par(mfrow=c(1,1),mai=c(0, 0, 0, 0),tcl=.25,mgp=c(3,.1,0),cex.axis=.8)
par(fig=c(.1,.99,.54,.98))
plot(1,1,type="l",ylim=c(1e-10,1e-4),log="y",xlim=c(0,55),xaxs="i",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
axis(2,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
abline(h=c(1e-10,1e-8,1e-6,1e-4),col="grey95")
segments(c(14,28),1e-11,y1=bck_means[c(132,272)],col=3,lty=c(3,1),lwd=3) #n2
segments(c(1,17,18),1e-11,y1=bck_means[c(2,162,172)],col=4,lty=c(3,3,1),lwd=3) #h2o
segments(c(20,40),1e-11,y1=bck_means[c(192,392)],col=5,lty=c(3,1),lwd=3) #ar
segments(c(16,32),1.5e-11,y1=bck_means[c(151,312)],col=6,lty=c(3,1),lwd=3) #o2
segments(c(12,16,44),1e-11,y1=bck_means[c(112,153,432)],col=8,lty=c(3,3,1),lwd=3) #co2
legend("topright",expression(N[2],H[2]*O,Ar,O[2],CO[2]),lty=1,lwd=2,col=c(3:6,8),cex=.7,bty="n")
lines(2:642/10+.8,bck_means,col="grey80") 
mtext("a) Background Period",line=-1,adj=0.02,cex=.8)

par(fig=c(.1,.99,.1,.54),new=TRUE)
plot(1,1,type="l",ylim=c(1e-10,1e-4),log="y",xlim=c(0,55),xaxs="i",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
axis(2,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
abline(h=c(1e-10,1e-8,1e-6,1e-4),col="grey95")
segments(2,1e-11,y1=ebl_means[12],lwd=3) #h2
#segments(34,1e-11,y1=ebl_means[332],col=7,lwd=2) #h2s
segments(c(12:17),1e-11,y1=ebl_means[c(112,122,132,142,152,162)],col=2,lty=c(3,3,3,2,1,3),lwd=3) #methane
segments(c(26:30),1e-11,y1=ebl_means[c(252,262,272,282,292)],col="purple",lty=c(3,3,1,2,3),lwd=3) #ethane
legend("topright",expression(CH[4],C[2]*H[6],H[2]),col=c(2,"purple",1),lty=1,cex=.7,bty="n",lwd=2)
lines(2:642/10+.8,ebl_means,col="grey80")
mtext("b) Peak Ebullition",line=-1,adj=.02,cex=.8)

axis(1)
mtext("m/z",1,1,font = 3)

dev.off()

# Pub Figure: Spectra in and out option 3: stars ####
setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/Spectra3.eps",height=4,width=5.5)

par(mfrow=c(1,1),mai=c(0, 0, 0, 0),tcl=.25,mgp=c(3,.1,0),cex.axis=.8)
par(fig=c(.1,.99,.54,.98))
plot(1,1,type="l",ylim=c(1e-10,1e-4),log="y",xlim=c(0,55),xaxs="i",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
axis(2,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
abline(h=c(1e-10,1e-8,1e-6,1e-4),col="grey95")
lines(2:642/10+.8,bck_means) 
points(c(14,28),bck_means[c(132,272)],col=3,pch=8) #n2
points(c(1,17,18),bck_means[c(2,162,172)],col=4,pch=8) #h2o
points(c(20,40),bck_means[c(192,392)],col=5,pch=8) #ar
points(c(16,32),bck_means[c(151,312)],col=6,pch=8) #o2
points(c(12,16,44),bck_means[c(112,153,432)],col=8,pch=8) #co2
legend("topright",expression(N[2],H[2]*O,Ar,O[2],CO[2]),pch=8,col=c(3:6,8),cex=.7,bty="n")
mtext("a) Background Period",line=-1,adj=0.02,cex=.8)

par(fig=c(.1,.99,.1,.54),new=TRUE)
plot(1,1,type="l",ylim=c(1e-10,1e-4),log="y",xlim=c(0,55),xaxs="i",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
axis(2,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
abline(h=c(1e-10,1e-8,1e-6,1e-4),col="grey95")
lines(2:642/10+.8,ebl_means)
points(2,ebl_means[12],pch=8) #h2
points(c(12:17),ebl_means[c(112,122,132,142,152,162)],col=2,pch=8) #methane
points(c(26:30),ebl_means[c(252,262,272,282,292)],col="purple",pch=8) #ethane
legend("topright",expression(CH[4],C[2]*H[6],H[2]),col=c(2,"purple",1),pch=8,cex=.7,bty="n")
mtext("b) Peak Ebullition",line=-1,adj=.02,cex=.8)

axis(1)
mtext("m/z",1,1,font = 3)

dev.off()


# Pub Figure: Spectra in and out option 4: black and red ####
setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/Spectra4.eps",height=3,width=5.5)

par(mfrow=c(1,1),mai=c(0, 0, 0, 0),tcl=.25,mgp=c(3,.1,0),cex.axis=.8)
par(fig=c(.1,.9,.15,.98))
plot(1,1,type="l",ylim=c(1e-10,1e-3),xlim=c(0,55),xaxs="i",log="y",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
abline(h=c(1e-10,1e-8,1e-6,1e-4),col="grey95")
axis(2,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4)*1.25+1.1e-10,
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
axis(4,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",4,1.5,cex=.8)
lines(2:642/10+.8,bck_means*1.25+1.1e-10)
lines(2:642/10+.8,ebl_means,col="red")
text(c(44,40,32,28,18),
     c(1e-7,5e-8,3e-8,2e-6,2e-6),
     c(expression(CO[2]),"Ar",expression(O[2]),expression(N[2]),expression(H[2]*O)))
text(c(2,14.5,28),c(8e-7,3e-4,3e-5),c(expression(H[2]),expression(CH[4]),expression(C[2]*H[6])),col="red")
arrows(c(12,26),c(1e-4,1e-5),x1=c(17,30),length = .06,col="red",code = 3)
axis(1)
mtext("m/z",1,1,font = 3)
legend("topright",c("Background","Peak Ebullition"),col=c(1,"red"),lty=1,bty="n",cex=.8)

dev.off()

# Pub Figure: Spectra in and out option 5: color labels ####
setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/Spectra5.eps",height=3,width=5.5)

par(mfrow=c(1,1),mai=c(0, 0, 0, 0),tcl=.25,mgp=c(3,.1,0),cex.axis=.8)
par(fig=c(.1,.99,.15,.98))
plot(1,1,type="l",ylim=c(1e-10,1e-3),log="y",xlim=c(0,55),xaxs="i",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
axis(2,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
abline(h=c(1e-10,1e-8,1e-6,1e-4),col="grey95")
lines(2:642/10+.8,bck_means,lty=3) 
lines(2:642/10+.8,ebl_means) 
text(c(44,40,32,28,20,18,16,16,14,12,2),
     c(8e-8,5e-8,2e-8,1.5e-6,6e-9,1.5e-6,3e-8,6e-8,3e-8,2e-9,5e-7),
     c(expression(CO[2]),"Ar",expression(O[2]),expression(N[2]),"Ar",expression(H[2]*O),"O","O","N","C",expression(H[2])),
     col=c(8,5,6,3,5,4,8,6,3,8,1),
     cex=.8)
text(c(14.5,28),
     c(1.5e-4,5e-6),
     c(expression(CH[4]),expression(C[2]*H[6])),
     col=c(2,"purple"),cex=.8)
arrows(c(12,26),c(8e-5,3e-6),x1=c(17,30),length = .06,col=c(2,"purple"),code = 3)
axis(1)
mtext("m/z",1,1,font = 3)
legend("topright",c("Background","Peak Ebullition"),lty=c(3,1),bty="n",cex=.8)

dev.off()

# Pub Figure: Spectra in and out option 6: stars and dots ####
setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/Spectra6.eps",height=3,width=5.5)

par(mfrow=c(1,1),mai=c(0, 0, 0, 0),tcl=.25,mgp=c(3,.1,0),cex.axis=.8)
par(fig=c(.1,.99,.15,.98))
plot(1,1,type="l",ylim=c(1e-10,1e-4),log="y",xlim=c(0,55),xaxs="i",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
axis(2,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
abline(h=c(1e-10,1e-8,1e-6,1e-4),col="grey95")
lines(2:642/10+.8,bck_means,col="grey60") 
lines(2:642/10+.8,ebl_means) 
points(2,ebl_means[12],pch=8) #h2
points(c(12:17),ebl_means[c(112,122,132,142,152,162)],col=2,pch=8) #methane
points(c(26:30),ebl_means[c(252,262,272,282,292)],col="purple",pch=8) #ethane
points(c(20,40),ebl_means[c(192,392)],col=5,pch=8) #ar
points(c(44),ebl_means[c(432)],col=8,pch=8) #co2
points(c(32),ebl_means[c(312)],col=6,pch=8) #o2
points(2,bck_means[12],pch=19) #h2
points(c(14,28),bck_means[c(132,272)],col=3,pch=19) #n2
points(c(17,18),bck_means[c(162,172)],col=4,pch=19) #h2o
points(c(20,40),bck_means[c(192,392)],col=5,pch=19) #ar
points(c(16,32),bck_means[c(151,312)],col=6,pch=19) #o2
points(c(12,16,44),bck_means[c(112,153,432)]+c(0,2e-9,0),col=8,pch=19) #co2

axis(1)
mtext("m/z",1,1,font = 3)
legend(49,1.5e-4,expression(N[2],H[2]*O,Ar,O[2],CO[2],CH[4],C[2]*H[6],H[2]),pch=8,col=c(3:6,8,2,"purple",1),cex=.8,bty="n")
legend(32,1.5e-4,c("Background","Peak Ebullition"),lty=c(1,1),col=c("grey60",1),bty="n",cex=.8)

dev.off()

# Pub Figure: Spectra in and out option 7: three plots ####
setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/Spectra7.eps",height=6,width=5.5)

par(mfrow=c(1,1),mai=c(0, 0, 0, 0),tcl=.25,mgp=c(3,.1,0),cex.axis=.8)
par(fig=c(.1,.99,.69,.99))
plot(1,1,type="l",ylim=c(1e-10,1e-4),log="y",xlim=c(0,55),xaxs="i",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
axis(2,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
abline(h=c(1e-10,1e-8,1e-6,1e-4),col="grey95")
lines(2:642/10+.8,bck_means,col="grey60") 
points(c(14,28),bck_means[c(132,272)],col=3,pch=8) #n2
points(c(17,18),bck_means[c(162,172)],col=4,pch=8) #h2o
points(c(20,40),bck_means[c(192,392)],col=5,pch=8) #ar
points(c(16,32),bck_means[c(151,312)],col=6,pch=8) #o2
points(c(12,16,44),bck_means[c(112,153,432)],col=8,pch=8) #co2
legend("topright",expression(N[2],H[2]*O,Ar,O[2],CO[2]),pch=8,col=c(3:6,8),cex=.7,bty="n")
mtext("a) Background Period",line=-1,adj=0.02,cex=.8)

par(fig=c(.1,.99,.39,.69),new=TRUE)
plot(1,1,type="l",ylim=c(1e-10,1e-4),log="y",xlim=c(0,55),xaxs="i",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
axis(2,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
abline(h=c(1e-10,1e-8,1e-6,1e-4),col="grey95")
lines(2:642/10+.8,ebl_means)
points(2,ebl_means[12],pch=8) #h2
points(c(12:17),ebl_means[c(112,122,132,142,152,162)],col=2,pch=8) #methane
points(c(26:30),ebl_means[c(252,262,272,282,292)],col="purple",pch=8) #ethane
legend("topright",expression(CH[4],C[2]*H[6],H[2]),col=c(2,"purple",1),pch=8,cex=.7,bty="n")
mtext("b) Ebullition",line=-1,adj=.02,cex=.8)

par(fig=c(.1,.99,.09,.39),new=TRUE)
plot(1,1,type="l",ylim=c(1e-10,1e-4),log="y",xlim=c(0,55),xaxs="i",
     frame.plot = TRUE,axes=FALSE,ann=FALSE)
axis(2,c((1:9)*10^-10,(1:9)*10^-9,(1:9)*10^-8,(1:9)*10^-7,(1:9)*10^-6,(1:9)*10^-5,10^-4),
     c(expression(10^{-10}),rep("",17),expression(10^{-8}),rep("",17),expression(10^{-6}),rep("",17),expression(10^{-4})),
     cex=.8,las=1)
mtext("Partial Pressure (Torr)",2,1.5,cex=.8)
abline(h=c(1e-10,1e-8,1e-6,1e-4),col="grey95")
lines(2:642/10+.8,wp_means)
points(2,wp_means[12],pch=8) #h2
points(c(12:17),wp_means[c(112,122,132,142,152,162)],col=2,pch=8) #methane
points(c(26:30),wp_means[c(252,262,272,282,292)],col="purple",pch=8) #ethane
legend("topright",expression(CH[4],C[2]*H[6],H[2]),col=c(2,"purple",1),pch=8,cex=.7,bty="n")
mtext("c) Peak Ebullition",line=-1,adj=.02,cex=.8)

axis(1)
mtext("m/z",1,1,font = 3)

dev.off()
