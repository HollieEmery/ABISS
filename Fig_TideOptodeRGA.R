# Pub Figure : Time Series of Optodes, Tides, RGA ####

library("tidyverse")
library("lubridate")

#data import ####

#ABISS data
ox1 <- read_csv("~/Dropbox/Harvard/ABISS/Analysis/FinalData/optode1.csv")
ox2 <- read_csv("~/Dropbox/Harvard/ABISS/Analysis/FinalData/optode2.csv")
rga <- read_csv("~/Dropbox/Harvard/ABISS/Analysis/FinalData/rga.csv")

#tide data
sensor <- read_csv("~/Dropbox/Harvard/ABISS/Analysis/FinalData/S_Hydrate_Ridge_Pressure_Tides.csv") %>%
  mutate(datetime=mdy_hms(Date_Time,tz = "UTC")) 
astoria <- read_csv("~/Dropbox/Harvard/ABISS/Analysis/FinalData/CO-OPS_9439040_met_Astoria.csv") %>%
  mutate(datetime=mdy_hms(paste(Date,Time_GMT),tz = "UTC")) #46.2067 -123.7683
southbeach <- read_csv("~/Dropbox/Harvard/ABISS/Analysis/FinalData/CO-OPS_9435380_met_SouthBeach.csv") %>%
  mutate(datetime=mdy_hms(paste(Date,Time_GMT),tz = "UTC")) #44.625 -124.045
garibaldi <- read_csv("~/Dropbox/Harvard/ABISS/Analysis/FinalData/CO-OPS_9437540_met_Garibaldi.csv") %>%
  mutate(datetime=mdy_hms(paste(Date,Time_GMT),tz = "UTC")) #45.555 -123.9183

tides <- sensor %>%
  mutate(datetime_hm = floor_date(datetime,"minute")) %>%
  group_by(datetime_hm) %>%
  summarise(Temp = mean(Temperature_C), Depth = mean(Depth_m))

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

# subsections ####
dat_wp <- dat_bk <- dat_eb <- dat_qu <- dat
dat_wp[-wp,-1] <- NA
dat_bk[-background,-1] <- NA
dat_eb[-ebullition,-1] <- NA
dat_qu[-quiet,-1] <- NA

dat_other <- dat
dat_other[-active_onsite,-1] <- NA
dat_other[wp[3:length(wp)-1],-1] <- NA
dat_other[ebullitionA[3:length(ebullitionA)-1],-1] <- NA
dat_other[ebullitionB[3:length(ebullitionB)-1],-1] <- NA
dat_other[ebullitionC[3:length(ebullitionC)-1],-1] <- NA
dat_other[quiet1[3:length(quiet1)-1],-1] <- NA
dat_other[quiet2[3:length(quiet2)-1],-1] <- NA
dat_other[quiet3[3:length(quiet3)-1],-1] <- NA

ox1[ox1$datetime < power_blip+600 | (ox1$datetime >= ROV_pu[2] & ox1$datetime<ROV_clear),-1] <- NA
ox2[ox2$datetime < power_blip+600 | (ox2$datetime >= ROV_pu[2] & ox2$datetime<ROV_clear),-1] <- NA

total = dat$total
total[dat$datetime < power_blip+600 | (dat$datetime >= ROV_pu[2] & dat$datetime<ROV_clear)] <- NA
  
#figure area
start = ymd_hm(201809151430)
end = ymd_hm(201809170900)
Time = ymd_hm(c(201809151200,201809151600,201809152000,201809160000,201809160400,201809160800,201809161200,201809161600,201809162000,201809170000,201809170400,201809170800))
time_names = c("12:00","16:00","20:00", "16-Sep","4:00","8:00","12:00","16:00","20:00","17-Sep","4:00","8:00")

arrow_starts = ymd_hm(201809151951,201809160945,201809170145)
arrow_ends = ymd_hm(201809160200,201809162145,201809170730)
diffs = arrow_ends-arrow_starts

setEPS()
#postscript("~/Dropbox/Harvard/ABISS/R_plots/TimeSeries.eps",height=10,width=7.5)
par(mfrow=c(1,1),mai=c(0, 0, 0, 0),tcl=.25,mgp=c(3,.1,0),cex.axis=.8,las=1)
par(fig=c(.08,.99,.75,.99))
plot(Verified_m~datetime,data=garibaldi,subset=datetime>power_blip,type="l",col=4,
     xlim=c(start,end),ylim=c(0,3.5),
     axes=FALSE,ann=FALSE,frame.plot=TRUE,xaxs="i")
lines(Verified_m~datetime,data=astoria,subset=datetime>power_blip,col=2)
lines(Verified_m~datetime,data=southbeach,subset=datetime>power_blip,col=3)
axis(2,line=-2,at=c(0,1,2,3))
mtext("Tide MLLW (m)",2,-1.5,las=0)
lines(Depth-792~datetime_hm,data=tides,subset=datetime_hm>power_blip+60*10,col=1)
axis(2,at=c(0,1,2,3),label=paste(792+c(0,1,2,3)))
mtext("OOI depth (m)",2,1.5,las=0)
legend("topright",c("Astoria","South Beach","Garibaldi","OOI station"," "),lwd=2,col=c(2:4,1,NA),ncol=2,cex=.8,bty="n")
mtext("a)",line=-1,adj=.02)

par(fig=c(.08,.99,.56,.74),new=TRUE)
plot(start,0,type="l",col=1,
     xlim=c(start,end),ylim=c(3.85,4.2),
     axes=FALSE,ann=FALSE,frame.plot=TRUE,xaxs="i")
lines(ox1$datetime,ox1$temp_C,col=1,lty=1)
lines(ox2$datetime,ox2$temp_C,col="grey",lty=1)
axis(2,c(3.9,4,4.1,4.2))
mtext("Temperature (°C)",2,1.5,las=0)
legend("topright",c("At seep","1.5 m above sf"),lty=1,col=c(1,"grey"),cex=.8,bty="n")
mtext("b)",line=-1,adj=.02)

par(fig=c(.08,.99,.36,.55),new=TRUE)
plot(start,0,type="l",col=1,
     xlim=c(start,end),ylim=c(7,15),#ylim=c(7,30),
     axes=FALSE,ann=FALSE,frame.plot=TRUE,xaxs="i")
lines(ox1$datetime,ox1$O2_uM_corr,col=1)
lines(ox2$datetime,ox2$O2_uM_corr,col="grey",lty=1)
axis(2)#axis(2,c(8,10,12,14))
mtext(expression(O[2]~(mu*M)),2,1.5,las=0)
legend("topright",c("At seep","1.5 m above sf"),lty=1,col=c(1,"grey"),cex=.8,bty="n")
mtext("c)",line=-1,adj=.02)

par(fig=c(.08,.99,.07,.35),new=TRUE)
plot(dat_other$datetime,apply(dat_other[,2:642],1,sum),type="l",
     xlim=c(start,end),ylim=c(5e-6,1e-3),log="y",
     axes=FALSE,ann=FALSE,frame.plot=TRUE,xaxs="i")
lines(dat_wp$datetime,apply(dat_wp[,2:642],1,sum),col="red")
lines(dat_bk$datetime,apply(dat_bk[,2:642],1,sum),col="grey")
lines(dat_eb$datetime,apply(dat_eb[,2:642],1,sum),col=4)
lines(dat_qu$datetime,apply(dat_qu[,2:642],1,sum),col=3)

arrows(arrow_starts,4e-4,x1=arrow_ends,code=3,length=.08)
text(arrow_starts+diffs/2,
     6e-4,"Ebullition")
arrows(ymd_hm(201809151815),2e-5,y1=1e-5,length=.08)
text(ymd_hm(201809151815),2.5e-5,"Background")
axis(1,at=Time,labels=time_names)
axis(2,at=c(1e-3,1e-4,1e-5,1e-6),
     labels=c(expression(10^{-3}),expression(10^{-4}),expression(10^{-5}),expression(10^{-6})))
mtext(expression(Sigma~Partial~Pressure~(Torr)),side=2,line=1.5,las=0)
mtext("d)",line=-1,adj=.02)

dev.off()
