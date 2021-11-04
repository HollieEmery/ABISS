library(tidyverse)
library(lubridate)

path <- "~/Dropbox/Harvard/ABISS/RV Falkor (September 2018; McArthur Ridge)/parsed_4_R/"

# fixing the time stamps ####
# Timestamps for ABISS data are apparently incorrect in two distinct ways. 
# Not sure why, suspect it might be related to instrument warm-up time in one case,
# and a time zone and/or daylight savings time issue in the other:
# logging data in OR (PDT) but programming the system (EDT/EST?) or analyzing the data (EDT) in MA.

# Details murky because the assignment of date-time stamps is a black box. 
# Per Dan - it happens in WHOI-programmed processing in the instruments' brains - he didn't know details.

## My (HEE) confidence that this timestamp is erroneously shifted by 4 hrs is very high:
## 1) ABISS optode logged <20 mM O2 and ~4*C (ie seafloor conditions) at "10:36 GMT",
##    but launch didn't happen until 13:00 GMT
## 2) Dan's notes show an accidental power-off and immediate power-on at ~17:35 GMT (~10:35 PDT), 
##    this corresponds to a break in the RGA and optode logs at timestamp 13:34:34 (and this is the only log break on 15-sept-2018)
## 3) Dan's notes describe periods where the lander and wand were being moved by the ROV;
##    shifting by four hours makes periods of high temps correspond exactly with periods the lander was on the move,
##    and it aligns the time where the wand was being positioned in the bubble stream to a burst of noise from the RGA

## My (HEE) confidence that the RGA timestamp is further shifted by ~10 mins 20 sec is also high:
## 1) Both of the optode and RGA logs BEGIN at the same times, down to the second: 14:35:49 and 17:34:34, 
##    but finish 10:18 and 10:25 apart: RGA at 17:23:02 and 07:14:47, optode at 17:33:20 and 07:25:12.
## 2) Optode log gap at the accidental shutoff (~17:35 in Dan's notes) is 1:14, RGA gap is 11:32.
## 3) There is no plausible reason for the RGA to stop logging 10 minutes before an *accidental* shutoff.
## 4) A more plausible scenario is that the RGA has a 10 minute start up delay before data is actually recorded. 
##    Perhaps time stamps are incremented from the startup time, not recorded in real time.
## 5) Thus the optode and RGA logs should be aligned such that they END at the same time. 

time_shift = 4*3600 # four hours
rga_shift = 10*60+21 # 10 minutes and 21 seconds (~avg of 10:18 and 10:25)

#optode compensation function ####
# compensates for pressure, salinity, and temperature

#AADI optode salinity and depth compensation ported to MATLAB (DRH) then R (HEE) based on AADI correction Excel sheet.

#Dan's notes:
#Note 1. Solubility and salinity compensation calculation based on Garcia and Gordon. 1992. 
#Oxygen solubility in seawater: Better fitting equations. Limnology and Oceanography: 37(6) :1307-1312.

#Note 2. Pressure compensation based on Hiroshi Uchida, Takeshi Kawano, Ikuo Kaneko and Masao Fukasawa. 
#In-Situ calibration of optode-based oxygen sensors. Journal of Atmospheric and Oceanic Technology December 2008

#Noted 3. For other calculations, refer to AADI Operating Manual TD218 and TD269

correct_o2 <- function(o2, temp, depth, salinity=34.35, op_set_sal=0){
  #default salinity = 34.35 (value in PSU)
  #op_set_sal = 0. set salinity on the optode, factory default is 0.
  #depth = depth (m) of measurements.
  #o2 and temp are optode output
  
  #pars
  B0 = -0.00624097
  B1 = -0.00693498
  B2 = -0.00690358
  B3 = -0.00429155
  C0 = -0.00000031168
  p_comp = 0.032
  p_comp_factor = ((abs(depth)/1000)*p_comp)+1
  
  scale_temp <- log((298.15-dat$temp_C) / (273.15+dat$temp_C))
  sal_comp <- exp((salinity-op_set_sal)*(B0+B1*scale_temp+B2*scale_temp^2+B3*scale_temp^3) + C0*(salinity^2-op_set_sal^2))
  
  comp_o2 = dat$O2_conc_uM*sal_comp*p_comp_factor
  
  return(comp_o2)
}

# optode 1 ####
fns <- dir(path, pattern = "^o2_1", full.names=FALSE, recursive=TRUE) 
dat <- tibble()
for(i in 1:2){ #length(fns)
  datestamp = substr(fns[i],6,13)
  timestamp = substr(fns[i],14,19)
  foo <- read_csv(paste0(path,fns[i])) %>%
    separate(Time,into = c('elpsd_seconds','units'),sep=" ") %>%
    mutate(datetime = ymd_hms(paste0(datestamp,timestamp))+as.numeric(elpsd_seconds)) %>%
    select(datetime,O2_conc_uM =`O2Concentration[uM]`,air_sat_pct = `AirSaturation[%%]`,temp_C = `Temperature[Deg.C]`)
  dat <- rbind(dat,foo)
}
dat$O2_uM_corr = correct_o2(dat$O2_conc_uM,dat$temp_C,832)
dat$datetime = dat$datetime + time_shift
write.csv(dat,"~/Dropbox/Harvard/ABISS/Analysis/FinalData/optode1.csv",row.names=FALSE)

# optode 2 ####
fns <- dir(path, pattern = "^o2_2", full.names=FALSE, recursive=TRUE) 
dat <- tibble()
for(i in 1:2){ #length(fns)
  datestamp = substr(fns[i],6,13)
  timestamp = substr(fns[i],14,19)
  foo <- read_csv(paste0(path,fns[i])) %>%
    separate(Time,into = c('elpsd_seconds','units'),sep=" ") %>%
    mutate(datetime = ymd_hms(paste0(datestamp,timestamp))+as.numeric(elpsd_seconds)) %>%
    select(datetime,O2_conc_uM =`O2Concentration[uM]`,air_sat_pct = `AirSaturation[%%]`,temp_C = `Temperature[Deg.C]`)
  dat <- rbind(dat,foo)
}
dat$O2_uM_corr = correct_o2(dat$O2_conc_uM,dat$temp_C,832)
dat$datetime = dat$datetime + time_shift
write.csv(dat,"~/Dropbox/Harvard/ABISS/Analysis/FinalData/optode2.csv",row.names=FALSE)

# rga ####
fns <- dir(path, pattern = "^ra", full.names=FALSE, recursive=TRUE) 
dat <- tibble()
for(i in 1:2){ #length(fns)
  datestamp = substr(fns[i],4,11)
  timestamp = substr(fns[i],12,17)
  foo <- read_csv(paste0(path,fns[i])) %>%
    separate(Time,into = c('elpsd_seconds','units'),sep=" ") %>%
    mutate(datetime = ymd_hms(paste0(datestamp,timestamp))+as.numeric(elpsd_seconds)) %>%
    select(datetime,contains("MassSpectrum"))
  dat <- rbind(dat,foo)
}
dat[,-1] <- dat[,-1]/1.72e12 # per Dan: this is according to instrument manual
dat$total <- rowSums( dat[,c(-1,-643)])
dat$datetime = dat$datetime + time_shift + rga_shift
write.csv(dat,"~/Dropbox/Harvard/ABISS/Analysis/FinalData/rga.csv",row.names=FALSE)

