library("readxl")
library("tidyverse")
library("lubridate")

#load henry's law function for easy calculation
source("henrys_law.R")

# what to go with... ####
# use only the closest two pressures (1000 psi, 1500 psi, approx depths 690, 1034 m) 
# - The most comparable, 832 is roughly in the middle of these (-143 & +202) 
# - There seems to be a substantial effect of pressure: a reason to exclude less-relevant ones
# - Looking at calibration day time series, there seems to be an AM "warm up": extra reason to exclude the first 1-2 "depths"

# use only 4 C
# - site is 4 with very little variation
# - still plenty of points in the curve
# - no real disagreements in the shape of the curve if you exclude or include other temps

# use only concentration ranges bounding field measurements
# - cuts down on possible non-linearity at high concentrations

# use ratio of index fragment to 18.0
# - normalizing to water deals with some of the systemic drift & noise - see Scott Wankel's GoM paper
# - some analytes are better with .1 or .2, but let's keep it simple and consistent
## note exception: o2 was simply the index fragment because it fit better - O2 is very weak anyway

# use linear fits
# - since only 1 temp and 2 pressures, I don't see a need for complex polynomial double fits
# - R2s are good: >.98 for all except oxygen which is only .86
## note exception: methane (described below)

# fit methane line to log/log transformations instead of linear
# - reduces the leverage of the higher, messier measurements
# - I dealt with the zeroes by setting them to 10
# -- 0s on day 6 clearly not measurements (all exactly "0"), it was probably either below the detection limit of the instrument or an assumption based on eg sparging
# -- empirically, 10 is best (curve R2 .9885 instead of .9556), also tried 0, 1, 5, 15, 8 etc

# exclude calibration day 1 for methane
# - data is kinda fishy, perhaps they were still figuring out technique?
# - removing it doesn't cut out any range target value on the curve


#methane calibration ####
cal_methane <- read_xlsx("~/Dropbox/Harvard/ABISS/cal_data.xlsx",sheet = "methane") %>%
  filter(Depth_approx_m > 500 & Depth_approx_m <1300) %>%
  filter(Calibration_day>1) %>%
  #filter(conc_GCMS_uM>0) %>%
  filter(Temp_C==4) %>%
  mutate(Calibration_day = as.character(Calibration_day),
         ch4=(ISMS_mz15+ISMS_mz15.1+ISMS_mz15.2)/3,
         h2o=(ISMS_mz18+ISMS_mz18.1+ISMS_mz18.2)/3,
         ISMS_15.18 = ch4/h2o,#ISMS_mz15/ISMS_mz18,
         conc_transform = log10(ifelse(conc_GCMS_uM<=10,10,conc_GCMS_uM)),
         #conc_transform = log1p(conc_GCMS_uM)
         #signal_transform = log10(ISMS_mz15),
         signal_transform = log10(ISMS_15.18))
  
ggplot(cal_methane) +
  geom_point(aes(y=ISMS_mz15,x=ifelse(conc_GCMS_uM<=10,10,conc_GCMS_uM),
                 col=P_psi,shape=as.character(Temp_C))) +
  scale_x_log10() +
  scale_y_log10()

ggplot(cal_methane) +
  geom_point(aes(y=ISMS_15.18,x=ifelse(conc_GCMS_uM<=10,10,conc_GCMS_uM),
                 col=P_psi,shape=as.character(Temp_C))) +
  scale_x_log10() +
  scale_y_log10()

#fit_methane <- lm(ISMS_mz15~conc_GCMS_uM,data=cal_methane)
fit_methane <- lm(signal_transform~conc_transform,data=cal_methane)

summary(fit_methane) 

# linear           n=48, R2=.8229 
# use 15/18 ratio: n=48, R2=.8547
# 4C only:         n=16, R2=.8423
# omit cal day 1:  n=42, R2=.8241
# 4C + omit day 1: n=14, R2=.8458

# log plot is better 

# dealing with zeros:
# omit zeros:      n=42, R2=.9406 
# log(x+1) to all: n=48, R2=.9476 
# set zeros to 1:  n=48, R2=.9481
# set zeros to 5:  n=48, R2=.9709
# set zeros to 10: n=48, R2=.9713

# omit cal day 1
# set zeros to 1:  n=42, R2=.9528
# set zeros to 5:  n=42, R2=.9788
# set zeros to 10: n=42, R2=.9811

# 4C only
# set zeros to 1:  n=16, R2=.9477
# set zeros to 5:  n=16, R2=.9707
# set zeros to 10: n=16, R2=.9706

# omit cal day 1; 4C only
# set zeros to 1:  n=14, R2=.9544
# set zeros to 5:  n=14, R2=.9813
# set zeros to 10: n=14, R2=.9833

#use 15/18 ratio
# set zeros to 1:  n=48, R2=.9424
# set zeros to 5:  n=48, R2=.965 
# set zeros to 10: n=48, R2=.9655 
# set zeros to 7:  n=48, R2=.9663 

# use 15/18 ratio; 4C only
# set zeros to 1:  n=16, R2=.9413 
# set zeros to 5:  n=16, R2=.9646
# set zeros to 10: n=16, R2=.9648 

# use 15/18 ratio; omit cal day 1
# set zeros to 1:  n=42, R2=.9547
# set zeros to 5:  n=42, R2=.9828
# set zeros to 10: n=42, R2=.9862

# use 15/18 ratio; omit cal day 1; 4C only
# set zeros to 1:  n=14, R2=.9556
# set zeros to 5:  n=14, R2=.9851
# set zeros to 10: n=14, R2=.9885 

#use averages of 15.0-15.2 and 18.0-18.2; use 15/18 ratio; omit cal day 1; 4C only: 
# set zeros to 10: n=14, R2=.9887 ******

#oxygen calibration ####
cal_oxygen <- read_xlsx("~/Dropbox/Harvard/ABISS/cal_data.xlsx",sheet = "oxygen") %>%
  filter(Depth_approx_m > 500 & Depth_approx_m <1300) %>%
  filter(Temp_C==4) %>%
  mutate(Calibration_day = as.character(Calibration_day),
         o2=(ISMS_mz32+ISMS_mz32.1+ISMS_mz32.2)/3,
         h2o=(ISMS_mz18+ISMS_mz18.1+ISMS_mz18.2)/3,
         ISMS_32.18 = o2/h2o,#ISMS_mz32/ISMS_mz18,
         conc_transform = log10(ifelse(conc_GCMS_uM<=1,1,conc_GCMS_uM)),
         signal_transform = log10(ISMS_mz32)) 

ggplot(cal_oxygen) +
  geom_point(aes(y=ISMS_mz32,x=ifelse(conc_GCMS_uM<=0,0,conc_GCMS_uM),
                 col=P_psi,shape=as.character(Temp_C))) #+

ggplot(cal_oxygen) +
  geom_point(aes(y=ISMS_32.18,x=conc_GCMS_uM,
                 col=P_psi,shape=as.character(Calibration_day)))

#fit_oxygen <- lm(signal_transform~conc_transform,data=cal_oxygen)
#fit_oxygen <- lm(ISMS_32.18~conc_GCMS_uM,data=cal_oxygen)
#fit_oxygen <- lm(ISMS_mz32~conc_GCMS_uM,data=cal_oxygen)
fit_oxygen <- lm(o2~conc_GCMS_uM,data=cal_oxygen)

plot(o2~conc_GCMS_uM,cal_oxygen)
abline(fit_oxygen,col="red")

summary(fit_oxygen) 
# all temps:
# log(x+1) to all: n=24, R2=.5506
# linear:          n=24, R2=.8226
# linear 32/18:    n=24, R2=.8598 

# 4C only:
# log(x+1) to all: n=8, R2=.7122
# linear:          n=8, R2=.8878
# linear 32/18:    n=8, R2=.77

# mean of 32.0-32.2; 4C only:
# linear:          n=8, R2=.8922 *********
# linear 32/18:    n=8, R2=.7511

#hydrogen calibration ####
cal_hydrogen <- read_xlsx("~/Dropbox/Harvard/ABISS/cal_data.xlsx",sheet = "hydrogen") %>%
  filter(Depth_approx_m > 500 & Depth_approx_m <1300) %>%
  filter(conc_GCMS_uM < 2000) %>%
  filter(Temp_C==4) %>%
  mutate(Calibration_day = as.character(Calibration_day),
         h2=(ISMS_mz2+ISMS_mz2.1+ISMS_mz2.2)/3,
         h2o=(ISMS_mz18+ISMS_mz18.1+ISMS_mz18.2)/3,
         ISMS_2.18 = h2/h2o,#ISMS_mz2/ISMS_mz18,
         conc_transform = log10(ifelse(conc_GCMS_uM<=1,1,conc_GCMS_uM)),
         signal_transform = log10(ISMS_mz2))

ggplot(cal_hydrogen) +
  geom_point(aes(y=ISMS_mz2,x=conc_GCMS_uM,col=P_psi,shape=as.character(Temp_C))) 
  
ggplot(cal_hydrogen) +
  geom_point(aes(y=ISMS_2.18,x=conc_GCMS_uM,col=P_psi,shape=as.character(Temp_C)))

#fit_hydrogen <- lm(ISMS_mz2~conc_GCMS_uM,data=cal_hydrogen)
fit_hydrogen <- lm(ISMS_2.18~conc_GCMS_uM,data=cal_hydrogen)
#fit_hydrogen <- lm(h2~conc_GCMS_uM,data=cal_hydrogen)

summary(fit_hydrogen) 
# all temps:
# linear:          n=36, R2=.9588
# linear 2/18:     n=36, R2=.9814

# 4C only:
# linear:          n=12, R2=.9768
# linear 2/18:     n=12, R2=.9913
# only up to 2 mM  n=8,  R2=.9914 ********

# 4C only, .1 fragments:
# linear 2/18:     n=12, R2=.9916 

# 4C only, .2 fragments:
# linear 2/18:     n=12, R2=.9903

# 4C only; mean of all fragments; only up to 2 mM:
# linear:          n=8, R2=.9898
# linear 2/18:     n=8, R2=.9889 

plot(ISMS_2.18~conc_GCMS_uM,cal_hydrogen)
abline(fit_hydrogen,col="red")

#co2 calibration ####
cal_co2 <- read_xlsx("~/Dropbox/Harvard/ABISS/cal_data.xlsx",sheet = "co2") %>%
  filter(Depth_approx_m > 500 & Depth_approx_m <1300) %>%
  filter(Temp_C==4) %>%
  filter(conc_GCMS_uM < 10000) %>%
  mutate(Calibration_day = as.character(Calibration_day),
         co2=(ISMS_mz44+ISMS_mz44.1+ISMS_mz44.2)/3,
         h2o=(ISMS_mz18+ISMS_mz18.1+ISMS_mz18.2)/3,
         ISMS_44.18 = co2/h2o,#ISMS_mz44/ISMS_mz18,
         conc_transform = log10(ifelse(conc_GCMS_uM<=1,1,conc_GCMS_uM)),
         signal_transform = log10(ISMS_mz44))

ggplot(cal_co2) +
  geom_point(aes(y=ISMS_mz44,x=conc_GCMS_uM,col=P_psi,shape=as.character(Temp_C))) 
  
ggplot(cal_co2) +
  geom_point(aes(y=ISMS_44.18,x=conc_GCMS_uM,col=P_psi,shape=as.character(Temp_C)))
  
fit_co2 <- lm(ISMS_mz44~conc_GCMS_uM,data=cal_co2)
fit_co2 <- lm(ISMS_44.18~conc_GCMS_uM,data=cal_co2)

summary(fit_co2) 
# all temps:
# linear:          n=36, R2=.9269
# linear 44/18:    n=36, R2=.9855

# 4C only:
# linear:          n=12, R2=.9643
# linear 44/18:    n=12, R2=.9947

# only up to 10 mM:
# linear:          n=6, R2=.9493
# linear 44/18:    n=6, R2=.9909

#mean of 44.0-44.2 fragments; only up to 10 mM; 4C only:
# linear 44/18:    n=6, R2=.9911 ****************

plot(ISMS_44.18~conc_GCMS_uM,cal_co2)
abline(fit_co2,col="red")

#ms data ####
dat <- read_csv("FinalData/rga.csv")
meth_sat. = henry(1e6,"methane",TC=4,z=832)*1000*1000 
dat1 <- dat
for(i in 2:644){
  dat1[,i] <- stats::filter(dat[,i],rep(1/7,7),sides = 2) 
}

conc_dat <- tibble(date_time = dat1$datetime,
                   #methane_uM = 10^((log10(pmax(0,dat1$MassSpectrum_142/dat1$MassSpectrum_172)) - fit_methane$coefficients[1])/fit_methane$coefficients[2]),
                   methane_uM = 10^((log10(pmax(0,rowMeans(dat1[,142:144])/rowMeans(dat1[,172:174]))) - fit_methane$coefficients[1])/fit_methane$coefficients[2]),
                   #hydrogen_uM = (dat1$MassSpectrum_12/dat1$MassSpectrum_172-fit_hydrogen$coefficients[1])/fit_hydrogen$coefficients[2],
                   hydrogen_uM = (rowMeans(dat1[,12:14])/rowMeans(dat1[,172:174])-fit_hydrogen$coefficients[1])/fit_hydrogen$coefficients[2],
                   #oxygen_uM = (dat1$MassSpectrum_312-fit_oxygen$coefficients[1])/fit_oxygen$coefficients[2],
                   oxygen_uM = (rowMeans(dat1[,312:314])-fit_oxygen$coefficients[1])/fit_oxygen$coefficients[2],
                   #co2_uM = (dat1$MassSpectrum_432/dat1$MassSpectrum_172-fit_co2$coefficients[1])/fit_co2$coefficients[2]
                   co2_uM = (rowMeans(dat1[,432:434])/rowMeans(dat1[,172:174])-fit_co2$coefficients[1])/fit_co2$coefficients[2],
                   ) %>%
  mutate(pct_sat = methane_uM/meth_sat.,
         free_gas = methane_uM>meth_sat.)

write_csv(conc_dat,"FinalData/gas_conc.csv")
