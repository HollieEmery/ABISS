library("tidyverse")
library("lubridate")

#data import ####

#ABISS data - data_process.R output
ox1 <- read_csv("FinalData/optode1.csv")
ox2 <- read_csv("FinalData/optode2.csv")
rga <- read_csv("FinalData/rga.csv")

#calibrated ch4 data - calibration.R output
conc_dat0 <- read_csv("FinalData/gas_conc.csv")

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

# inferred by looking at the data
wp <- c(13332:13451,13870:13943) #weird peaks
quiet1 <- 4506:4893
quiet2 <- 6110:6686
quiet3 <- 12101:12867
quiet <- c(quiet1,quiet2,quiet3)

wp_start = rga$datetime[c(13332,13870)]
wp_end = rga$datetime[c(13451,13943)]
quiet_start = rga$datetime[c(4506,6110,12101)]
quiet_end = rga$datetime[c(4893,6686,12867)]

# indices of phases of deployment 
background = which(rga$datetime >= power_blip+20*60 & rga$datetime < ROV_pu[2]) #lander in the mud while ROV searched for site, starting 20 mins after restart
active_onsite = which(rga$datetime >= wand_in_pos)
positioning = which((rga$datetime >= ROV_pu[1] & rga$datetime < ROV_pd[1]) | (rga$datetime >= ROV_pu[2] & rga$datetime < wand_pos[2]))
ebullitionA = which(rga$datetime >= wand_in_pos & rga$datetime < ymd_hm(201809160200))
ebullitionB = which(rga$datetime >= ymd_hm(201809160945) & rga$datetime < ymd_hm(201809162145))
ebullitionC = which(rga$datetime >= ymd_hm(201809170145) & rga$datetime < ymd_hm(201809170730))
ebullitionC[ebullitionC %in% wp] = NA
ebullitionC <- na.omit(ebullitionC)
ebullition = c(ebullitionA,ebullitionB,ebullitionC)


# calculate moving average ####
dat = rga # 1 minute simple moving average
for(i in 2:644){
  dat[,i] <- stats::filter(dat[,i],rep(1/7,7),sides = 2) 
}


# objects with only the good parts ####
dat1 <- dat 
dat1[-active_onsite,-1] <- NA
dat1[background,] <- dat[background,]

conc_dat <- conc_dat0
conc_dat[is.na(dat1$total),-1] <- NA

ox1[ox1$datetime < power_blip+600 | (ox1$datetime >= ROV_pu[2] & ox1$datetime<ROV_clear),-1] <- NA
ox2[ox2$datetime < power_blip+600 | (ox2$datetime >= ROV_pu[2] & ox2$datetime<ROV_clear),-1] <- NA

#take out weird peaks
dat2 <- dat1 
dat2[wp,-1] <- NA

# pull out relevent vectors for gases  ####
datetime <- dat1$datetime#[-1:-100]
ch4 <- rowMeans(dat1[,142:144])
ch4_cal <- conc_dat$methane_uM
total <- dat1$total

ox_rnd <- ox1 %>%
  mutate(mintime = round_date(datetime,unit="minute")) %>%
  group_by(mintime) %>%
  summarize(o2_uM = mean(O2_uM_corr))

rga_rnd <- tibble(datetime,ch4,ch4_cal,total) %>%
  mutate(mintime = round_date(datetime,unit="minute")) %>%
  group_by(mintime) %>%
  summarize(methane = mean(ch4_cal), rga = mean(total)) %>%
  left_join(ox_rnd,by="mintime") %>%
  mutate(phase = case_when(mintime >= power_blip+20*60 & mintime < ROV_pu[2] ~ "background",
                           mintime >= wp_start[1] & mintime < wp_end[1] ~ "free_gas",
                           mintime >= wp_start[2] & mintime < wp_end[2] ~ "free_gas",
                           mintime >= ymd_hm(201809151951) & mintime < ymd_hm(201809160200) ~ "ebullitionA",
                           mintime >= ymd_hm(201809160945) & mintime < ymd_hm(201809162145) ~ "ebullitionB",
                           mintime >= ymd_hm(201809170145) & mintime < ymd_hm(201809170730) ~ "ebullitionC",
                           mintime >= quiet_start[1] & mintime < quiet_end[1] ~ "quiet",
                           mintime >= quiet_start[2] & mintime < quiet_end[2] ~ "quiet",
                           mintime >= quiet_start[3] & mintime < quiet_end[3] ~ "quiet",
                           TRUE ~ "other"
                           ))

rga_rnd %>% 
  filter(phase!="background") %>%
  filter(phase!="free_gas") %>%
  arrange(phase) %>%
  ggplot(aes(o2_uM,methane)) +
  geom_smooth(method=lm,color="black") +
  #scale_y_log10() +
  geom_point(aes(o2_uM,methane,color=phase))

rga_rnd %>% 
  filter(phase!="background") %>%
  filter(phase!="free_gas") %>%
  arrange(phase) %>%
  ggplot(aes(o2_uM,methane,color=phase)) +
  geom_point() +
  geom_smooth(method=lm)

summary(lm(o2_uM~methane,rga_rnd,subset = phase!="free_gas" & phase!="background"))  



