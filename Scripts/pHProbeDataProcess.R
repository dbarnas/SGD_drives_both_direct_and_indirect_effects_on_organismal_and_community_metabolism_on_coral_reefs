## process pH for Orion
# Created by Nyssa Silbiger
# Edited on 07/27/2021



#############################
### LOAD LIBRARIES
#############################
library(tidyverse)
library(seacarb)
library(broom)
library(here)
library(lubridate)


#############################
### READ IN DATA
#############################

## bring in pH calibration files and raw data files
pHcalib<-read_csv(here("Data","Tris_Calibration_3_10_2023.csv")) %>%
  mutate(TrisCalDate = mdy(TrisCalDate))
pHData<-read_csv(here("Data","Water_Sampling.csv"))%>%
  mutate(TrisCalDate = mdy(TrisCalDate),
         Date = mdy(Date))



#############################
### PROCESS DATA
#############################

## take the mV calibration files by each date and use them to calculate pH
pHSlope<-pHcalib %>%
  nest_by(TrisCalDate)%>%
  mutate(fitpH = list(lm(mVTris~TTris, data = data))) %>% # linear regression of mV and temp of the tris
  summarise(broom::tidy(fitpH)) %>% # make the output tidy
  select(TrisCalDate, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%# put slope and intercept in their own column
  right_join(.,pHData) %>% # join with the pH sample data
  mutate(mVTris = temp.orion*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  mutate(pH = pH(Ex=mV.orion,Etris=mVTris,S=sal.ysi,T=temp.orion)) %>% # calculate pH of the samples using the pH seacarb function
  select(Date, Time, Site, Location,Day_Night, Tide, sal.ysi, TempInSitu, pH)


#############################
### SAVE FILE
#############################

#write_csv(x = pHSlope, file = here("Data","pHProbe_Data_calculated.csv"))


#############################
### SUMMARISE
#############################

library(plotrix)
pHSlope %>% 
  group_by(Location) %>% 
  summarise(meanSal = mean(sal.ysi),
            sesal = std.error(sal.ysi),
            minSal = min(sal.ysi),
            maxSal = max(sal.ysi),
            meanTemp = mean(TempInSitu),
            setemp = std.error(TempInSitu),
            minTemp = min(TempInSitu),
            maxTemp = max(TempInSitu),
            meanpH = mean(pH),
            sepH = std.error(pH),
            maxpH = max(pH),
            minpH = min(pH))
