##########################################################################
##########################################################################
#### Conductivity Calibration script for HOBO Conductivity-Temperature logger data
#### Brings in raw .csv files by Serial number and exports a tidy file

#### Reference: https://hasenmuellerlab.weebly.com/uploads/3/1/8/7/31874303/2019_shaughnessy_et_al_ema.pdf
#### Reference: https://www.aqion.de/site/112
#### Reference: https://oceanobservatories.org/wp-content/uploads/2015/10/1341-00040_Data_Product_SPEC_PRACSAL_OOI.pdf

# Author: Danielle Barnas  
# created: 9-23-2020
# modified: 5-2-2022

##########################################################################
##########################################################################

### Brings in calibration start and stop times and logger start and stop times in one large csv log
### Avoids need to input a particular serial number unless otherwise desired

###################################
### Load Libraries
###################################

#devtools::install_github("dbarnas/mooreasgd") # if package has updated since last run
library(tidyverse)
library(lubridate)
library(gsw)
library(here)
library(gridExtra)
#library(mooreasgd)


###################################
### File Paths
###################################

### Input
# Path to folder storing logger .csv files
path.log<-here("Data","CTDfiles","Treatment_CTDs", "Spring2023", "Cond") # Logger in situ file path (CT and Water Level files)
#path.WL<-here("Data","May2021","Depth")
file.date <- Sys.Date() # date used in naming output file(s)
hobo.csv <- FALSE # TRUE if csv has been processed and calibrated through HOBOware
csv.pattern <- "csv$" # file identifier at path.log (ex. "csv$")
#ct.serial <- "332" # if isolating one CT logger

### Output
# Path to store logger files
path.output<-here("Data","CTDfiles","Treatment_CTDs","QC_files") # Output file path


###################################
### Logger Launch and Retrieval dates
###################################

# Log dates
start.date <- ymd_hm('2023-02-04 14:00')
end.date <- ymd_hm('2023-03-28 12:00')

mid.cal.date <- ymd_hm('2023-03-09 14:00')


###################################
### Import calibration and launch records
###################################

# Read in files that are updated with calibration and launch information
calibration.log<-read_csv(here("Data","CT_Calibration_Log.csv")) # Calibration time logs
launch.log<-read_csv(here("Data","Launch_Log.csv")) # Launch time logs


###################################
### Pressure data
###################################

### Pressure at Mo'orea sites changes salinity nominally, so using consistent Absolute Pressure (dbar) for Salinity calculation
Pres_dbar<-10 # surface pressure in decibar


#################################################################################
# DO NOT CHANGE ANYTHING BELOW HERE ----------------------------------
#################################################################################

############################################################
### Read in Calibration and Logger Files
### Temperature compensation (source: https://www.aqion.de/site/112)
############################################################

# CT_cleanup function pulled from 'mooreasgd' package
# Reads in all raw csv's and returns multiple tidied csv's

# Conductivity Calibration files, if different file path from in situ logs
# condCal<-CT_cleanup(data.path = path.cal, path.pattern = c(file.date,"csv$"), tf.write = F)

# In Situ Conductivity files
source(here("Scripts", "mooreasgd_R", "CT_cleanup.R"))
condLog<-CT_cleanup(data.path = path.log, output.path = path.output, 
                    path.pattern = csv.pattern, tf.write = FALSE, tf.recursive = FALSE)

# focus on treatment locations (not seep)
logger.remove <- str_detect(condLog$LoggerID, "353")
condLog <- condLog %>% 
  cbind(logger.remove) %>% 
  filter(logger.remove == FALSE)

# Run if temperature may be in farenheit
condLog <- condLog %>%
  mutate(TempInSitu = if_else(TempInSitu > 50, ((TempInSitu - 32) *5/9), TempInSitu)) # crude F -> C; mutate if temp is higher than anticipated celcius in field

############################################################
### Parse date and time
############################################################

### Parse dates into date and type vector types

## Calibration
## unite date to time columns and parse to POSIXct datetime
calibration.log <- calibration.log %>% 
  filter(LoggerID == "323" | LoggerID == "320") %>% # filter out the best loggers now
  unite(col = 'time_in', date,time_in, sep = " ", remove = F) %>% # unite while maintaining separate date column
  unite(col = 'time_out', date,time_out, sep = " ", remove = F) %>% 
  mutate(date = mdy(date),
         time_in = mdy_hms(time_in),
         time_out = mdy_hms(time_out),
         LoggerID = as.character(LoggerID)) %>% 
  filter(time_in == start.date & pre_post == "pre" | 
           time_out == end.date & pre_post == "post") %>%  # filter only current pre- and post-calibration dates
  select(-notes)


# if selecting single CT from same calibration date
if(exists('ct.serial') == T){
  calibration.log <- calibration.log %>% 
    filter(LoggerID == ct.serial)
}

## in situ
launch.log <- launch.log %>%
  mutate(Date_launched = mdy_hm(time_start), # parse to date-time format
         Date_retrieved = mdy_hm(time_end)) %>% 
  filter(Date_launched >= start.date & Date_retrieved <= end.date) %>% 
  separate(Serial,into = c(NA,'LoggerID'), sep ="_", remove = F) %>% 
  filter(LoggerID == "323" | LoggerID == "320") %>% # filter out the best loggers now
  right_join(calibration.log, by = 'LoggerID') %>%  # select only logger ID's in calibration.log
  select(Serial, LoggerID, Date_launched, Date_retrieved, Log_Type)
  
#filter(Date_launched == start.date & Date_retrieved == end.date)# %>% # filter only current launch dates
  #unite(col = "Time_launched",Date_launched,Time_launched, sep = " ", remove = F) %>% # reunite date and time columns
  #unite(col = "Time_retrieved",Date_retrieved,Time_retrieved, sep = " ", remove = F) 
# launch.log <- launch.log %>% 
#   filter(CowTagID != "Sled") #%>% 
  #filter(LoggerID != "351") # Varari 3/18/2022


# if selecting single CT from same calibration date
if(exists('ct.serial') == T){
  launch.log <- launch.log %>% 
    filter(LoggerID == ct.serial)
}

if(hobo.csv == FALSE){
############################################################
### Calibration
############################################################

source(here("Scripts", "mooreasgd_R", "CT_one_cal.R"))
source(here("Scripts", "mooreasgd_R", "CT_two_cal.R"))
  
# get length of vector with all serial numbers for date(s) specified
n1 <- calibration.log %>%
  distinct(LoggerID) %>%
  nrow() %>%
  as.numeric()

# create empty df to add date-filtered data
CalLog<-tibble(date = as.POSIXct(NA),
               LoggerID = as.character(),
               FullLoggerID = as.character(),
               TempInSitu = as.numeric(),
               ECond.mS.cm = as.numeric(),
               Salinity_psu = as.numeric())

# Filter out calibration date and time and return dataframe with all logger calibration logs
for(i in 1:n1) {
  
  sn.a<-calibration.log %>% # vector of all serial numbers
    distinct(LoggerID) #%>% # pull each distinct Serial number
  #separate(col = 'LoggerID', into = c(NA,'LoggerID'), sep = "_")
  sn.a<-as.character(sn.a[i,]) # i'th serial number
  
  # filter by the i'th serial number
  C1.cal<-condLog %>% # C1: full in situ log; will be reduced to calibration logs
    filter(LoggerID == stringr::str_subset(condLog$LoggerID, pattern = sn.a))
  C2.cal<-calibration.log %>% # C2: shows calibration log relevant to LoggerID and datetime of calibration
    filter(LoggerID == stringr::str_subset(calibration.log$LoggerID, pattern = sn.a)) %>% 
    arrange(time_in) # order data by increasing date and time
  
  # pull out calibration times
  startHigh <- C2.cal %>%
    filter(cond_HL == 'H') %>% 
    select(time_in)
  endHigh <- C2.cal %>%
    filter(cond_HL == 'H') %>% 
    select(time_out)
  startLow <- C2.cal %>%
    filter(cond_HL == 'L') %>% 
    select(time_in)
  endLow <- C2.cal %>%
    filter(cond_HL == 'L') %>% 
    select(time_out)
  
  
    
    # specify whether calibration reference is an electrical conductivity or specific conductance (temp-compensated) value
    if(C2.cal$EC_SC == 'EC'){
      ec.cal = TRUE
    } else {ec.cal = FALSE}
    
    # single time point calibration
    calibration<-CT_one_cal(data = C1.cal,
                            cal.ref = C2.cal$cond_uS[1],
                            cal.ref.temp = C2.cal$temp_C[1],
                            startCal = startHigh,
                            endCal = endHigh,
                            date = C1.cal$date,
                            temp = C1.cal$TempInSitu,
                            EC.logger = C1.cal$E_Conductivity,
                            EC.cal = ec.cal) %>%
      select(date, LoggerID, TempInSitu_Cal, EC_Cal) %>% # only keep calibrated values and what we need
      rename(ECond.mS.cm = EC_Cal, # rename electrical conductivity column as .1 for first/only time point
             TempInSitu = TempInSitu_Cal)  # rename calibrated temperature readings, calibrated to secondary probe, if available
    

  
#####################################################
  
  
  # Calculate Practical Salinity using gsw package with PSS-78 equation
  calibration <- calibration %>%
    mutate(Salinity_psu = gsw_SP_from_C(C = ECond.mS.cm*0.001, t = TempInSitu, p = Pres_dbar))
  
  # Create column for FullLoggerID vs LoggerID
  calibration <- calibration %>% 
    rename(FullLoggerID = LoggerID) %>% 
    mutate(LoggerID = sn.a) %>% 
    select(date, LoggerID, FullLoggerID, TempInSitu, ECond.mS.cm, Salinity_psu)
  
  CalLog <- CalLog %>% 
    rbind(calibration) # add i'th logger's data to running dataframe
}

} else { # end of calibration if hobo.csv = FALSE
  CalLog <- condLog %>%
    separate(col = 'LoggerID', into = c(NA,'ID'), sep = "_", remove = FALSE) %>%  # needs to be standardized - relies on consistent file naming
    rename(LoggerID = ID,
           FullLoggerID = LoggerID,
           ECond.mS.cm = E_Conductivity)
} # End Initial Calibration


CalLog %>% 
  filter(Salinity_psu > 25) %>% 
  ggplot(aes(x = date, y = Salinity_psu, color = LoggerID)) +
  geom_point()+
  facet_wrap(~LoggerID)

###############################################################
 ### DO SECONDARY CALIBRATION FROM MID-POINT CAL
###############################################################
## where are we seeing this disconnect in each graph?
removeLog <- CalLog %>% 
  filter(between(date, ymd_hms("2023-03-24 15:00:00"), ymd_hms("2023-04-03 00:00:00")))
redCalLog <- CalLog %>% 
  anti_join(removeLog) %>% 
  filter(Salinity_psu > 25)

#### RECALIBRATE MID-LAUNCH WHEN VALUES JUMP FROM HANDLING IN THE FIELD
redCalLog %>% 
  filter(Salinity_psu < 40) %>% 
  #filter(between(date, ymd_hms("2023-03-02 00:00:00"), ymd_hms("2023-03-06 00:00:00"))) %>% 
  ggplot(aes(x = date, y = Salinity_psu, color = LoggerID)) +
  geom_point() +
  facet_wrap(~LoggerID, nrow = 2)



# reverse order of data before 3/9 for both loggers
preMidCal <- redCalLog %>% 
  filter(date < mid.cal.date)

reCalLog <- redCalLog %>% 
  filter(date >= mid.cal.date)

reCalLog <- rbind(preMidCal, reCalLog)

## read in calibration data for midcal
calibration.log.full<-read_csv(here("Data","CT_Calibration_Log.csv")) # Calibration time logs

calibration.log <- calibration.log.full %>% 
  filter(LoggerID == "323" | LoggerID == "320") %>% # filter out the best loggers now
  unite(col = 'time_in', date,time_in, sep = " ", remove = F) %>% # unite while maintaining separate date column
  mutate(date = mdy(date),
         time_in = mdy_hms(time_in),
         LoggerID = as.character(LoggerID)) %>% 
  select(-notes)

# temps taken at same time as EC after water sampling

calibration.log <- calibration.log %>% 
  mutate(pre_post = "post")
  


source(here("Scripts", "mooreasgd_R", "CT_one_cal.R"))
source(here("Scripts", "mooreasgd_R", "CT_two_cal.R"))

# get length of vector with all serial numbers for date(s) specified
n1 <- calibration.log %>%
  distinct(LoggerID) %>%
  nrow() %>%
  as.numeric()

# create empty df to add date-filtered data
midCalLog<-tibble(date = as.POSIXct(NA),
               LoggerID = as.character(),
               FullLoggerID = as.character(),
               TempInSitu = as.numeric(),
               ECond.mS.cm = as.numeric(),
               Salinity_psu = as.numeric())

# Filter out calibration date and time and return dataframe with all logger calibration logs
for(i in 1:n1) {
  
  sn.a<-calibration.log %>% # vector of all serial numbers
    distinct(LoggerID) #%>% # pull each distinct Serial number
  #separate(col = 'LoggerID', into = c(NA,'LoggerID'), sep = "_")
  sn.a<-as.character(sn.a[i,]) # i'th serial number
  
  # filter by the i'th serial number
  C2.cal<-calibration.log %>% # C2: shows calibration log relevant to LoggerID and datetime of calibration
    filter(LoggerID == stringr::str_subset(calibration.log$LoggerID, pattern = sn.a)) %>% 
    arrange(time_in) # order data by increasing date and time
  
  postC1.cal<-reCalLog %>% # C1: full in situ log; will be reduced to calibration logs
    filter(LoggerID == stringr::str_subset(reCalLog$LoggerID, pattern = sn.a)) %>% 
    filter(date >= mid.cal.date)
  postC2.cal <- C2.cal[2,]
  
  preC1.cal <- reCalLog %>% # C1: full in situ log; will be reduced to calibration logs
    filter(LoggerID == stringr::str_subset(reCalLog$LoggerID, pattern = sn.a)) %>% 
    filter(date <= mid.cal.date)
  preC2.cal <- C2.cal[1,]

  # pull out calibration times
  startHigh <- C2.cal %>%
    select(time_in)
  
  
  # single time point calibration (only for post ... pre was completed above)
  # Offset between the calibration reference and the logger reading
  ## POST-
  cal.ref <- postC2.cal$sal_psu # USING SALINITY FOR OFFSET INSTEAD OF EC 9/12/2023
  mean.ec <- postC1.cal %>% 
    filter(date == startHigh[2,1]) %>% 
    summarise(Salinity_psu) %>% # USING SALINITY FOR OFFSET INSTEAD OF EC 9/12/2023
    as.numeric()
  offset.ec<-cal.ref - mean.ec
  
  # Apply offset to logger data
  postC1.cald <- postC1.cal %>% 
    filter(date >= mid.cal.date) %>% 
    dplyr::mutate(Salinity_psu = Salinity_psu + offset.ec) # USING SALINITY FOR OFFSET INSTEAD OF EC 9/12/2023
    #mutate(Salinity_psu = gsw_SP_from_C(C = ECond.mS.cm*0.001, t = TempInSitu, p = Pres_dbar))
  postC1.cald %>% 
    ggplot(aes(x = date, y = Salinity_psu, color = LoggerID)) +
    geom_point() +
    facet_wrap(~LoggerID, nrow = 2)
  

  
 
  ### Drift Compensation up to march 9, 2023 14:00:00
  EC_Cal.1 <- C2.cal$sal_psu[1] # USING SALINITY FOR OFFSET INSTEAD OF EC 9/12/2023
  EC_Cal.2 <- C2.cal$sal_psu[2] # USING SALINITY FOR OFFSET INSTEAD OF EC 9/12/2023
  
  preC1.cal <- preC1.cal %>% 
    filter(between(date,startHigh[1,],startHigh[2,]))  # selects a calibration time point
  
  EC.drift.1 = EC_Cal.1 - preC1.cal$Salinity_psu[1] # offset from 1st value
  EC.drift.2 = EC_Cal.2 - preC1.cal$Salinity_psu[length(preC1.cal$Salinity_psu)] # offset from last value
  drift.cor.ec = (EC.drift.2 - EC.drift.1) / length(preC1.cal$Salinity_psu) # total units drifted over time / total units
  

  ## CALCULATE CORRECTED SALINITY WITH DRIFT
  
  preC1.cald <- preC1.cal %>% 
    filter(between(date,startHigh[1,],startHigh[2,])) %>% # selects all data between calibration time points
    arrange(date) %>%
    mutate(EC.drift.init = EC.drift.1,
           drift.cor = drift.cor.ec, # establish a column filled with the drift correction values
           drift.correction=cumsum(drift.cor), # fill the drift correction column with sequentially larger drift corrections from correlation value to full drift
           EC.drift.fin = EC.drift.init + drift.correction) %>%  # add the drift correction value to EC.drift.1 sequentially until we reach EC.drift.2 in final time point
    select(-drift.cor) %>% 
    mutate(Salinity_psu_corr = Salinity_psu + EC.drift.fin) %>%  # USING SALINITY FOR OFFSET INSTEAD OF EC 9/12/2023
    #mutate(ECond.mS.cm = EC_Cal.1 + drift.correction.ec) # USING SALINITY FOR OFFSET INSTEAD OF EC 9/12/2023
    #mutate(Salinity_psu = gsw_SP_from_C(C = ECond.mS.cm*0.001, t = TempInSitu, p = Pres_dbar))
    select(-Salinity_psu) %>% 
    rename(Salinity_psu = Salinity_psu_corr) %>% 
    filter(date != startHigh[2,1]) # remove final time point due to overlap with post-calibration file
  
  
  
  #####################################################
  
  calibration <- full_join(preC1.cald, postC1.cald) %>% 
    arrange(LoggerID, date)
  
  midCalLog <- midCalLog %>% 
    rbind(calibration) # add i'th logger's data to running dataframe
}

midCalLog <- midCalLog %>% 
  select(date,LoggerID, FullLoggerID, TempInSitu, ECond.mS.cm, Salinity_psu)

View(midCalLog)

rename <- c("320" = "Low SGD - 320", "323" = "High SGD - 323")

midCalLog %>% 
  ggplot(aes(x = date, y = Salinity_psu, color = LoggerID)) +
  geom_point() +
  facet_wrap(~LoggerID, labeller = as_labeller(rename))





############################################################
### In Situ Logger Data
############################################################

n2 <- launch.log %>%
  distinct(LoggerID) %>%
  nrow() %>%
  as.numeric()

# create empty df to add date-filtered data
Log <-tibble(date = as.POSIXct(NA),
             LoggerID = as.character(),
             FullLoggerID = as.character(),
             TempInSitu = as.numeric(),
             ECond.mS.cm = as.numeric(),
             Salinity_psu = as.numeric())

# create list for storing plots
p <- list()

# Pull out in situ logger data
for(i in 1:n2) {

  sn.b <- launch.log %>% # vector of all serial numbers
    distinct(LoggerID) #%>% 
  #separate(col = 'LoggerID', into = c(NA,'LoggerID'), sep = "_")
  sn.b <- as.character(sn.b[i,]) # i'th serial number
  
  # filter by the i'th serial number
  C1.log<-midCalLog %>% 
    filter(LoggerID == str_subset(midCalLog$LoggerID, pattern = sn.b)) #%>% 
  #mutate(LoggerID = paste0("CT_",sn.b)) # make serial the same for easy join
  C2.log<-launch.log %>% 
    filter(LoggerID == str_subset(launch.log$LoggerID, pattern = sn.b)) %>% 
    distinct()
  
  
  launch.start<-C2.log %>% 
    select(Date_launched)
  launch.end<-C2.log %>% 
    select(Date_retrieved)
  
  
  # pull out original CT serial name
  sn.b <- midCalLog %>% 
    filter(LoggerID == str_subset(midCalLog$LoggerID, pattern = sn.b)) %>% # filter for i'th sn
    select(LoggerID)
  sn.b <- as.character(sn.b[1,]) # character string of full serial name
  
  C1.log <- C1.log %>% 
    mutate(LoggerID = sn.b)
  
  
  
  C1.log <- C1.log %>% 
     filter(between(date, launch.start[1,1], launch.end[1,1]))
  
  # Create Plot and save to list p
  p[[i]] <- C1.log %>% 
    ggplot(aes(x = date, y = Salinity_psu, color = TempInSitu)) + 
    geom_point() + 
    theme_bw() +
    labs(x = "Date", color = "Temperature (C)", y = "Salinity (psu)") + 
    ylim(33,40) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7))+
    ggtitle(sn.b)
  
  
  #write.csv(C1.log, paste0(path.output,"/QC_",sn.b,".csv")) # write csv file for each Logger
  
  Log<-Log %>% 
    rbind(C1.log) %>% 
    distinct()
  
}


#############################
### SAVE FILE AND PLOTS
#############################

# write_csv(Log, paste0(path.output,"/Full_CT_",file.date,".csv")) # write csv file for full set of logger data
# 
# # Save all plots in a single dated pdf
# pdf(paste0(path.output,"/",file.date,"_CT_plots.pdf"), onefile = TRUE)
# for (i in seq(length(p))) {
#   tplot <- p[[i]]
#   print(tplot)
# }
# dev.off()


