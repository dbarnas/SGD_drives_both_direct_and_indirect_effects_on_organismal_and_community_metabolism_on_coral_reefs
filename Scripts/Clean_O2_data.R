### Processing raw oxy-10 O2 files
### removes data points that are not used for analysis
### Created by Danielle Barnas
### Created on May 3, 2023


#############################
### LOAD LIBRARIES
#############################
library(tidyverse)
library(here)
library(lubridate)


#############################
### CREATE FUNCTIONS ###
#############################

# function to read in files, standard
readfilefun <- function(mypath, funfile){
  df <- read_csv(paste0(mypath,"/", funfile), skip = 1) %>% 
    unite(col = "datetime", Date, Time, remove = F) %>% 
    mutate(datetime = mdy_hms(datetime))
  return(df)
}

# function to read in file and filter by datetime
run2fun <- function(run2file){
  df <- read_csv(here("Data", "RespoFiles", "RawO2", "20230326", "RUN2LIGHT", run2file), skip = 1) %>% 
    unite(col = "datetime", Date, Time, remove = F) %>% 
    mutate(datetime = mdy_hms(datetime)) %>% 
    filter(datetime >= ymd_hms("2023-03-26 11:48:00"))
  return(df)
}

rawfilefun <- function(rawpath, cleanpath, funfile){
  df <- read_csv(paste0(rawpath, "/", funfile), skip = 1) %>% 
    unite(col = "datetime", Date, Time, remove = F) %>% 
    mutate(datetime = mdy_hms(datetime))
  #write_csv(df, paste0(cleanpath,"/", funfile))
}

# Broad write path
wpath <- "Data/RespoFiles/CleanO2"



#############################
### READ IN DATA AND PROCESS INDIVIDUALLY ###
#############################


### RUN 2 LIGHT
# ignore channel 6 in this run due to later rerun
list2L <- list.files(path = "Data/RespoFiles/RawO2/20230326/RUN2LIGHT")

data2.1 <- run2fun(list2L[7])
data2.2 <- run2fun(list2L[4])
data2.3 <- run2fun(list2L[3]) %>%
  filter(datetime < ymd_hms("2023-03-26 12:30:00"))
data2.4 <- run2fun(list2L[2])
data2.5 <- run2fun(list2L[6])
data2.7 <- run2fun(list2L[8]) %>% 
  filter(datetime < ymd_hms("2023-03-26 12:47:00"))
data2.8 <- run2fun(list2L[9])
data2.9 <- run2fun(list2L[1]) %>% 
  filter(datetime > ymd_hms("2023-03-26 12:26:00"))

# write files
# write_csv(data2.1, paste0(wpath,"/20230326/RUN2LIGHT/",list2L[7]))
# write_csv(data2.2, paste0(wpath,"/20230326/RUN2LIGHT/",list2L[4]))
# write_csv(data2.3, paste0(wpath,"/20230326/RUN2LIGHT/",list2L[3]))
# write_csv(data2.4, paste0(wpath,"/20230326/RUN2LIGHT/",list2L[2]))
# write_csv(data2.5, paste0(wpath,"/20230326/RUN2LIGHT/",list2L[6]))
# write_csv(data2.7, paste0(wpath,"/20230326/RUN2LIGHT/",list2L[8]))
# write_csv(data2.8, paste0(wpath,"/20230326/RUN2LIGHT/",list2L[9]))
# write_csv(data2.9, paste0(wpath,"/20230326/RUN2LIGHT/",list2L[1]))


### RUN 2 DARK
path2d <- "Data/RespoFiles/RawO2/20230326/RUN2DARK"
list2d <- list.files(path = path2d)

data2.1d<- readfilefun(path2d, list2d[7]) %>%  # a little wonky until 14:26
  filter(datetime > ymd_hms("2023-03-26 14:27:00"))
data2.5d <- readfilefun(path2d, list2d[6]) %>% # keep beginning signal
  filter(between(datetime, ymd_hms("2023-03-26 14:09:20"), ymd_hms("2023-03-26 14:22:58")))
data2.7d<- readfilefun(path2d, list2d[8]) %>%  # cut off end points, not a clean curve...
  filter(datetime < ymd_hms("2023-03-26 14:47:10"))
data2.8d<- readfilefun(path2d, list2d[9]) %>%  # blank chamber a little wonky until 14:26
  filter(datetime > ymd_hms("2023-03-26 14:27:00"))
data2.9d <- readfilefun(path2d, list2d[1]) %>% 
  filter(between(datetime, ymd_hms("2023-03-26 14:48:56"), ymd_hms("2023-03-26 15:02:30")))


# write files
# write_csv(data2.1d, paste0(wpath,"/20230326/RUN2DARK/",list2d[7]))
# write_csv(data2.5d, paste0(wpath,"/20230326/RUN2DARK/",list2d[6]))
# write_csv(data2.7d, paste0(wpath,"/20230326/RUN2DARK/",list2d[8]))
# write_csv(data2.8d, paste0(wpath,"/20230326/RUN2DARK/",list2d[9]))
# write_csv(data2.9d, paste0(wpath,"/20230326/RUN2DARK/",list2d[1]))

# unmodified files
rawfilefun(path2d, paste0(wpath,"/20230326/RUN2DARK"),list2d[2])
rawfilefun(path2d, paste0(wpath,"/20230326/RUN2DARK"),list2d[3])
rawfilefun(path2d, paste0(wpath,"/20230326/RUN2DARK"),list2d[4])


### RUN 3 LIGHT
path3L <- "Data/RespoFiles/RawO2/20230326/RUN3LIGHT"
list3l <- list.files(path = path3L)

data3.1 <- readfilefun(path3L, list3l[9])  # blank, steep o2 drop then levels off ~ 17:40
  #filter(datetime > ymd_hms("2023-03-26 17:15:00"))  # not sure which segment to keep. talk to Nyssa?
data3.2 <- readfilefun(path3L, list3l[2]) %>%  # levels off
  filter(datetime < ymd_hms("2023-03-26 17:30:00"))
data3.3 <- readfilefun(path3L, list3l[3]) %>%  # cut off beginning ~10min and leveling off at the end
  filter(between(datetime, ymd_hms("2023-03-26 17:14:00"), ymd_hms("2023-03-26 17:52:00"))) 
data3.5 <- readfilefun(path3L, list3l[7]) %>%  # broken seal partway through run, ch5
  filter(datetime < ymd_hms("2023-03-26 17:21:0"))  #| datetime > ymd_hms("2023-03-26 17:50:00")) # keep the clean curve at the beginning
data3.7 <- readfilefun(path3L, list3l[4]) %>%  # clean curve at start, then starts to level off
  filter(datetime < ymd_hms("2023-03-26 17:40:00"))



# write files
# write_csv(data3.1, paste0(wpath,"/20230326/RUN3LIGHT/",list3l[9]))
# write_csv(data3.2, paste0(wpath,"/20230326/RUN3LIGHT/",list3l[2]))
# write_csv(data3.3, paste0(wpath,"/20230326/RUN3LIGHT/",list3l[3]))
# write_csv(data3.5, paste0(wpath,"/20230326/RUN3LIGHT/",list3l[7]))
# write_csv(data3.7, paste0(wpath,"/20230326/RUN3LIGHT/",list3l[4]))

# unmodified files
rawfilefun(path3L, paste0(wpath,"/20230326/RUN3LIGHT"),list3l[1])
rawfilefun(path3L, paste0(wpath,"/20230326/RUN3LIGHT"),list3l[5])
rawfilefun(path3L, paste0(wpath,"/20230326/RUN3LIGHT"),list3l[6])
rawfilefun(path3L, paste0(wpath,"/20230326/RUN3LIGHT"),list3l[8])

  
### RUN 3 DARK
# all is well!

path3d <- "Data/RespoFiles/RawO2/20230326/RUN3DARK"
list3d <- list.files(path = path3d)
# write unmodified files
rawfilefun(path3d, paste0(wpath,"/20230326/RUN3DARK"),list3d[1])
rawfilefun(path3d, paste0(wpath,"/20230326/RUN3DARK"),list3d[2])
rawfilefun(path3d, paste0(wpath,"/20230326/RUN3DARK"),list3d[3])
rawfilefun(path3d, paste0(wpath,"/20230326/RUN3DARK"),list3d[4])
rawfilefun(path3d, paste0(wpath,"/20230326/RUN3DARK"),list3d[5])
rawfilefun(path3d, paste0(wpath,"/20230326/RUN3DARK"),list3d[6])
rawfilefun(path3d, paste0(wpath,"/20230326/RUN3DARK"),list3d[7])
rawfilefun(path3d, paste0(wpath,"/20230326/RUN3DARK"),list3d[8])
rawfilefun(path3d, paste0(wpath,"/20230326/RUN3DARK"),list3d[9])


### RUN 4 LIGHT
path4L <- "Data/RespoFiles/RawO2/20230327/RUN4LIGHT"
list4l <- list.files(path = path4L)

data4.7 <- readfilefun(path4L, list4l[3]) %>%  # get rid of the very beginning when seal was breaking and then end when seal was weak again
  filter(between(datetime, ymd_hms("2023-03-27 11:44:00"), ymd_hms("2023-03-27 12:15:00")))

# write files
#write_csv(data4.7, paste0(wpath,"/20230327/RUN4LIGHT/",list4l[3]))

# unmodified files
rawfilefun(path4L, paste0(wpath,"/20230327/RUN4LIGHT"),list4l[1])
rawfilefun(path4L, paste0(wpath,"/20230327/RUN4LIGHT"),list4l[2])
rawfilefun(path4L, paste0(wpath,"/20230327/RUN4LIGHT"),list4l[4])
rawfilefun(path4L, paste0(wpath,"/20230327/RUN4LIGHT"),list4l[5])
rawfilefun(path4L, paste0(wpath,"/20230327/RUN4LIGHT"),list4l[6])
rawfilefun(path4L, paste0(wpath,"/20230327/RUN4LIGHT"),list4l[7])
rawfilefun(path4L, paste0(wpath,"/20230327/RUN4LIGHT"),list4l[8])
rawfilefun(path4L, paste0(wpath,"/20230327/RUN4LIGHT"),list4l[9])


### RUN 4 DARK
path4d <- "Data/RespoFiles/RawO2/20230327/RUN4DARK"
list4d <- list.files(path = path4d)

data4.6d <- readfilefun(path4d, list4d[2]) %>%  # need to remove middle seal break, ch6
  filter(datetime < ymd_hms("2023-03-27 13:49:20"))  # | datetime > ymd_hms("2023-03-27 13:59:00")) # keep the clean curve at the start
data4.7d <- readfilefun(path4d, list4d[3]) %>%  # need to remove middle seal weakening, ch7
  filter(datetime < ymd_hms("2023-03-27 13:53:00"))  # | datetime > ymd_hms("2023-03-27 13:59:00")) %>% # keep the clean curve at the start

# write files
# write_csv(data4.6d, paste0(wpath,"/20230327/RUN4DARK/",list4d[2]))
# write_csv(data4.7d, paste0(wpath,"/20230327/RUN4DARK/",list4d[3]))

# unmodified files
rawfilefun(path4d, paste0(wpath,"/20230327/RUN4DARK"),list4d[1])
rawfilefun(path4d, paste0(wpath,"/20230327/RUN4DARK"),list4d[4])
rawfilefun(path4d, paste0(wpath,"/20230327/RUN4DARK"),list4d[5])
rawfilefun(path4d, paste0(wpath,"/20230327/RUN4DARK"),list4d[6])
rawfilefun(path4d, paste0(wpath,"/20230327/RUN4DARK"),list4d[7])
rawfilefun(path4d, paste0(wpath,"/20230327/RUN4DARK"),list4d[8])
rawfilefun(path4d, paste0(wpath,"/20230327/RUN4DARK"),list4d[9])



### RUN 5 LIGHT
# no channel 8 in this run
path5L <- "Data/RespoFiles/RawO2/20230327/RUN5LIGHT"
list5l <- list.files(path = path5L)

data5.6 <- readfilefun(path5L, list5l[2]) %>%  # levels off near end
  filter(datetime < ymd_hms("2023-03-27 18:22:00"))
data5.1 <- readfilefun(path5L, list5l[3]) %>%  # remove random point and leveling at end
  filter(Value < 350) %>% 
  filter(datetime < ymd_hms("2023-03-27 18:26:00"))
data5.2 <- readfilefun(path5L, list5l[8]) %>%  # remove random point, ch2, blank showed dec then inc signals
  filter(Value < 250) 


# write files
# write_csv(data5.6, paste0(wpath,"/20230327/RUN5LIGHT/",list5l[2]))
# write_csv(data5.1, paste0(wpath,"/20230327/RUN5LIGHT/",list5l[3]))
# write_csv(data5.2, paste0(wpath,"/20230327/RUN5LIGHT/",list5l[8]))

# unmodified files
rawfilefun(path5L, paste0(wpath,"/20230327/RUN5LIGHT"),list5l[1])
rawfilefun(path5L, paste0(wpath,"/20230327/RUN5LIGHT"),list5l[4])
rawfilefun(path5L, paste0(wpath,"/20230327/RUN5LIGHT"),list5l[5])
rawfilefun(path5L, paste0(wpath,"/20230327/RUN5LIGHT"),list5l[6])
rawfilefun(path5L, paste0(wpath,"/20230327/RUN5LIGHT"),list5l[7])



### RUN 5 DARK
# no channel 8 in this run
path5d <- "Data/RespoFiles/RawO2/20230327/RUN5DARK"
list5d <- list.files(path = path5d)

data5.4d <- readfilefun(path5d, list5d[1]) %>% # remove beginning before clear signal
  filter(datetime > ymd_hms("2023-03-27 20:05:00"))
data5.6d <- readfilefun(path5d, list5d[2])%>% # remove beginning before clear signal
  filter(datetime > ymd_hms("2023-03-27 20:05:00"))
data5.3d <- readfilefun(path5d, list5d[4]) %>% # remove beginning before clearer signal
  filter(datetime > ymd_hms("2023-03-27 20:00:00")) 
data5.2d <- readfilefun(path5d, list5d[8])  %>% # remove beginning before clearer signal
  filter(datetime > ymd_hms("2023-03-27 20:00:00"))

# write files
# write_csv(data5.4d, paste0(wpath,"/20230327/RUN5DARK/",list5d[1]))
# write_csv(data5.6d, paste0(wpath,"/20230327/RUN5DARK/",list5d[2]))
# write_csv(data5.3d, paste0(wpath,"/20230327/RUN5DARK/",list5d[4]))
# write_csv(data5.2d, paste0(wpath,"/20230327/RUN5DARK/",list5d[8]))

# unmodified files
rawfilefun(path5d, paste0(wpath,"/20230327/RUN5DARK"),list5d[3])
rawfilefun(path5d, paste0(wpath,"/20230327/RUN5DARK"),list5d[5])
rawfilefun(path5d, paste0(wpath,"/20230327/RUN5DARK"),list5d[6])
rawfilefun(path5d, paste0(wpath,"/20230327/RUN5DARK"),list5d[7])


### RUN 6 LIGHT
# no channel 7 in this run due to light malfunction
path6L <- "Data/RespoFiles/RawO2/20230328/RUN6LIGHT"
list6l <- list.files(path = path6L)

data6.3 <- readfilefun(path6L, list6l[5]) %>%  # starts to level off
  filter(datetime < ymd_hms("2023-03-28 15:20:00"))
data6.2 <- readfilefun(path6L, list6l[6]) %>%  # multiple seal breaks, ch2
  filter(between(datetime, ymd_hms("2023-03-28 14:57:00"), ymd_hms("2023-03-28 15:27:00"))) # longest stretch is between two seal breaks mid-run
data6.1 <- readfilefun(path6L, list6l[7]) %>%  # small blip near the end
  filter(datetime < ymd_hms("2023-03-28 15:33:00"))

# write files
# write_csv(data6.3, paste0(wpath,"/20230328/RUN6LIGHT/",list6l[5]))
# write_csv(data6.2, paste0(wpath,"/20230328/RUN6LIGHT/",list6l[6]))
# write_csv(data6.1, paste0(wpath,"/20230328/RUN6LIGHT/",list6l[7]))

# unmodified files
rawfilefun(path6L, paste0(wpath,"/20230328/RUN6LIGHT"),list6l[1])
rawfilefun(path6L, paste0(wpath,"/20230328/RUN6LIGHT"),list6l[2])
rawfilefun(path6L, paste0(wpath,"/20230328/RUN6LIGHT"),list6l[3])
rawfilefun(path6L, paste0(wpath,"/20230328/RUN6LIGHT"),list6l[4])
rawfilefun(path6L, paste0(wpath,"/20230328/RUN6LIGHT"),list6l[8])


### RUN 6 DARK
# no channel 7 in this run due to light malfunction in Light Run
# all is well!
path6d <- "Data/RespoFiles/RawO2/20230328/RUN6DARK"
list6d <- list.files(path = path6d)

# unmodified files
rawfilefun(path6d, paste0(wpath,"/20230328/RUN6DARK"),list6d[1])
rawfilefun(path6d, paste0(wpath,"/20230328/RUN6DARK"),list6d[2])
rawfilefun(path6d, paste0(wpath,"/20230328/RUN6DARK"),list6d[3])
rawfilefun(path6d, paste0(wpath,"/20230328/RUN6DARK"),list6d[4])
rawfilefun(path6d, paste0(wpath,"/20230328/RUN6DARK"),list6d[5])
rawfilefun(path6d, paste0(wpath,"/20230328/RUN6DARK"),list6d[6])
rawfilefun(path6d, paste0(wpath,"/20230328/RUN6DARK"),list6d[7])
rawfilefun(path6d, paste0(wpath,"/20230328/RUN6DARK"),list6d[8])


