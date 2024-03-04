### Process CT from Treatment Locations
### Supplemental Figure 1: Timeseries of salinity and temperature from SGD experimental locations
### Created by Danielle Barnas
### Created on 10/20/2022
### Modified on 4/21/2023

#############################
#### LOAD LIBRARIES ####
#############################
library(tidyverse)
library(here)
library(lubridate)
#library(mooreasgd)
library(timetk)
library(stats)
library(agricolae)
library(emmeans)
library(PNWColors)
library(patchwork)
library(plotrix)
#devtools::install_github("Dbarnas/mooreasgd")


#############################
#### READ IN DATA ####
#############################
ctLog <- read_csv(here("Data","CTDfiles","Treatment_CTDs","QC_files","Full_CT_2023-09-11.csv"))

wl870a <- read_csv(here("Data","CTDfiles","Treatment_CTDs","Spring2023","Depth","Depth_870_HighSGD_20230302.csv"), skip = 1) %>% 
  rename(date = contains("Date"), Pres.kPa = contains("Abs"), TempInSitu = contains("Temp"), Depth.m = contains("Level")) %>% 
  mutate(LoggerID = "870", Treatment = "High")
wl870b <- read_csv(here("Data","CTDfiles","Treatment_CTDs","Spring2023","Depth","Depth_870_HighSGD_20230330.csv"), skip = 1) %>% 
  rename(date = contains("Date"), Pres.kPa = contains("Abs"), TempInSitu = contains("Temp"), Depth.m = contains("Level")) %>% 
  mutate(LoggerID = "870", Treatment = "High")
wl872a <- read_csv(here("Data","CTDfiles","Treatment_CTDs","Spring2023","Depth","Depth_872_LowSGD_20230302.csv"), skip = 1) %>% 
  rename(date = contains("Date"), Pres.kPa = contains("Abs"), TempInSitu = contains("Temp"), Depth.m = contains("Level")) %>% 
  mutate(LoggerID = "872", Treatment = "Low")
wl872b <- read_csv(here("Data","CTDfiles","Treatment_CTDs","Spring2023","Depth","Depth_872_LowSGD_20230330.csv"), skip = 1) %>% 
  rename(date = contains("Date"), Pres.kPa = contains("Abs"), TempInSitu = contains("Temp"), Depth.m = contains("Level")) %>% 
  mutate(LoggerID = "872", Treatment = "Low")


#############################
#### PROCESS WATER LEVEL DATA ####
#############################
# start and end date (full deployment)
sdate <- ymd_hms("2023-02-05 00:00:00")
edate <- ymd_hms("2023-03-24 06:00:00") 

mypal <- pnw_palette("Starfish", n = 7)[c(1, 4)]

#bind and rename column headings for WL logs
wlfile <- rbind(wl870a, wl870b,
                wl872a, wl872b) %>% 
  mutate(date = mdy_hms(date)) %>% 
  unite(Treatment, LoggerID, col = "wlLogger", sep = "-", remove = F) %>% 
  select(-c('#', TempInSitu, Pres.kPa, LoggerID)) %>% 
  filter(between(date, sdate, edate)) %>% 
  mutate(Depth.m = -Depth.m)


#### PLOT WATER LEVEL AT 2 LOCATIONS ####
wlplot <- wlfile %>%  
  filter(between(date, sdate, edate)) %>% 
  ggplot(aes(x = date,
             y = Depth.m,
             color = wlLogger)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = mypal)+
  labs(color = "Depth Logger") +
  ylim(min = -1.25, max = 0) +
  facet_wrap(~wlLogger, scales = "fixed")


#############################
#### PROCESS CT DATA ####
#############################

# assign treatment category and select most reliable CT loggers from each location

ctfile <- ctLog %>% 
  mutate(Treatment = if_else(LoggerID == 349, "High SGD Exposure", 
                             if_else(LoggerID == 323, "High SGD Exposure",
                                     if_else(LoggerID == 320, "Low SGD Exposure", "Low SGD Exposure")))) %>% 
  filter(between(date, sdate,edate)) %>%
  filter(LoggerID != 332,
         LoggerID != 349) %>% 
  unite(Treatment, LoggerID, col = "LoggerID", sep = "-", remove = F) %>% 
  filter(Salinity_psu > 24) %>% 
  select(-FullLoggerID) %>% 
  mutate(Treatment = factor(Treatment, levels = c("Low SGD Exposure", "High SGD Exposure")))

# CT 323 & 349 High
# CT 320 & 332 Low
# ended up only using 323 and 320 as reliable loggers

# remove anomalous / unreliable data points from Low treatment CT
excludect <- ctfile %>%
  filter(Treatment == "Low SGD Exposure") %>% 
  filter(between(date, ymd_hms("2023-03-03 19:10:00"), ymd_hms("2023-03-03 19:30:00")) |
           between(date, ymd_hms("2023-03-05 00:40:00"), ymd_hms("2023-03-05 01:45:00")) |
           between(date, ymd_hms("2023-03-15 18:35:00"), ymd_hms("2023-03-15 18:35:00")) |
           between(date, ymd_hms("2023-03-17 20:05:00"), ymd_hms("2023-03-17 20:55:00")))

exclude.date <- ctfile %>% 
  filter(between(date, ymd_hms("2023-03-02 09:00:00"), ymd_hms("2023-03-05 06:00:00")))

ctfile <- ctfile %>% 
  anti_join(excludect) %>% 
  anti_join(exclude.date)


#### JOIN CT AND WATER LEVEL DATA ####
fulldata <- full_join(ctfile, wlfile) %>% 
  filter(between(date,sdate,edate))
  

#############################
#### VISUALIZE CTD DATA ####
#############################
ctfile %>% 
  filter(between(date,sdate,edate)) %>% 
  group_by(LoggerID) %>% 
  summarise(meanct = mean(Salinity_psu),
            se = std.error(Salinity_psu),
            minct = min(Salinity_psu),
            maxct = max(Salinity_psu),
            meanT = mean(TempInSitu),
            minT = min(TempInSitu),
            maxT = max(TempInSitu),
            seT = std.error(TempInSitu))


mypal2 <- pnw_palette("Starfish", n = 7)[c(7, 3)]


## salinity plot
ctplot <- ctfile %>%
  filter(between(date,sdate,edate)) %>% 
  ggplot(aes(x = date, y = Salinity_psu, color = Treatment)) + 
  geom_point(size=0.5) + 
  theme_bw() +
  theme(strip.background = element_rect(fill= "white"),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(x = "", y = "Salinity (psu)") + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  facet_wrap(~Treatment, ncol = 2) +
  scale_color_manual(values = c("#0072B2", "#D55E00"))

## temperature plot
tempplot <- ctfile %>%
  filter(between(date,sdate,edate)) %>% 
  ggplot(aes(x = date, y = TempInSitu, color = Treatment)) + 
  geom_point(size = 0.5) + 
  theme_bw() +
  theme(strip.background = element_rect(fill= "white"), 
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12))+
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank()) +
  labs(x = "Date", y = expression("Temperature ("*~degree*C*")")) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  facet_wrap(~Treatment, ncol = 2) +
  scale_color_manual(values = c("#0072B2", "#D55E00"))

## depth plot
d.rename <- c(`High-870` = "High SGD Exposure", `Low-872` = "Low SGD Exposure")

wlplot <- wlfile %>%  
  filter(between(date,sdate,edate)) %>% 
  filter(between(date, sdate, edate)) %>% 
  ggplot(aes(x = date,
             y = Depth.m,
             color = wlLogger)) +
  geom_point(size = 0.5) +
  theme_bw() +
  theme(strip.background = element_rect(fill= "white"), 
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_color_manual(values = mypal2)+
  labs(x = "Date",
       y = "Depth (m)") +
  ylim(min = -1.25, max = 0) +
  facet_wrap(~wlLogger, scales = "fixed", ncol = 2,
             labeller = as_labeller(d.rename))



#############################
#### SUPPLEMENTAL FIGURE 1: SGD ENVIRONMENTAL DATA ####
#############################

salPlot <- ctplot / tempplot + # / wlplot +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = "a")

salPlot




#ggsave(here("Output", "PaperFigures","Supp_Salinity_plot.png"),salPlot, device = "png", width = 8, height = 6)



