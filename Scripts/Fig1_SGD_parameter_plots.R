#### Figure 1 and Supplemental Figure 2: Environmental characteristics of SGD experimental locations
#### Created by Danielle Barnas


##########################################################
### Figure 1C-E
##########################################################

#############################
### LOAD LIBRARIES
#############################
library(tidyverse)
library(here)
library(lubridate)
library(timetk)
library(PNWColors)
library(patchwork)
library(plotrix)
library(stats)
library(agricolae)
library(emmeans)




#############################
### READ IN DATA
#############################
ctLog <- read_csv(here("Data","CTDfiles","Treatment_CTDs","QC_files","Full_CT_2023-09-11.csv"))



#### PROCESS CT DATA ####
# start and end date (full deployment)
sdate <- ymd_hms("2023-02-05 00:00:00")
edate <- ymd_hms("2023-03-24 06:00:00") 



# assign treatment category and select most reliable CT loggers from each location

ctfile <- ctLog %>% 
  mutate(Treatment = if_else(LoggerID == 349, "High SGD Exposure", 
                             if_else(LoggerID == 323, "High SGD Exposure",
                                     if_else(LoggerID == 320, "Low SGD Exposure", "Low SGD Exposure")))) %>% 
  filter(between(date, sdate,edate)) %>%
  filter(LoggerID != 332,
         LoggerID != 349) %>% 
  unite(Treatment, LoggerID, col = "LoggerID", sep = "-", remove = F) %>% 
  filter(Salinity_psu > 24) %>% # excludes data points when air-exposed
  select(-FullLoggerID)

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



#############################
### VISUALIZE CTD AND NUTRIENT DATA 
#############################

mypal2 <- pnw_palette("Starfish", n = 7)[c(7, 3)]

## CALCULATE SIMPLE STATS
sum_ct <- ctfile %>% 
  filter(between(date,sdate,edate)) %>% 
  separate(date, into = c("date", "time"), sep = " ") %>% 
  group_by(Treatment, date) %>% 
  summarise(mean_sal = mean(Salinity_psu, na.rm = TRUE),
            max_sal = max(Salinity_psu),
            min_sal = min(Salinity_psu),
            range_sal = max_sal - min_sal,
            mean_temp = mean(TempInSitu, na.rm = TRUE),
            max_temp = max(TempInSitu),
            min_temp = min(TempInSitu),
            range_temp = max_temp - min_temp)%>% 
  ungroup() %>% 
  mutate(Treatment = if_else(Treatment == "High SGD Exposure", "High SGD", "Low SGD"),
         Treatment = factor(Treatment, levels = c("Low SGD", "High SGD"))) # refactor for graphical order

## average daily salinity range
anova(lm(data = sum_ct,
         range_sal ~ Treatment))

mean_sal_range_plot <- sum_ct %>% 
  group_by(Treatment) %>% 
  summarise(mean_range = mean(range_sal, na.rm = TRUE),
            se = plotrix::std.error(range_sal, na.rm = TRUE)) %>% 
  ggplot(aes(x = Treatment, y = mean_range, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_range - se, ymax = mean_range + se), width = .1) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_line(color = "aliceblue"), # present but barely visible
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12)) +
  labs(y = "Mean daily salinity range (psu)",
       x = "Environmental Exposure Treatment") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(values = c("#0072B2", "#D55E00"))  #"blue", "red"))

## minimum salinity
# points
mean_min_sal_plot <- sum_ct %>% 
  group_by(Treatment) %>% 
  summarise(mean_min_sal = mean(min_sal, na.rm = TRUE),
            se = plotrix::std.error(min_sal, na.rm = TRUE)) %>% 
  ggplot(aes(x = Treatment, y = mean_min_sal, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_min_sal - se, ymax = mean_min_sal + se), width = .1) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_line(color = "aliceblue"), # present but barely visible
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12)) +
  labs(y = "Mean minimum salinity (psu)",
       x = "Environmental Exposure Treatment") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(values = c("#0072B2", "#D55E00"))  #"blue", "red"))


## average temperature
anova(lm(data = sum_ct,
         mean_temp ~ Treatment))

mean_temp_plot <- sum_ct %>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(mean_temp, na.rm = TRUE),
            se = plotrix::std.error(mean_temp, na.rm = TRUE)) %>% 
  ggplot(aes(x = Treatment, y = mean, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .1) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_line(color = "aliceblue"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12)) +
  labs(y = expression("Mean temperature ("*~degree*C*")"),
       x = "Environmental \nExposure Treatment") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(values = c("#0072B2", "#D55E00"))  #"blue", "red"))

mean_min_sal_plot + mean_sal_range_plot + mean_temp_plot


chem <- read_csv(here("Data", "Nutrients_experimental_sites.csv"))

sum.chem <- chem %>% 
  separate(Sample_ID, into = c("Site", "Tide"), sep = "_", remove = TRUE) %>% 
  select(Site, Tide, NN_umolL) %>% 
  filter(Site != "Seep") %>% 
  group_by(Site) %>% 
  summarise(mean = mean(NN_umolL, na.rm = TRUE),
            se = plotrix::std.error(NN_umolL, na.rm = TRUE)) %>% 
  mutate(Site = if_else(Site == "LowSGD", "Low SGD", "High SGD"),
         Site = factor(Site, levels = c("Low SGD", "High SGD")))

mean_nn_plot <- sum.chem %>% 
  ggplot(aes(x = Site, y = mean, color = Site)) +
  geom_point() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .1) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_line(color = "aliceblue"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12)) +
  labs(y = expression("Mean Nitrate+Nitrite ("*mu*"mol L-1)"),
       x = "Environmental Exposure Treatment") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(values = c("#0072B2", "#D55E00"))  #"blue", "red"))


#############################
### PATCHWORK PLOTS
#############################



Fig1plotC <- (mean_min_sal_plot) +
  plot_annotation(tag_levels = list(c('c')),
                  tag_prefix = '(',
                  tag_suffix = ')') +
  theme(plot.tag = element_text(face = 'italic'))

Fig1plotD <- (mean_temp_plot) +
  plot_annotation(tag_levels = list(c('d')),
                  tag_prefix = '(',
                  tag_suffix = ')') +
  theme(plot.tag = element_text(face = 'italic'))

Fig1plotE <- (mean_nn_plot) +
  plot_annotation(tag_levels = list(c('e')),
                  tag_prefix = '(',
                  tag_suffix = ')') +
  theme(plot.tag = element_text(face = 'italic'))


Fig1plotCDE_axis_title <- wrap_elements(plot_spacer()) +
  labs(tag = expression("Environmental Exposure Treatment")) +
  theme(
    plot.tag = element_text(size = rel(1.3)),
    plot.tag.position = "top"
  )



# ggsave(here("Output", "PaperFigures", "Figure_1c.png"), Fig1plotC, device = "png", height = 3.6, width = 2.9)
# ggsave(here("Output", "PaperFigures", "Figure_1d.png"), Fig1plotD, device = "png", height = 3.6, width = 3)
# ggsave(here("Output", "PaperFigures", "Figure_1e.png"), Fig1plotE, device = "png", height = 3.6, width = 2.9)
# ggsave(here("Output", "PaperFigures", "Figure_1cde_axis_title.png"), Fig1plotCDE_axis_title, device = "png", height = 4, width = 9)




##########################################################
### Supplemental Figure 2: Raw nutrient data from SGD experimental locations
##########################################################

#############################
### READ IN DATA
#############################

chem <- read_csv(here("Data", "Nutrients_experimental_sites.csv")) 
ta <- read_csv(here("Data", "Water_Sample_rawTA.csv")) 


#############################
### PROCESS DATA
#############################

chem <- chem %>% 
  separate(Sample_ID, into = c("Site", "Tide"), sep = "_", remove = TRUE) %>% 
  rename(`Nitrate+Nitrite` = NN_umolL, `Phosphate` = Phosphate_umolL, `Silicate` = Silicate_umolL) %>% 
  filter(Site != "Seep") %>% 
  mutate(Site = if_else(Site == "LowSGD", "Low SGD", "High SGD"),
         Site = factor(Site, levels = c("Low SGD", "High SGD"))) %>% 
  select(-c(SLAB_ID, Ammonia_umolL))

ta <- ta %>% 
  filter(Location != "SEEP") %>% 
  rename(Site = Location) %>% 
  mutate(Site = if_else(Site == "HighSGD", "High SGD", "Low SGD")) %>% 
  select(Site, DT, TA) %>% 
  rename(`Total Alkalinity` = TA,
         Tide = DT)


sum.chem2 <- chem %>% 
  full_join(ta) %>% 
  pivot_longer(cols = Phosphate:`Total Alkalinity`, names_to = "Parameters", values_to = "Values")


#############################
### VISUALIZE
#############################

sum.chem <- chem %>% 
  full_join(ta) %>% 
  pivot_longer(cols = Phosphate:`Total Alkalinity`, names_to = "Parameters", values_to = "Values") %>% 
  group_by(Site, Parameters) %>% 
  summarise(mean = mean(Values, na.rm = TRUE),
            se = plotrix::std.error(Values, na.rm = TRUE)) %>% 
  ungroup()

supp.chem.plot.fun <- function(param = "Nitrate+Nitrite"){
  
  data1 <- sum.chem %>% 
    filter(Parameters == param)
  data2 <- sum.chem2 %>% 
    filter(Parameters == param)
  
  supp.chem.plot <- data1 %>% 
    ggplot(aes(x = Site, 
               y = mean,
               fill = Site)) +
    geom_jitter(data = data2,
                aes(x = Site, 
                    y = Values,
                    fill = Site), 
                shape = 21,
                color = "black",
                size = 1.8,
                alpha = 0.2,
                position = position_jitterdodge(dodge.width = 1)) +
    geom_errorbar(data = data1,
                  aes(ymin = mean-se, ymax = mean+se),
                  position = position_dodge(width = 1),
                  width = 0.2) +
    geom_point(data = data1,
               shape = 21,
               size = 3.2,
               color = "black",
               position = position_dodge(width = 1),
               aes(fill = Site)) + 
    theme_bw()+
    #geom_text_repel(aes(label = SampleID)) +
    theme(strip.background = element_rect(fill = "white"))+
    scale_fill_manual(values = c("#0072B2", "#D55E00"))  +
    facet_wrap(~ Parameters, scales = "free") +
    theme(axis.title.y = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 10),
          strip.text = element_text(size = 12),
          panel.grid = element_blank())
  
  return(supp.chem.plot)
}


#############################
### PATCHWORK PLOTS
#############################

suppa <- supp.chem.plot.fun(param = "Nitrate+Nitrite") +
  labs(fill = "Experimental \nTreatment",
       y = expression("N+N ("*mu*"mol L-1)"))

suppb <- supp.chem.plot.fun(param = "Phosphate") +
  labs(fill = "Experimental \nTreatment",
       y = expression("PO"[4]*" ("*mu*"mol L-1)"))

suppc <- supp.chem.plot.fun(param = "Silicate") +
  labs(fill = "Experimental \nTreatment",
       y = expression("SiO"[4]*" ("*mu*"mol L-1)"))

suppd <- supp.chem.plot.fun(param = "Total Alkalinity") +
  labs(fill = "Experimental \nTreatment",
       y = expression("TA ("*mu*"mol kg-1)"))

supp.plot2.full <- ((suppa + suppb) / (suppc + suppd)) + 
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") +
  theme(plot.tag = element_text(face = 'italic')) +
  plot_layout(guides = 'collect')


supp.plot2.full


# ggsave(here("Output","PaperFigures","SuppFig2_nutrients.png"), supp.plot2.full, device = "png", width = 8, height = 5)