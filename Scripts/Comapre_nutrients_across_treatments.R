###########################################################
### NUTRIENTS ACROSS TREATMENTS
###########################################################


#############################
### LOAD LIBRARIES
#############################
library(tidyverse)
library(lubridate)
library(here)
library(car) # test data normality
library(curl) # pull data from url



#############################
### READ IN DATA
#############################
AugChemData<-read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv')) %>% 
  mutate(Season = "Summer2021") %>% 
  filter(Location == "Varari") # study focal site
newchem <- read_csv(here("Data", "Nutrients_experimental_sites.csv"))


#############################
### CLEAN DATA
### 2021 summer water sampling
#############################


### REMOVE UNRELATED DATA POINTS

#### Filter out redundant low tide day sample, first 'low tide' was super high
removeSite1 <- AugChemData %>%
  filter(Season == "Summer2021",
         Tide == "Low",
         Day_Night == "Day",
         Date == ymd("2021-08-06"))

removeSite2 <- AugChemData %>%
  filter(CowTagID %in% c("VSPRING", "Varari_Well"))

#### Anti-join and select only relevant columns
AugChem <- AugChemData %>%
  anti_join(removeSite1) %>%
  anti_join(removeSite2) %>% 
  select(CowTagID,
         Tide,
         Day_Night,
         Salinity:NN_umolL,
         Season)

### SELECT RELEVANT SITES

#### Within August data, select sites V13 and V17 relating to the two 2023 experimental locations
ReducedChemData <- AugChem %>% 
  relocate(Season, .after = Day_Night) %>% 
  filter(CowTagID == "V13" | # low SGD
           CowTagID == "V17" | # high SGD
           CowTagID == "VSEEP") %>% 
  mutate(Site = if_else(CowTagID == "V17", "HighSGD", 
                        if_else(CowTagID == "V13", "LowSGD", "Seep"))) %>% 
  select(Tide:NN_umolL, Site) %>% 
  drop_na() %>%  # don't have some nutrient data at V13 Summer2021 season
  select(-c(Salinity:pH))
  


#############################
### CLEAN DATA
### 2023 spring water sampling
#############################

#### create columns to join with August chem data
NewRedChem <- newchem %>% 
  separate(Sample_ID, into = c("Site", "DN_Tide", "Date", "Time"), sep = "_") %>% 
  mutate(Day_Night = if_else(grepl("D", DN_Tide), "Day", "Night"),
         Tide = if_else(grepl("L", DN_Tide), "Low", 
                        if_else(grepl("H", DN_Tide), "High", "Mid"))) %>% 
  select(-c(SLAB_ID, Ammonia_umolL, Date, Time, DN_Tide)) %>% 
  mutate(Season = "Spring2023")



#### Join sampling data and remove Seep
fullchem <- ReducedChemData %>% 
  full_join(NewRedChem) %>% 
  filter(Site != "Seep")


#############################
### SUMMARIZE DATA
#############################

### Summarize across Day_Night samples
fullSumChem_DN <- fullchem %>% 
  pivot_longer(cols = Phosphate_umolL:NN_umolL, 
               names_to = "Parameters", 
               values_to = "Values") %>% 
  group_by(Site, Season, Tide, Parameters) %>% 
  summarise(Mean = mean(Values, na.rm = TRUE),
            SD = sd(Values, na.rm = TRUE),
            SE = plotrix::std.error(Values, na.rm = TRUE),
            Min = min(Values, na.rm = TRUE),
            Max = max(Values, na.rm = TRUE),
            Range = Max-Min,
            Var = var(Values, na.rm = TRUE),
            CV = SD/Mean)


#############################
### STATISTICS
#############################

### CHECK ASSUMPTIONS
mod1 <- lm(data = fullchem %>% filter(Site != "Seep"), 
           Phosphate_umolL ~ Site)
plot(mod1)
qqp(mod1)
leveneTest(mod1)

mod2 <- lm(data = fullchem %>% filter(Site != "Seep"), 
           NN_umolL ~ Site)
plot(mod2)
qqp(mod2)
leveneTest(mod2)



### STANDARD DEVIATION ~ INTERACTION OF SITE AND TIDE

#### Summer2021, Spring2023
anova(lm(data = fullSumChem_DN %>% filter(Parameters=="NN_umolL"), # significant: tide p=0.046 and interaction p=0.028
         Var ~ Site*Tide))
anova(lm(data = fullSumChem_DN %>% filter(Parameters=="NN_umolL"), # significant interaction p=0.044
         SD ~ Site*Tide))

#### Summer2021, Spring2023
anova(lm(data = fullSumChem_DN %>% filter(Parameters=="Phosphate_umolL"), # significant: tide p=0.017 and nearly interaction p=0.051
         Range ~ Site*Tide))
anova(lm(data = fullSumChem_DN %>% filter(Parameters=="Phosphate_umolL"), # significant: tide p=0.029 and interaction p=0.049
         SD ~ Site*Tide))
         
         
         ### Standard deviation between day and night samples is significant across tides and sites 
         ### for both NN and Phosphate across Summer 2021 and Spring 2023 data
         ### Therefore, there is significant dispersal about the mean for both NN and Phosphate across sites and tides


#############################
### VISUALIZE
#############################
fullchem %>% 
  pivot_longer(cols = Phosphate_umolL:NN_umolL, 
               names_to = "Parameters", 
               values_to = "Values") %>% 
  ggplot(aes(x = Site, y = Values, color = Tide, shape = Day_Night)) +
  geom_point(size = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  facet_wrap(~Parameters, scales = "free_y")
