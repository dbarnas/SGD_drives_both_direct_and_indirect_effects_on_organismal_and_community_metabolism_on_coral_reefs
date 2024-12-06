# CALCULATE NET ECOSYSTEM CALCIFICATION OF ASSEMBLAGES IN SGD TREATMENTS
# Created by Danielle Barnas
# Created on June 3, 2023

#############################
### LOAD LIBRARIES
#############################

library(tidyverse)
library(here)
library(patchwork)
library(PNWColors)
library(plotrix)




#############################
### READ IN DATA
#############################

rawTA <- read_csv(here("Data", "RespoFiles", "RawTA_2023.csv"))
rVolume <- read_csv(here("Data", "RespoFiles", "RespoVolume.csv"))
Rmeta <- read_csv(here("Data", "RespoFiles", "RespoMetadata.csv"))
Ameta <- read_csv(here("Data", "RespoFiles", "AssemblageMetadata_calc.csv"))




#############################
### CLEAN DATA
#############################

# relabel SampleID's in Rmeta
Rmeta <- Rmeta %>% 
  filter(run_block != "RUN1") %>% 
  mutate(SampleID = if_else(SampleID == "A11_H", "A1_H",
                    if_else(SampleID == "A11_L", "A1_L",
                    if_else(SampleID == "A12_H", "A2_H",
                    if_else(SampleID == "A12_L", "A2_L",
                    if_else(SampleID == "B11_H", "B1_H",
                    if_else(SampleID == "B11_L", "B1_L",
                    if_else(SampleID == "B12_H", "B2_H",
                    if_else(SampleID == "B12_L", "B2_L", SampleID)))))))))

# get treatment - chamber_channel joining df
treat_ch <- Rmeta %>% 
  filter(run_block != "RUN1") %>% 
  filter(AT != "BLANK") %>% 
  select(SampleID, run_block, AT, ET, chamber_channel) %>% 
  distinct()
  

# clean rawTA
rawTA <- rawTA %>% 
  filter(run_block != "RUN1") %>% # remove RUN1 from all analyses. assemblages were re-run later
  select(-c(BottleID, Day_Night, SGD, Tide, Time, Salinity))

# bind Initial values as separate column
Init <- rawTA %>% 
  filter(Group == "INITIAL") %>% 
  pivot_wider(names_from = Group, values_from = TA)
myTA <- rawTA %>% 
  filter(Group != "INITIAL") %>% 
  full_join(Init) %>% 
  drop_na(TA) %>%  # remove runs that did not exist (should delete rows from "Initial" Group)
  left_join(treat_ch)

# calculate assemblage volume
myVol <- Ameta %>% 
  separate(SampleID, into = c("Group", NA), sep = "_") %>% 
  select(Group:Volume.ml) %>% 
  mutate(ch.volume = 6000-Volume.ml) # subtract assemblage volume from total chamber volume to get total water volume

# assembalge afdw (for normalization)
myAFDW <- Ameta %>% 
  select(-Volume.ml)
  
myAFDW


#############################
### CALCULATE NEC RATE
#############################


# Blank delta TA
Blank_TA <- myTA %>% 
  filter(Group == "BLANK") %>% 
  mutate(Blank_delTA = INITIAL - TA) %>% 
  select(run_block, Group, Blank_delTA) %>% 
  pivot_wider(names_from = Group, values_from = Blank_delTA) %>% 
  rename(Blank_delTA = BLANK)

# calculate time duration
myMeta <- Rmeta %>% 
  filter(light_dark == "LIGHT") %>% 
  select(Date:chamber_channel) %>% 
  separate(SampleID, into = c("Group", NA, NA), sep = "_") %>% # create joining variable for TA
  mutate(time.sec = as.numeric(stop_time - start_time),
         time.min = time.sec/60,
         time.hr = time.min / 60) %>% 
  select(-c(time.sec, time.min, start_time, stop_time, Date))


# calculate delta TA and correct for blank values
delTA <- myTA %>%
  select(-Date) %>% 
  filter(Group != "BLANK") %>% 
  left_join(Blank_TA) %>% 
  mutate(del_TA = INITIAL - TA,
         del_TA.corr = del_TA - Blank_delTA)

# join TA to respo metadata and calculate NEC
myNEC <- delTA %>% 
  left_join(myMeta) %>% 
  left_join(myVol) %>% # join chamber volumes
  left_join(myAFDW) %>%  # join assemblage AFDW to normalize calcification rates
  
  # using AFDW (biomass) to normalize all orgs
  mutate(NEC.umol.g.hr = (del_TA.corr/2) * 1.023 * (ch.volume/AFDW.g) * (1 / time.hr) * (1/1000)) 


# remove duplicate value from scratched chamber run: B3_H from Run 2 in chamber 6 when probe was removed during run
# reran in Run 5
removeNEC <- myNEC %>% 
  filter(SampleID == "B3_H" & chamber_channel == 6)
myNEC <- myNEC %>% 
  anti_join(removeNEC)

#############################
### VISUALIZE
#############################

my_pal <- pnw_palette(name="Starfish",n=2,type="discrete")

NECplot <- myNEC %>% 
  drop_na() %>% 
  ggplot(aes(x=ET, 
             y=NEC.umol.g.hr,
             color = AT))+
  geom_boxplot(aes(color = AT),outlier.shape = NA)+ 
  geom_jitter(aes(color = AT), position = position_jitterdodge()) +
  theme_bw()+
  #geom_text_repel(aes(label = SampleID)) +
  theme(strip.background = element_rect(fill = "white"))+
  labs(x = "Environmental Treatment (SGD Level)",
       color = "Assemblage \n Treatment",
       y = expression("Rate ("*mu*"mol CaCO"[3]*" g-1 hr-1)")) +
  scale_color_manual(values=my_pal)
NECplot
# ggsave(here("Output", "PaperFigures","NECrate.png"), NECplot, device = "png", width = 6, height = 6)



# percent difference across treatments
myNEC %>% 
  mutate(P_R = "NEC") %>% 
  group_by(AT, ET, P_R) %>% 
  select(-c(run_block,chamber_channel,TA,INITIAL,Blank_delTA:AFDW.g)) %>% 
  summarise(meanPR = mean(NEC.umol.g.hr)
            #sdPR = sd(mmol.gram.hr),
            #sePR = sdPR / sqrt(nrow(.))
  ) %>% 
  ungroup() %>% 
  group_by(P_R, ET) %>% 
  pivot_wider(names_from = AT, values_from = meanPR) %>% 
  mutate(ATdiff = 100*(LOW-HIGH)/HIGH) %>% 
  arrange(P_R)

### Patch NEC with R, GP, NP
RespoR_Normalized_Full <- read_csv(here("Data","RespoFiles","Respo_RNormalized_AllRates.csv"))

my_pal <- pnw_palette(name="Starfish",n=2,type="discrete")

myFullRespo <- RespoR_Normalized_Full %>% 
  rename(Values = mmol.gram.hr) %>%  # just for joining df
  select(SampleID:ET, P_R, Values)
myNEC_NEP <- myNEC %>% 
  mutate(P_R = "NC") %>% 
  rename(Values = NEC.umol.g.hr) %>%  # just for joining df
  select(SampleID:ET, P_R, Values) %>% 
  rbind(myFullRespo) %>% 
  mutate(ET = factor(ET, levels = c("LOW", "HIGH")),
         AT = factor(AT, levels = c("LOW", "HIGH")))




#############################
### SAVE FILE
#############################

#write_csv(myNEC_NEP, here("Data", "RespoFiles", 'All_EcoMet_Rates.csv'))



#############################
### VISUALIZE
#############################

myNEC_NEP <- read_csv(here("Data", "RespoFiles", 'All_EcoMet_Rates.csv'))

meanRates <- myNEC_NEP %>% 
  group_by(AT,ET,P_R) %>% 
  summarise(meanVal = mean(Values),
            sd = sd(Values),
            se = std.error(Values))

RatesPlotFun <- function(myValue, renameVal){
  data <- meanRates %>% filter(P_R == {{myValue}}) %>% mutate(P_R = renameVal)
  data2 <- myNEC_NEP %>% filter(P_R == {{myValue}}) %>% mutate(P_R = renameVal)
  RatesPlot <- data %>% 
    ggplot(aes(x=ET, 
               y=meanVal,
               fill = AT))+
    geom_jitter(data = data2,
                aes(x = ET, 
                    y = Values,
                    fill = AT), 
                shape = 21,
                color = "black",
                size = 1.5,
                alpha = 0.5,
                position = position_jitterdodge(dodge.width = 1)) +
    geom_errorbar(data = data,
                  position = position_dodge(width = 1),
                  aes(ymin = meanVal-se, ymax = meanVal+se), width = 0.2) +
    geom_point(data = data,
               shape = 21,
               size = 3,
               color = "black",
               position = position_dodge(width = 1),
               aes(fill = AT))+ 
    theme_bw()+
    #geom_text_repel(aes(label = SampleID)) +
    theme(strip.background = element_rect(fill = "white"))+
    labs(x = "SGD Exposure Treatment",
         fill = "Assemblage \n Treatment",
         #y =""# "Rate O2 change (mmol O2 g-1 hr-1)"
    ) +
    scale_fill_manual(values=my_pal)  +
    facet_wrap(~ P_R, scales = "fixed") +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = 10),
          strip.text = element_text(size = 10),
          panel.grid = element_blank())
  return(RatesPlot)
}


gpp<-RatesPlotFun(myValue = "GP", renameVal = "Gross photosynthesis") + labs(y = expression("Rate (mmol O"[2]*" g-1 hr-1)"))
npp<-RatesPlotFun(myValue = "NP", renameVal = "Net photosynthesis")  + labs(y = expression("Rate (mmol O"[2]*" g-1 hr-1)"))
rp<-RatesPlotFun(myValue = "R", renameVal = "Respiration")  + labs(y = expression("Rate (mmol O"[2]*" g-1 hr-1)"))
ncp<-RatesPlotFun(myValue = "NC", renameVal = "Net calcification") + labs(y = expression("Rate ("*mu*"mol CaCO"[3]*" g-1 hr-1)"))

EcoFunPlot <- (gpp + npp) / (rp + ncp) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = 'collect')  
  # grid::textGrob('SGD Exposure Treatment')
EcoFunPlot


#############################
### STATISTICAL TEST
#############################

library(agricolae) # HSD.test()

# models
model1 <- lm(data = myNEC, NEC.umol.g.hr ~ AT*ET)
anova(model1)
summary(model1)
HSD.test(model1, "AT", console=TRUE)

summary(lm(data = myNEC_NEP %>% filter(P_R == "R") %>% filter(ET == "HIGH"), Values~AT))
myNEC_NEP %>% 
  mutate(P_R = "NP") %>% 
  group_by(AT, ET, P_R) %>% 
  #select(-c(run_block,chamber_channel,TA,INITIAL,Blank_delTA:AFDW.g)) %>% 
  summarise(meanPR = mean(Values)
            #sdPR = sd(mmol.gram.hr),
            #sePR = sdPR / sqrt(nrow(.))
  ) %>% 
  ungroup() %>% 
  group_by(P_R, ET) %>% 
  pivot_wider(names_from = AT, values_from = meanPR) %>% 
  mutate(ATdiff = 100*(LOW-HIGH)/HIGH) %>% 
  arrange(P_R)


model1 <- lm(data = myNEC, NEC.umol.g.hr ~ AT)


### CHECK ASSUMPTIONS

plot(model1)

library(car)
qqp(model1)

modData <- myNEC_NEP %>% 
  pivot_wider(values_from = Values, names_from = P_R)


mod1 <- lm(data = modData, NC ~ AT*ET)
plot(mod1)
qqp(mod1)
leveneTest(mod1)
lowdat <- modData %>% filter(AT == "LOW")
highdat <- modData %>% filter(AT == "HIGH")
shapiro.test(lowdat$NC)
shapiro.test(highdat$NC)

modData %>% 
  select(NC, AT) %>% 
  mutate(NC_trans = NC*(1/3)) %>% 
  ggplot(aes(x = NC_trans))+
  geom_histogram(aes(fill = AT))

mod2 <- lm(data = modData, NP ~ AT)
plot(mod2)
qqp(mod2)
leveneTest(mod2)

mod3 <- lm(data = modData, GP ~ AT)
plot(mod3)
qqp(mod3)
leveneTest(mod3)

mod4 <- lm(data = modData, R ~ AT*ET)
plot(mod4)
qqp(mod4)
leveneTest(mod4)

