###########################################################
### GROWTH
###########################################################


#############################
### LOAD LIBRARIES
#############################
library(tidyverse)
library(here)
library(patchwork)
library(stats)
library(ggrepel)
library(lme4)
library(lmerTest) #Need this to get anova results from lmer
library(agricolae) # HSD.test()
library(car)
library(plotrix)



#############################
### READ IN DATA
#############################
# bw <- read_csv(here("Data","Growth","buoyant_weights.csv"))
# ww <- read_csv(here("Data", "Growth", "wet_weights_disp.csv"))
# disp <- read_csv(here("Data", "Growth", "water_displacement.csv"))
allspecies <- read_csv(here("Data", "RespoFiles", "SpeciesMetadata_calculated_perday.csv")) %>% 
  drop_na(delWeight.mg_norm_day)



#############################
### CLEAN DATA
#############################

### GROWTH ACROSS TREATMENTS

#### remove unpaired values
anti_pair <- allspecies %>% 
  group_by(Sp,SpRep, AT) %>% 
  count(SpRep) %>% 
  filter(n == 1) %>% 
  select(-n)
species <- allspecies %>% 
  anti_join(anti_pair)

# for TO length values
to_anti_pair <- species %>% 
  filter(Sp == "TO") %>% 
  drop_na(pLength) %>% 
  group_by(Sp,SpRep, AT) %>% 
  count(SpRep) %>% 
  filter(n == 1) %>% 
  select(Sp,SpRep, AT)
species <- species %>% 
  anti_join(to_anti_pair) %>% 
  filter(SpeciesID != "PR3HL") # must have mis-weighed because no growth shown.


#### to use even unpaired orgs (as random effects)
species <- allspecies %>% 
  filter(Sp != "LC") %>% # full species atrophe
  filter(Sp != "TO") # should look into which kept their full structure, remove all for now
# remove outliers due to breakage or non-treatment tissue loss
broken <- species %>% 
  filter(Sp == "PR" & pWeight < 2 |
         Sp == "GS" & pWeight < -5 | 
         Sp == "HO" & pWeight < 0)
species <- species %>% 
  anti_join(broken)

# make ET a factor with levels Low then High
species <- species %>% 
  mutate(ET = factor(ET, levels = c("LOW", "HIGH")))


#############################
### SAVE FILE
#############################

write_csv(species, here("Data", "Growth", "All_Weight_pChange.csv"))





#############################
### VISUALIZE
#############################

#### percent change in weight
weightPlot <- species %>% 
  unite(Sp, AT, col = "Sp_AT", remove = F, sep = "_") %>% 
  unite(Sp, SpRep, AT, col = "SpRep_AT", remove = F, sep = "_") %>% 
  ggplot(aes(x = ET, y = pWeight, fill = ET)) +
  geom_boxplot(alpha = 0.6, show.legend = FALSE) +
  geom_point(position = position_jitterdodge(),
             shape = 21,
             size = 3,
             color = "black",
             aes(fill = ET),
             show.legend = FALSE) +
  facet_wrap(~Sp, scales = "free_y") +
  theme_bw() +
  geom_line(aes(group = SpRep_AT, color = AT),
            size = 1) +
  scale_fill_manual(values = c("grey", "grey4")) +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15)) +
  labs(x = "SGD Environmental Treatment",
       y = "% Weight change",
       color = "Assemblage \nTreatment")
#### weight/cm2/day
growthPlot <- species %>% 
  unite(Sp, AT, col = "Sp_AT", remove = F, sep = "_") %>% 
  unite(Sp, SpRep, AT, col = "SpRep_AT", remove = F, sep = "_") %>% 
  ggplot(aes(x = ET, y = delWeight.mg_norm_day, fill = ET)) +
  geom_boxplot(alpha = 0.6, show.legend = FALSE) +
  geom_point(position = position_jitterdodge(),
             shape = 21,
             size = 3,
             color = "black",
             aes(fill = ET),
             show.legend = FALSE) +
  facet_wrap(~Sp, scales = "free_y") +
  theme_bw() +
  geom_line(aes(group = SpRep_AT, color = AT),
            size = 1) +
  scale_fill_manual(values = c("grey", "grey4")) +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15)) +
  labs(x = "SGD Environmental Treatment",
       y = "Growth",
       color = "Assemblage \nTreatment")


species <- species %>% 
  mutate(PartSp = if_else(Sp == "DN", "D. nummiforme",
                  if_else(Sp == "GS", "Porifera unk", 
                  if_else(Sp == "HO", "H. opuntia", 
                  if_else(Sp == "LK", "L. kotschyanum",
                  if_else(Sp == "ME", "M. grisea",
                  if_else(Sp == "PA", "P. acuta",
                  if_else(Sp == "PR", "P. rus",
                  if_else(Sp == "VF", "V. fastigiata", "")))))))))

meanspecies <- species %>% 
  unite(Sp, AT, col = "Sp_AT", remove = F, sep = "_") %>% 
  unite(Sp, SpRep, AT, col = "SpRep_AT", remove = F, sep = "_") %>% 
  group_by(Sp,PartSp, ET) %>% 
  summarise(meanVal = mean(delWeight.mg_norm_day),
            sd = sd(delWeight.mg_norm_day),
            se = std.error(delWeight.mg_norm_day))
  # summarise(meanVal = mean(pWeight),
  #           sd = sd(pWeight),
  #           se = std.error(pWeight))
meanspecies %>% 
  ggplot(aes(x = ET, y = meanVal, fill = ET)) +
  geom_point(data = meanspecies,
             shape = 21,
             size = 2,
             color = "black",
             aes(fill = ET),
             show.legend = FALSE) +
  geom_errorbar(data = meanspecies,
                aes(ymin = meanVal-se, ymax = meanVal+se), width = 0.2) +
  geom_point(data = species,
             position = position_jitterdodge(),
             shape = 21,
             size = 2,
             color = "black",
             aes(x = ET, y = delWeight.mg_norm_day,
             # aes(x = ET, y = pWeight,
                 fill = ET),
             show.legend = FALSE,
             alpha = 0.2) +
  facet_wrap(~PartSp, scales = "free_y") +
  theme_bw() +
  #scale_fill_manual(values = c("black", "black")) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "italic"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15)) +
  labs(x = "SGD Exposure Treatment",
       y = "growth",
       color = "Assemblage \nTreatment")

anova(lmer(data = species %>% filter(Sp == "PR"), delWeight.mg_norm_day ~ ET + (1|SpRep)))
anova(lmer(data = species %>% filter(Sp == "VF"), delWeight.mg_norm_day ~ ET + (1|SpRep)))
anova(lmer(data = species %>% filter(Sp == "HO"), delWeight.mg_norm_day ~ ET + (1|SpRep)))
