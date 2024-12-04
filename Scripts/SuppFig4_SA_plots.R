#### Supplemental Figure 3: Surface area standard regression
#### Created by Danielle Barnas


##########################################################
### Supplemental Figure 3
##########################################################

#############################
### LOAD LIBRARIES
#############################
library(tidyverse)
library(here)
library(patchwork)


#############################
### READ IN DATA
#############################
foilSA <- read_csv(here("Data", "Growth", "foil_SA_curve.csv"))
waxSA <- read_csv(here("Data", "Growth", "wax_dip_curve.csv"))


#############################
### CREATE PLOT
#############################
wax.plot <- waxSA %>% 
  ggplot(aes(x = dowel_SA_cm2, y = wax_weight)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = expression("Dowel surface area (cm"^2*")"),
       y = "Wax weight (g)")

foil.plot <- foilSA %>% 
  ggplot(aes(x = dowel_SA_cm2, y = aluminum_weight)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = expression("Dowel surface area (cm"^2*")"),
       y = "Aluminum foil weight (g)")

sa.plots <- wax.plot + foil.plot

summary(lm(data = waxSA, wax_weight ~ dowel_SA_cm2))
summary(lm(data = foilSA, aluminum_weight ~ dowel_SA_cm2))
  
#############################
### SAVE PLOT
#############################
# ggsave(here("Output", "PaperFigures", "SuppFig3.png"), sa.plots, device = "png", height = 5, width = 8)
