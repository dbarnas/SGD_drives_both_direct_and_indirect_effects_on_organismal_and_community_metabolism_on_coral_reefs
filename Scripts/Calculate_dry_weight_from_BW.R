#### Calculate Species Weights from Buoyant Weight (Davies 1989)
### Created by Danielle Barnas

## dry weight of object = weight in water + ((weight in air * Density of water) / Density of object)

## Density of aragonite = 2.93 g/cm-3 (Jokiel et al. 1978)
## Density of CCA = 2.71 g cm-3 (Comeau et al. 2014) "converted to dry weight"
## SA: aluminum foil technique (Marsh 1970)
## avg sea water density = 1.023 g cm-3



#############################
### LOAD LIBRARIES
#############################
library(tidyverse)
library(here)
library(seacarb)



#############################
### READ IN DATA
#############################
bw <- read_csv(here("Data","Growth","buoyant_weights.csv"))


#############################
### CALCULATE DRY WEIGHT
#############################
#### Using function rho from seacarb package
## rho(S = 35, T = 25, P = 0)

bwd <- bw %>% 
  select(date:Salinity) %>% 
  separate(SpeciesID, into = c("ID", "Rep"), sep = " ", remove = FALSE) %>% 
  filter(ID != "DI" &
           ID != "Air" &
           ID != "SW") %>% 
  select(-c(ID, Rep)) %>% 
  mutate(sw_dens = rho(S = Salinity, T = Temp.c, P = 0), # calculate density of seawater
         sw_dens = sw_dens * 0.001) %>%  # convert from kg cm-3 to g cm-3
  drop_na(sw_dens) %>%
  mutate(n = 1:nrow(.))

# isolate LK for skeletal density constant
nLK <- grep(bwd$SpeciesID, pattern = "LK")
nLK <- bwd %>%
  filter(n %in% nLK) %>% 
  mutate(Sp = "LK")

# rejoin LK data to include Sp column and add skeletal densities
bwd <- bwd %>% 
  anti_join(nLK) %>% 
  mutate(Sp = "Coral") %>% 
  rbind(nLK) %>%  # deals with possible join issue
  select(-n) %>% 
  mutate(skel_dens = if_else(Sp == "LK", 2.71, 2.93)) %>% # add respective skeletal densities
  mutate(dry_weight.g = BW.g / (1 - (sw_dens/skel_dens))) %>% 
  select(-Sp)

# dry weight of object = weight in water / (1 - (Density of water / Density of object)) (Jokiel et al. 1978)


# check if any drift in measurements pre- and post-deployment
bwd_standard <- bw %>% 
  select(date:Salinity) %>% 
  separate(SpeciesID, into = c("ID", "Rep"), sep = " ", remove = FALSE) %>% 
  filter(ID == "DI" |
           ID == "Air" |
           ID == "SW") %>% 
  select(-c(AT, ET))
# no significant difference between pre- and post-deployment
bwd_standard %>% 
  group_by(pre_post, ID) %>% 
  summarise(bw_mean = mean(BW.g, na.rm = TRUE),
            bw_se = plotrix::std.error(BW.g, na.rm = TRUE))

#############################
### SAVE FILE
#############################

#write_csv(bwd, here("Data", "Growth", "Skeletal_Dry_Weight_Calc.csv"))



  
