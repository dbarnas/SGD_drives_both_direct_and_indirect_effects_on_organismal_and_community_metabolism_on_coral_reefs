#### Calculate assemblage bio data for respirometry: volume and total biomass
#### Created by Danielle Barnas
#### Created on 9/12/2022

#############################
### LOAD LIBRARIES
#############################

library(tidyverse)
library(here)


#############################
### READ IN DATA
#############################

species <- read_csv(here("Data", "RespoFiles", "SpeciesMetadata.csv"))
#bw_weight <- read_csv(here("Data", "Growth", "buoyant_weights.csv"))
bw_weight <- read_csv(here("Data", "Growth", "Skeletal_Dry_Weight_Calc.csv")) # calculated dry weights
ww_weight <- read_csv(here("Data", "Growth", "wet_weights_disp.csv"))
coralSA <- read_csv(here("Data", "Growth","wax_dip_corals.csv"))
ccaSA <- read_csv(here("Data", "Growth","foil_SA_CCA.csv"))
respo_vol <- read_csv(here("Data","RespoFiles", "RespoVolume.csv"))
afdw <- read_csv(here("Data","RespoFiles","AFDW.csv"))


#############################
#### Calculate Species Biometrics ####
#############################

# remove dummy columns
species <- species %>% 
  select(-c(Weight.g:AFDW.g))

metaSpecies <- respo_vol %>% 
  full_join(species) %>% 
  select(SampleID, SpeciesID, AT, ET)

# isolate species vector to only join species used for respo
respo_species <- species %>% 
  full_join(respo_vol) %>% 
  filter(Respo == "yes") %>% 
  mutate(RespoVolume.ml = postVol - preVol) %>% 
  select(-c(postVol, preVol))

# clean surface area df
coralSA <- coralSA %>% 
  mutate(SA_cm2 = 38.6 * Wax_weight + 2.46) %>% # calculate outside of excel
  select(SpeciesID, SA_cm2)
ccaSA <- ccaSA %>% 
  mutate(SA_cm2 = 241 * foil_weight - 2.48) %>% # calculate outside of excel
  select(SpeciesID, SA_cm2)
# join
surfacearea <- rbind(coralSA, ccaSA)


# bring in dates of bw and ww to get growth per time
bweight.dates <- bw_weight %>% 
  select(date, pre_post, SpeciesID)
wweight.dates <- ww_weight %>% 
  select(date, pre_post, SpeciesID)

# calculate buoyant weights
# normalize the change in weight to the initial weight
bw <- bw_weight %>% 
  left_join(species) %>% 
  select(SpeciesID, Sp, SpRep, AT, ET, pre_post, dry_weight.g) %>% 
  drop_na(Sp) %>% # remove all calibration values, leaving species
  separate(col = SpeciesID, into = c("SpeciesID", NA), sep = "_") %>% 
  group_by(SpeciesID, pre_post) %>% 
  mutate(dry_weight.g = mean(dry_weight.g)) %>% # for replicated measurements, take average weight
  ungroup() %>% 
  distinct() %>% # remove duplicated rows from averaging
  pivot_wider(names_from = pre_post, values_from = dry_weight.g) %>% 
  drop_na(post) %>% # remove dead organisms
  mutate(delWeight.g = post - pre,
         delWeight.g_biomnorm = delWeight.g/pre, # normalize to biomass
         pWeight = delWeight.g/pre*100) %>%  # calculate weight difference over soak period
  select(SpeciesID, Sp, SpRep, AT, ET, delWeight.g, delWeight.g_biomnorm, pWeight) %>%  # prep df for joining to species df
  filter(SpeciesID != "PR3HL") # mis-weighed


# normalize weight to surface area and time
# bw <- full_join(bw, surfacearea) %>%
#   mutate(delWeight.g_norm = delWeight.g / SA_cm2) %>%
#   select(-SA_cm2)


# standardize to days between measurements
bw <- bw %>% 
  left_join(bweight.dates) %>% 
  mutate(date = mdy(date)) %>% 
  mutate(numdays = yday(date)) %>% 
  distinct() %>%  # two coral values weirdly duplicated
  pivot_wider(names_from = pre_post, values_from = numdays, id_cols = SpeciesID) %>% 
  mutate(DaysInSitu = post-pre) %>% 
  select(SpeciesID, DaysInSitu) %>% 
  full_join(bw) %>% 
  mutate(delWeight.mg_biomnorm_day = delWeight.g_biomnorm / DaysInSitu*1000) # normalize to days in field/days between measurements



# calculate wet weights
ww <- ww_weight %>%
  left_join(species) %>%
  left_join(afdw) %>% 
  # mutate(Dry.g = DrywTin.g - Tin.g,
  #        AFDW.g = Dry.g - (AFDWwTin.g-Tin.g)) %>% 
  select(SpeciesID, Sp, SpRep, AT, ET, pre_post, WW.g) %>%
  pivot_wider(names_from = pre_post, values_from = WW.g) %>%
  drop_na(post) %>% # for any unprocessed or lost data # lost 39 organisms
  mutate(delWeight.g = post - pre,
         delWeight.g_biomnorm = delWeight.g/pre, # normalize to biomass
         pWeight = delWeight.g/pre*100) %>%  # calculate weight difference over soak period
  select(SpeciesID, Sp, SpRep, AT, ET, delWeight.g, delWeight.g_biomnorm, pWeight)  
  
# standardize to days between measurements
ww <- ww %>% 
  left_join(wweight.dates) %>% 
  mutate(date = mdy(date)) %>% 
  mutate(numdays = yday(date)) %>% 
  distinct() %>%  # two coral values weirdly duplicated
  pivot_wider(names_from = pre_post, values_from = numdays, id_cols = SpeciesID) %>% 
  mutate(DaysInSitu = post-pre) %>% 
  select(SpeciesID, DaysInSitu) %>% 
  full_join(ww) %>% 
  mutate(delWeight.mg_biomnorm_day = delWeight.g_biomnorm / DaysInSitu*1000) # normalize to days in field/days between measurements


# calculate TO lengths
length <- ww_weight %>% 
  left_join(species) %>% 
  select(SpeciesID, Sp, SpRep, AT, ET, pre_post, TopLength.cm) %>%
  pivot_wider(names_from = pre_post, values_from = TopLength.cm) %>%
  drop_na(post) %>% 
  mutate(delTopLength.cm = post - pre,
         pLength = delTopLength.cm/pre*100) %>% 
  select(SpeciesID, Sp, SpRep, AT, ET, delTopLength.cm, pLength)

# calculate individual volumes from displacements
# high priority input data
volume <- ww_weight %>% 
  left_join(species) %>%  
  mutate(Volume.ml  = postVol - preVol) %>% 
  select(SpeciesID, Sp, SpRep, AT, ET, pre_post, Volume.ml) %>% 
  drop_na(Volume.ml)

rvol <- respo_vol %>% 
  left_join(species) %>% 
  drop_na(postVol) %>% 
  mutate(pre_post = "post",
         Volume.ml = postVol - preVol) %>% 
  select(SpeciesID, Sp, SpRep, AT, ET, pre_post, Volume.ml) %>% 
  filter(SpeciesID != "LC2HH" & SpeciesID != "PA3HL") %>% # currently showing duplicate volumes - need to manage
  filter(Sp != "DN") %>%  # more accurate volume taken in ww_weight
  filter(Sp != "DB") # not used for growth

# volume %>% rbind(rvol) %>% 
#   group_by(SpeciesID, pre_post) %>% 
#   count() %>% 
#   filter(n > 1)

volume <- volume %>% 
  rbind(rvol) %>% 
  distinct() %>% 
  pivot_wider(names_from = pre_post, values_from = Volume.ml) %>%
  drop_na(post) %>%
  mutate(delVolume.ml = post - pre,
         pVolume = delVolume.ml/pre*100) %>%  # calculate weight difference over soak period
  select(-c("pre","post"))
  

# write csv with calculated values
species_cal <- rbind(bw, ww) %>%  # only bw until I have other species
  full_join(length) %>% 
  full_join(volume) %>%
  drop_na(delWeight.g) %>% 
  left_join(surfacearea)
  
# write_csv(species_cal, here("Data", "RespoFiles", "SpeciesMetadata_calculated_perday.csv"))

#############################
##### Assemblage Metrics #####
#############################

species <- species %>% 
  select(SpeciesID:SpRep)

AssemVolume <- respo_vol %>% 
  left_join(species) %>% 
  mutate(Volume.ml = postVol - preVol) %>% 
  select(SampleID, SpeciesID, Sp, SpRep, AT, ET, Volume.ml) %>% 
  # add in temporary volume for missing value for now
  mutate(Volume.ml = if_else(SpeciesID == "PR3HL", 25, Volume.ml)) %>% 
  group_by(SampleID) %>% 
  summarise(Volume.ml = sum(Volume.ml)) %>%
  distinct(SampleID, Volume.ml) # unique values for assemblages only

# PR9HH needed two tins

cor_afdw <- afdw %>% 
  mutate(Dry.g = DrywTin.g - Tin.g,
         AFDW.g = Dry.g - (AFDWwTin.g-Tin.g))

pr9hhafdw <- cor_afdw %>% 
  filter(SpeciesID == "PR9HH_2" | 
           SpeciesID == "PR9HH") %>% 
  mutate(Dry.g = sum(Dry.g),
         AFDW.g = sum(AFDW.g)) %>%  # combine dry and ash free dry weights of PR9HH data - adding for the same coral
  filter(SpeciesID == "PR9HH") # remove now redundant second weight

AssemAFDW <- cor_afdw %>% 
  filter(SpeciesID != "PR9HH" &
           SpeciesID != "PR9HH_2") %>% # remove old values of pr9hh
  rbind(pr9hhafdw) %>% # bind summed values of pr9hh
  left_join(metaSpecies) %>%
  # add in temporary afdw for missing value for now
  mutate(AFDW.g = if_else(SpeciesID == "LC10HH", 0.7, AFDW.g)) %>% 
  drop_na(SampleID) %>% # remove any orgs not used in respo
  group_by(SampleID, AT, ET) %>% 
  summarise(AFDW.g = sum(AFDW.g))
  


Respo_biometrics <- AssemVolume %>% 
  full_join(AssemAFDW) %>% 
  relocate(Volume.ml, .after = ET)

#############################
### Export csv
#############################
#write_csv(Respo_biometrics, here("Data","RespoFiles","AssemblageMetadata_calc.csv"))



