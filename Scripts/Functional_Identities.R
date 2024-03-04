# FUNCTIONAL IDENTITIES OF ASSEMBLAGE TYPES
# Created by Danielle Barnas
# Created on June 4, 2023

#############################
### LOAD LIBRARIES
#############################

library(tidyverse)
library(here)
library(patchwork)
library(PNWColors)




#############################
### READ IN DATA
#############################

fe <- read_csv(here("Data", "FunDiv", "Distinct_FE.csv"))
ft <- read_csv(here("Data", "FunDiv", "Distinct_Taxa.csv"))
species.fe <- read_csv(here("Data", "FunDiv", "Species_FE.csv"))
meta <- read_csv(here("Data", "RespoFiles", "SpeciesMetadata.csv"))
metaCalc <- read_csv(here("Data", "RespoFiles", "SpeciesMetadata_calculated.csv"))
rVolume <- read_csv(here("Data", "RespoFiles", "RespoVolume.csv"))

#############################
### CLEAN DATA
#############################

# species used in assemblages
meta <- meta %>% 
  mutate(FullSp = if_else(FullSp == "Grey sponge", "Grey Sponge", FullSp)) %>% 
  mutate(FullSp = if_else(FullSp == "Montipora efflorescans", "Montipora efflorescens", FullSp)) %>% 
  filter(FullSp != "Halimeda incrassata")
species <- meta %>% 
  select(FullSp) %>% 
  distinct()
species <- species$FullSp
meta <- meta %>% 
  filter(FullSp %in% species) %>% 
  select(FullSp, Sp) %>% 
  distinct()
metaSp <- meta$Sp

# get species ID's used in respo
metaCalc <- metaCalc %>% 
  select(SpeciesID, Sp) %>% 
  filter(Sp %in% metaSp) %>% 
  left_join(meta)

# assign species ID with assemblage ID
as.species <- rVolume %>% 
  select(SampleID, AT, ET, SpeciesID)
# need to grab DB separately
DB <- rVolume[grep("DB", rVolume$SpeciesID), ]

as.species <- metaCalc %>% 
  left_join(as.species)
DB <- DB %>% 
  select(SampleID, AT, ET, SpeciesID) %>% 
  mutate(Sp = "DB", 
         FullSp = "Dictyota bartayresiana")
as.species <- as.species %>% 
  rbind(DB)
  
# join species and fe in assemblages
as.species.fe <- as.species %>% 
  rename(Taxa = FullSp) %>% 
  left_join(species.fe) %>% 
  left_join(fe) %>% 
  left_join(ft) %>% 
  drop_na(SampleID) # remove SpeciesID not used in respo

# only keep traits used in analysis
red.as.species.fe <- as.species.fe %>% 
  select(SpeciesID:ER)
View(red.as.species.fe)

# view total FE per AT
red.as.species.fe %>% 
  group_by(AT, ET) %>% 
  count(FE) %>%
  mutate(n = 1) %>%  # now get FE richness
  summarise(n = sum(n))

# view distinct FE per AT
View(red.as.species.fe %>% 
  distinct(AT, FE, Taxon_Group, Morph2, Calc, ER) %>% 
  arrange(AT,FE))

# view redundancy of FE in AT
View(red.as.species.fe %>% 
  group_by(AT, SampleID) %>% 
  count(FE))

# view species for each FE in AT
View(red.as.species.fe %>% 
       distinct(Taxa,FE,AT))

