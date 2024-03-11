#### Table 1: High and low SGD assemblage treatments, with corresponding species and  
####          functional traits. Novel species from each treatment are indicated in bold.
#### Created by Danielle Barnas


##########################################################
### Table 1
##########################################################

#############################
### LOAD LIBRARIES
#############################
library(tidyverse)
library(here)
library(kableExtra)
#devtools::install_github("renkun-ken/formattable")
library(formattable)




#############################
### READ IN DATA
#############################
spfedf <- read_csv(here("Data", "FunDiv", "Species_FE.csv"))
comp <- read_csv(here("Data", "FunDiv", "Species_Abundances_wide.csv"))
rawcomp <- read_csv(here("Data", "FunDiv", "Species_Composition_2022.csv"))


#############################
### REDUCE DATA TO ASSEMBLAGE SPECIES
#############################
spfe <- spfedf %>% 
  mutate(Taxa = if_else(Taxa == "Grey Sponge", "Porifera 1", Taxa)) %>% # rename to porifera for publication fig
  rename(Species = Taxa) %>% 
  separate(FE, into = c("Phyla", 
                        "Morphology", 
                        "Calcification Type", 
                        "Energetic Resource"),
           sep = ",")


# distinct exposure areas
highsite <- c("V14", "V15", "V11", "V18")
lowsite <- c("V2", "V3", "V4", "V16")

# identified taxa that wouldn't be used for assemblages
unused_sp <- tibble(Taxa = c("Turf", "Crustose Corallines", "Heteractis magnifica", "Padina boryana")) #PB died in all treatments

abund <- rawcomp %>% 
  filter(Location == "Varari") %>% 
  filter(Taxa %in% spfedf$Taxa) %>% 
  filter(CowTagID %in% highsite | CowTagID %in% lowsite) %>% 
  select(CowTagID, Taxa, SpeciesCounts) %>% 
  mutate(treatment = if_else(CowTagID %in% highsite, "High", "Low")) %>% 
  group_by(treatment,Taxa) %>% 
  summarise(counts = sum(SpeciesCounts)) %>% 
  ungroup() %>% 
  group_by(treatment) %>% 
  mutate(totalCounts = sum(counts),
         cover = counts / totalCounts * 100) %>% 
  select(treatment, Taxa, cover) %>% 
  anti_join(unused_sp) %>% 
  arrange(treatment, desc(cover))

# select top 9 species from each exposure area
highAb <- abund %>% 
  mutate(Taxa = if_else(Taxa == "Grey Sponge", "Porifera 1", Taxa)) %>% 
  filter(treatment == "High") %>% 
  mutate(cover = round(cover,1)) %>% 
  filter(Taxa != "Halimeda opuntia") # ~same pcover as LC, chose LC for genera diversity
highAb <- highAb[1:8,]

lowAb <- abund %>% 
  mutate(Taxa = if_else(Taxa == "Grey Sponge", "Porifera 1", Taxa)) %>% 
  filter(treatment == "Low") %>% 
  mutate(cover = round(cover,1)) 
lowAb <- lowAb[1:8,]

species_rank <- rbind(lowAb, highAb) %>% 
  rename(`Benthic % Cover` = cover,
         Species = Taxa,
         `SGD Exposure` = treatment)



#############################
### EXTEND FUNCTIONAL ENTITY NAMES
#############################
spfe.long <- spfe %>% 
  mutate(Morphology = if_else(Morphology == "Br", "Branching",
                              if_else(Morphology == "Mas", "Massive",
                                      if_else(Morphology == "Enc", "Encrusting",
                                              if_else(Morphology == "Cushion", "Cushion-like",
                                                      if_else(Morphology == "Poly", "Polypoid", Morphology)))))) %>% 
  mutate(`Calcification Type` = if_else(`Calcification Type` == "Herm", "Heratypic",
                                        if_else(`Calcification Type` == "NC", "Non-calcifying",
                                                if_else(`Calcification Type` == "AC", "Articulated",
                                                        if_else(`Calcification Type` == "Non-AC", "Non-articulated",`Calcification Type`))))) %>% 
  mutate(`Energetic Resource` = if_else(`Energetic Resource` == "Mix", "Mixotrophic", 
                                        if_else(`Energetic Resource` == "Auto", "Autotrophic",
                                                if_else(`Energetic Resource` == "Het", "Heterotrophic", `Energetic Resource`))))



#############################
### JOIN SP AND FE -> TABLE
#############################

species_rank_table <- species_rank %>% 
  left_join(spfe.long) %>% 
  relocate(`SGD Exposure`, .before = Species) %>% 
  relocate(`Benthic % Cover`, .after = `Energetic Resource`) %>% 
  arrange(desc(`SGD Exposure`), desc(`Benthic % Cover`)) %>% 
  mutate(Species = if_else(Species == "Porifera 1", "Porifera unknown", Species)) # rename Poifera 1 as Porifera unknown


formattable(species_rank_table, list(
  `SGD Exposure` = color_tile("transparent", "aliceblue"),
  Species = color_bar("white"),
  Phyla = color_bar("white"),
  Morphology = color_bar("white"),
  `Calcification Type` = color_bar("white"), 
  `Energetic Resource` = color_bar("white"),
  `Benthic % Cover` = color_bar("lightblue")
))

# remove first column (group later)
sp_rank_table <- species_rank_table[2:7]

SpeciesTable <- sp_rank_table %>% 
  kbl() %>% 
  kable_classic(html_font = "Times New Roman",
                font_size = 15) %>% 
  row_spec(0, italic = TRUE, bold = TRUE) %>%  # header row
  row_spec(5, bold = TRUE) %>% # novel Low species
  row_spec(6, bold = TRUE) %>% 
  row_spec(8, bold = TRUE) %>% 
  row_spec(12, bold = TRUE) %>% # novel High species
  row_spec(15, bold = TRUE) %>% 
  row_spec(16, bold = TRUE) %>% 
  column_spec(1, italic = TRUE) %>% 
  column_spec(6, bold = TRUE,color = case_when(sp_rank_table$`Benthic % Cover` > 8 ~ 'white',
                                               TRUE ~ 'black')) %>% 
  column_spec(6, background = case_when(#sp_rank_table$`Benthic % Cover` > 19 ~ '#015798', # '#145DA0', 
    sp_rank_table$`Benthic % Cover` > 19 ~ '#0277BD', 
    sp_rank_table$`Benthic % Cover` >= 10 ~ '#0288D1', 
    sp_rank_table$`Benthic % Cover` >= 9 ~ '#039BE5', 
    sp_rank_table$`Benthic % Cover` >= 6 ~ '#03A9F4', 
    sp_rank_table$`Benthic % Cover` >= 4 ~ '#29B6F6', 
    sp_rank_table$`Benthic % Cover` > 3.5 ~ '#4FC3F7', 
    sp_rank_table$`Benthic % Cover` >= 3 ~ '#81D4FA', 
    #sp_rank_table$`Benthic % Cover` >= 2 ~ '#B3E5FC', 
    sp_rank_table$`Benthic % Cover` >= 2 ~ '#99EDFF', 
    sp_rank_table$`Benthic % Cover` >1.5 ~ '#B2FCFF', 
    #sp_rank_table$`Benthic % Cover` >1.5 ~ '#B3E5FC', 
    sp_rank_table$`Benthic % Cover` <1.5 ~ '#CCFFFF', 
    TRUE ~ 'black')) %>% 
  pack_rows("Low SGD Assemblage", 1, 8) %>% 
  pack_rows("High SGD Assemblage", 9, 16)

SpeciesTable


# SpeciesTable %>%
#   as_image(file = here("Output", "PaperFigures", "SpeciesFunctionTable.png"))
