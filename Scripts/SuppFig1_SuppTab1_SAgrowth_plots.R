#### Supplemental Figure 1 and Supplemental Table 1: surface area-normalized growth
#### Created by Danielle Barnas


##########################################################
### Supplemental Figure 1
##########################################################

#############################
### LOAD LIBRARIES
#############################
library(tidyverse)
library(here)
library(patchwork)
library(stats)
library(lme4)
library(lmerTest) #Need this to get anova results from lmer
library(agricolae) # HSD.test()
library(car)
library(emmeans)
library(kableExtra)



#############################
### READ IN DATA
#############################

species <- read_csv(here("Data","RespoFiles","SpeciesMetadata_calculated_perday.csv"))
meta <- read_csv(here("Data","RespoFiles","SpeciesMetadata.csv")) %>% select(SpeciesID, FullSp)




#############################
### PROCESS SPECIES NAMES
#############################
species <- species %>% 
  select(-delTopLength.cm, -pLength, -delVolume.ml, -pVolume)

meta <- species %>% 
  left_join(meta) %>% 
  distinct() %>%
  mutate(Sp = if_else(Sp == "ME", "MG", 
                      if_else(Sp == "GS", "Por1", Sp))) %>% 
  mutate(FullSp = if_else(Sp == "MG", "Montipora grisea", 
                          if_else(Sp == "Por1", "Porifera unknown", FullSp)))

species <- species %>%
  mutate(Sp = if_else(Sp == "ME", "MG", 
                      if_else(Sp == "GS", "Por1", Sp))) %>% 
  left_join(meta)

species <- species %>% 
  mutate(PartSp = if_else(Sp == "DN", "D. nummiforme",
                          if_else(Sp == "Por1", "Porifera unknown", 
                                  if_else(Sp == "HO", "H. opuntia", 
                                          if_else(Sp == "LK", "L. kotschyanum",
                                                  if_else(Sp == "MG", "M. grisea",
                                                          if_else(Sp == "PA", "P. acuta",
                                                                  if_else(Sp == "PR", "P. rus",
                                                                          if_else(Sp == "VF", "V. fastigiata", ""))))))))) %>% 
  mutate(ET = factor(ET, levels = c("LOW", "HIGH")))

species <- species %>% 
  unite(Sp, AT, col = "Sp_AT", remove = F, sep = "_") %>% 
  unite(Sp, SpRep, AT, col = "SpRep_AT", remove = F, sep = "_") %>% 
  group_by(Sp, PartSp, ET) %>% 
  # normalize to days between measurements
  mutate(delWeight.g_day = delWeight.g/DaysInSitu)

# taxon-specific normalization metric
sa.sp <- c("PR", "PA", "MG", "LK")

species <- species %>% 
  filter(Sp %in% sa.sp) %>% 
  relocate(delWeight.g_day, .after = delWeight.g) %>% 
  mutate(growth.mg.cm2.d = delWeight.g_day/SA_cm2*1000) %>% 
  select(SpeciesID:ET, FullSp:growth.mg.cm2.d)
      
meanspecies <- species %>% 
  summarise(meanVal = mean(growth.mg.cm2.d, na.rm = TRUE), # already grouped by Sp and ET
            sd = sd(growth.mg.cm2.d, na.rm = TRUE),
            se = plotrix::std.error(growth.mg.cm2.d, na.rm = TRUE)) %>% 
  mutate(ET = factor(ET, levels = c("LOW", "HIGH"))) %>% 
  drop_na(meanVal)



#############################
### VISUALIZE
#############################
## create above plot in patchwork for labeling
weightPlotFun <- function(myfilter){
  data <- meanspecies %>% 
    filter(Sp == {{myfilter}})
  data2 <- species %>%
    filter(Sp == {{myfilter}})
  weightPlot <- data %>% 
    ggplot(aes(x = ET, y = meanVal, fill = ET)) +
    geom_point(data = data,
               shape = 21,
               size = 2,
               color = "black",
               aes(fill = ET),
               show.legend = FALSE) +
    geom_errorbar(data = data,
                  aes(ymin = meanVal-se, ymax = meanVal+se), width = 0.2) +
    geom_point(data = data2,
               position = position_jitterdodge(),
               shape = 21,
               size = 2,
               color = "black",
               aes(x = ET, y = growth.mg.cm2.d,
                   fill = ET),
               show.legend = FALSE,
               alpha = 0.1) +
    facet_wrap(~PartSp, scales = "free_y") +
    theme_bw() +
    scale_fill_manual(values = c("#0072B2", "#D55E00")) +
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 10,face = "italic"),
          panel.grid = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_blank(),
          legend.position = "none") +
    labs(#x = "SGD Exposure Treatment",
      color = "Assemblage \nTreatment")
  return(weightPlot)
}
a<-weightPlotFun("PA") + theme(legend.position = "right")
b<-weightPlotFun("PR")
c<-weightPlotFun("MG")
d<-weightPlotFun("LK")

layout <- '
AB
CD
'

weightPatch <- a + b + c + d +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") +
  plot_layout(guides = 'collect',
              design = layout) +
  theme(plot.tag = element_text(face = 'italic'))


weightPatch.2 <- wrap_elements(weightPatch) +
  labs(tag = expression("Growth (mg cm"^"-2"*"day"^"-1"*")")) +
  theme(
    plot.tag = element_text(size = rel(1.3), angle = 90),
    plot.tag.position = "left"
  )

wrap_elements(weightPatch.2) +
  labs(tag = "SGD Exposure Treatment") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom"
  )
weightPatch.2

# ggsave(here("Output","PaperFigures","SuppFig1.png"), weightPatch.2, width = 8, height = 5)




##########################################################
### Supplemental Table 1: Type III Analysis of Variance table displaying percent 
### change in weight of species pairs in either high or low SGD exposure (ET).
##########################################################


#############################
### REPEATED MEASURES ANOVA / PAIRING AS RANDOM EFFECT
#############################
# repeated measures anova
# models testing ET significance (not taking into account pairing...)
modelData <- species

# renumber pairs for 1:26 as necessary
modelData <- species %>%
  ungroup() %>% 
  select(SpeciesID:ET) %>% 
  count(AT, Sp,SpRep) %>%
  unite(Sp, SpRep, col = "UniqueSp", sep = "_", remove = F) %>% 
  group_by(Sp) %>% 
  mutate(newRep = seq_along(UniqueSp)) %>% 
  arrange(UniqueSp) %>% 
  select(Sp, SpRep, AT, UniqueSp, newRep) %>% 
  right_join(modelData)


Growth_Anova_Table <- tibble(Species = as.character(),
                             Source = as.character())

mySource <- as_tibble(c("~ ET + (1|Pairs)")) %>% rename(Source = value)

modTib <- tibble()

for(i in 1:4){
  mySp <- unique(modelData$Sp)[i] # get species name
  
  moddata <- modelData %>% filter(Sp == mySp) # get individual species data only
  
  # # whether AT nested or not
  # spcount <- moddata %>% count(Sp)
  # spcount <- as.numeric(spcount$n[1])
  # 
  # if(spcount > 26){
  #   mymodel <- lmer(data = moddata, pWeight ~ ET + (1|AT:SpRep)) # species in both assemblage types
  # } else {
  mymodel <- lmer(data = moddata, growth.mg.cm2.d ~ ET + (1|newRep))
  # }
  
  mod <- Growth_Anova_Table %>% 
    full_join(as_tibble(cbind(mySource,anova(mymodel)[1]))) %>% 
    full_join(as_tibble(cbind(mySource,anova(mymodel)[2]))) %>% 
    full_join(as_tibble(cbind(mySource,anova(mymodel)[3]))) %>% 
    full_join(as_tibble(cbind(mySource,anova(mymodel)[4]))) %>% 
    full_join(as_tibble(cbind(mySource,anova(mymodel)[5]))) %>% 
    full_join(as_tibble(cbind(mySource,anova(mymodel)[6]))) %>% 
    mutate(Species = mySp)
  
  colnames(modTib) <- colnames(mod) 
  modTib <- modTib %>% 
    rbind(mod)
}


sp.names <- meta %>% select(Sp, FullSp) %>% distinct()

fullmod <- modTib %>% 
  # use full species names
  rename(Sp = Species) %>% 
  left_join(sp.names) %>% 
  mutate(Sp = factor(Sp, levels = c("PA", "PR", "MG", "LK"))) %>% 
  arrange(Sp) %>%
  relocate(FullSp, .before = Source) %>% 
  
  # prep for table
  rename(`Growth ~ ET + (1|Pairs)`=FullSp,
         SS = 'Sum Sq',
         MS = 'Mean Sq',
         Fvalue = 'F value',
         p = 'Pr(>F)') %>% 
  select(-Source) %>% 
  mutate_at(.vars = vars(SS:DenDF), 
            .funs = ~signif(., 2)) %>% # all to 2 sig figs or round to 2
  
  mutate_at(.vars = vars(Fvalue), 
            .funs = ~signif(., 4)) %>% # all to 4 sig figs or round to 2
  
  mutate_at(.vars = vars(p), 
            .funs = ~signif(., 3)) %>% # all to 4 sig figs or round to 2
  
  mutate_at(.vars = vars(SS:MS), 
            .funs = ~if_else(. > 1000,format(.,scientific=T),format(., scientific = F))) %>% 
  
  mutate(Fvalue = if_else(is.na(Fvalue), "-", as.character(Fvalue))) %>% 
  mutate(p = if_else(is.na(p), "-", as.character(p))) %>% 
  rename('F' = Fvalue) %>% 
  mutate(star = if_else(Sp == "PR", "**", "")) %>% 
  #                       if_else(Sp == "VF", "***",
  #                               if_else(Sp == "HO", "*", "")))) %>% 
  unite(p, star, sep = " ", col = "p") %>%
  select(-Sp)


#############################
### CREATE KABLE TABLE
#############################

GrowthTable <- fullmod %>%
  kbl(align = c("l","r","r","r","r", "c", "l")) %>%
  #kbl() %>% 
  kable_classic(html_font = "Times New Roman",
                font_size = 20) %>% 
  row_spec(0, italic = TRUE, bold = TRUE) %>%  # header row
  column_spec(1, italic = TRUE)

GrowthTable # manual save as Supp Table 1
