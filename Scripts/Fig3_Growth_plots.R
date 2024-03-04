#### Figure 3 & Supplemental Table 1: Mean percent change in weight with 
####            standard error and raw associated data points for species pairs 
####            in two SGD exposure treatments.
#### Created by Danielle Barnas


##########################################################
### Figure 3
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

species <- read_csv(here("Data", "Growth", "All_Weight_pChange.csv"))
meta <- read_csv(here("Data","RespoFiles","SpeciesMetadata.csv"))



#############################
### PROCESS SPECIES NAMES
#############################
meta <- meta %>% 
  select(FullSp, Sp) %>% 
  distinct() %>%
  mutate(Sp = if_else(Sp == "ME", "MG", 
                      if_else(Sp == "GS", "Por1", Sp))) %>% 
  mutate(FullSp = if_else(Sp == "MG", "Montipora grisea", 
                          if_else(Sp == "Por1", "Porifera 1", FullSp)))

species <- species %>%
  mutate(Sp = if_else(Sp == "ME", "MG", 
                      if_else(Sp == "GS", "Por1", Sp))) %>% 
  left_join(meta)

species <- species %>% 
  mutate(PartSp = if_else(Sp == "DN", "D. nummiforme",
                          if_else(Sp == "Por1", "Porifera 1", 
                                  if_else(Sp == "HO", "H. opuntia", 
                                          if_else(Sp == "LK", "L. kotschyanum",
                                                  if_else(Sp == "MG", "M. grisea",
                                                          if_else(Sp == "PA", "P. acuta",
                                                                  if_else(Sp == "PR", "P. rus",
                                                                          if_else(Sp == "VF", "V. fastigiata", ""))))))))) %>% 
  mutate(ET = factor(ET, levels = c("LOW", "HIGH")))

meanspecies <- species %>% 
  unite(Sp, AT, col = "Sp_AT", remove = F, sep = "_") %>% 
  unite(Sp, SpRep, AT, col = "SpRep_AT", remove = F, sep = "_") %>% 
  group_by(Sp,PartSp, ET) %>% 
  summarise(meanVal = mean(pWeight),
            sd = sd(pWeight),
            se = plotrix::std.error(pWeight)) %>% 
  mutate(ET = factor(ET, levels = c("LOW", "HIGH")))



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
               aes(x = ET, y = pWeight,
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
      #y = "% Weight change",
      color = "Assemblage \nTreatment")
  return(weightPlot)
}
a<-weightPlotFun("PA") + theme(legend.position = "right")
b<-weightPlotFun("PR")
c<-weightPlotFun("MG")
d<-weightPlotFun("VF") #+ theme(axis.title.y = element_text(size = 14)) + labs(y=expression("% "*Delta*" Weight"))
e<-weightPlotFun("HO")
f<-weightPlotFun("LK")
g<-weightPlotFun("Por1")
h<-weightPlotFun("DN") #+ theme(axis.title.x = element_text(size = 14)) + labs(x = "SGD Exposure Treatment")

layout <- '
ABCD
EFGH
'

weightPatch <- a + b + c + d + e + f + g + h +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") +
  plot_layout(guides = 'collect',
              design = layout) +
  theme(plot.tag = element_text(face = 'italic'))


weightPatch.2 <- wrap_elements(weightPatch) +
  labs(tag = expression("% "*Delta*" Weight")) +
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

# ggsave(here("Output","PaperFigures","Fig3_Weight_Change_long.png"), weightPatch.2, width = 8, height = 5)




##########################################################
### Supplemental Table 1: Type III Analysis of Variance table displaying changes in 
### growth of species pairs in either high or low SGD exposure (ET).
##########################################################


#############################
### REPEATED MEASURES ANOVA / PAIRING AS RANDOM EFFECT
#############################
# repeated measures anova
# models testing ET significance (not taking into account pairing...)
modelData <- species

# renumber pairs for 1:26 as necessary
modelData <- species %>% 
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

for(i in 1:8){
  mySp <- unique(modelData$Sp)[i] # get species name
  
  moddata <- modelData %>% filter(Sp == mySp) # get individual species data only
  
  # # whether AT nested or not
  # spcount <- moddata %>% count(Sp)
  # spcount <- as.numeric(spcount$n[1])
  # 
  # if(spcount > 26){
  #   mymodel <- lmer(data = moddata, pWeight ~ ET + (1|AT:SpRep)) # species in both assemblage types
  # } else {
  mymodel <- lmer(data = moddata, pWeight ~ ET + (1|newRep))
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



fullmod <- modTib %>% 
  # use full species names
  rename(Sp = Species) %>% 
  left_join(meta) %>% 
  mutate(Sp = factor(Sp, levels = c("PA", "PR", "MG", "VF", "HO", "LK", "Por1", "DN"))) %>% 
  arrange(Sp) %>%
  relocate(FullSp, .before = Source) %>% 
  
  # prep for table
  rename(`Growth ~ ET + (1|Pairs)`=FullSp,
         SS = 'Sum Sq',
         MS = 'Mean Sq',
         Fvalue = 'F value',
         p = 'Pr(>F)') %>% 
  select(-Source) %>% 
  mutate_at(.vars = vars(SS:p), 
            .funs = ~signif(., 2)) %>% # all to 2 sig figs or round to 2
  
  mutate_at(.vars = vars(SS:MS), 
            .funs = ~if_else(. > 1000,format(.,scientific=T),format(., scientific = F))) %>% 
  
  mutate(Fvalue = if_else(is.na(Fvalue), "-", as.character(Fvalue))) %>% 
  mutate(p = if_else(is.na(p), "-", as.character(p))) %>% 
  rename('F' = Fvalue) %>% 
  mutate(star = if_else(Sp == "PR", "**",
                        if_else(Sp == "VF", "***",
                                if_else(Sp == "HO", "*", "")))) %>% 
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

GrowthTable

# GrowthTable %>% 
#   as_image(file = here("Output", "Thesis_Figures_Output", "GrowthAnovaTable.png"))


