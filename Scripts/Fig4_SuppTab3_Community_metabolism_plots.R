#### Figure 4 & Supplemental Table 3: Mean rates of ecosystem metabolism parameters (GP, NP, R, and NC) 
####            with standard error and raw associated data for each assemblage type 
####            in two SGD exposure treatments.
#### Created by Danielle Barnas


##########################################################
### Figure 4
##########################################################

#############################
### LOAD LIBRARIES
#############################

library(tidyverse)
library(here)
library(patchwork)
library(PNWColors)
library(agricolae) # HSD.test()
library(kableExtra)
library(plotrix)
library(emmeans)



#############################
### READ IN DATA
#############################

myNEC_NEP <- read_csv(here("Data", "RespoFiles", 'All_EcoMet_Rates.csv'))

myNEC_NEP <- myNEC_NEP %>%
  mutate(ET = factor(ET, levels = c("LOW", "HIGH")),
         AT = factor(AT, levels = c("LOW", "HIGH")))

meanRates <- myNEC_NEP %>% 
  group_by(AT,ET,P_R) %>% 
  summarise(meanVal = mean(Values),
            sd = sd(Values),
            se = std.error(Values))

#############################
### VISUALIZATION
#############################

# palette
my_pal <- pnw_palette(name="Starfish",n=2,type="discrete")

# plot function


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
                alpha = 0.2,
                position = position_jitterdodge(dodge.width = 1)) +
    geom_errorbar(data = data,
                  position = position_dodge(width = 1),
                  aes(ymin = meanVal-se, ymax = meanVal+se), width = 0.2) +
    geom_point(data = data,
               shape = 21,
               size = 3.1,
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


# apply function
gpp<-RatesPlotFun(myValue = "GP", renameVal = "Gross photosynthesis") + labs(y = expression("Rate (mmol O"[2]*" g"^-1*" hr"^-1*")"))
npp<-RatesPlotFun(myValue = "NP", renameVal = "Net photosynthesis")  + labs(y = expression("Rate (mmol O"[2]*" g"^-1*" hr"^-1*")"))
rp<-RatesPlotFun(myValue = "R", renameVal = "Respiration")  + labs(y = expression("Rate (mmol O"[2]*" g"^-1*" hr"^-1*")"))
ncp<-RatesPlotFun(myValue = "NC", renameVal = "Net calcification") + labs(y = expression("Rate ("*mu*"mol CaCO"[3]*" g"^-1*" hr"^-1*")"))

# patch plots
EcoFunPlot <- (gpp + npp) / (rp + ncp) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") +
  plot_layout(guides = 'collect')  


#############################
### STATISTICAL TEST
#############################


EM_Anova_Table <- tibble(Parameter = as.character(),
                         Source = as.character())

mySource <- as_tibble(c("AT", "ET", "AT x ET Interaction", "Residuals")) %>% rename(Source = value)

modTib <- tibble()
myContrasts <- tibble()

for(i in 1:4){
  myParam <- unique(myNEC_NEP$P_R)[i] # get parameter name
  
  moddata <- myNEC_NEP %>% filter(P_R == myParam) # get individual parameter data only
  
  mymodel <- lm(data = myNEC_NEP %>% filter(P_R == myParam), Values ~ AT*ET)
  
  # HSD.test(mymodel, "AT", console = TRUE)
  tempContrasts <- as_tibble(emmeans(mymodel, pairwise ~ AT*ET, adjust="Tukey")$contrasts) %>% 
    mutate(P_R = myParam)
  myContrasts <- myContrasts %>% 
    rbind(tempContrasts)
  
  mod <- EM_Anova_Table %>% 
    full_join(as_tibble(cbind(mySource,anova(mymodel)[1]))) %>% # df
    full_join(as_tibble(cbind(mySource,anova(mymodel)[2]))) %>% # SS
    full_join(as_tibble(cbind(mySource,anova(mymodel)[3]))) %>% # MS
    full_join(as_tibble(cbind(mySource,anova(mymodel)[4]))) %>% # F
    full_join(as_tibble(cbind(mySource,anova(mymodel)[5]))) %>% # p
    mutate(Parameter = myParam)
  
  #summary(mod)
  #HSD.test(mod, "AT", console=TRUE)
  
  colnames(modTib) <- colnames(mod) 
  modTib <- modTib %>% 
    rbind(mod)
}



fullmod <- modTib %>% 
  rename(df = 'Df',
         SS = 'Sum Sq',
         MS = 'Mean Sq',
         Fvalue = 'F value',
         p = 'Pr(>F)') %>% 
  mutate_at(.vars = vars(SS:p), .funs = ~signif(., 2)) %>% 
  mutate(Fvalue = if_else(is.na(Fvalue), "-", as.character(Fvalue))) %>% 
  mutate(star = case_when(p < 0.001 ~ '***',
                          p < 0.01 ~ '**',
                          p < 0.05 ~ '*',
                          p >= 0.05 ~ '')) %>% 
  mutate(p = if_else(is.na(p), " - ", as.character(p))) %>% 
  unite(p, star, sep = " ", col = "p", na.rm = TRUE) %>% 
  rename('F' = Fvalue) %>% 
  mutate(Parameter = factor(Parameter, levels = c("GP", "NP", "R", "NC"))) %>% 
  arrange(Parameter)


# plot with contrasts


# view contrasts from tukey test
myContrasts %>% 
  select(P_R, contrast, p.value)


# apply function
gpp<-RatesPlotFun(myValue = "GP", renameVal = "Gross photosynthesis") + 
  labs(y = expression("Rate (mmol O"[2]*" g"^-1*" hr"^-1*")")) +
  geom_text(aes(label=c("a", "a", "b", "b")), 
            position = position_dodge(width=1), 
            vjust=c(-2.1, -2.0,-1.5, -1.5))

npp<-RatesPlotFun(myValue = "NP", renameVal = "Net photosynthesis")  + 
  labs(y = expression("Rate (mmol O"[2]*" g"^-1*" hr"^-1*")")) +
  geom_text(aes(label=c("a", "a", "b", "b")), 
            position = position_dodge(width=1), 
            vjust=c(-2.8, -2.1,-1.9, -1.6))

rp<-RatesPlotFun(myValue = "R", renameVal = "Respiration")  + 
  labs(y = expression("Rate (mmol O"[2]*" g"^-1*" hr"^-1*")")) +
  geom_text(aes(label=c("a", "ab", "a", "ac")), 
            position = position_dodge(width=1), 
            vjust=c(-2.5, -3.0,-2.1, -1.8))

ncp<-RatesPlotFun(myValue = "NC", renameVal = "Net calcification") + 
  labs(y = expression("Rate ("*mu*"mol CaCO"[3]*" g"^-1*" hr"^-1*")")) +
  geom_text(aes(label=c("a", "a", "b", "b")), 
            position = position_dodge(width=1), 
            vjust=c(-2.1, -2.1,-1.5, -1.5))

# patch plots
EcoFunPlot <- (gpp + npp) / (rp + ncp) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") +
  plot_layout(guides = 'collect') 

EcoFunPlot


# ggsave(here("Output","PaperFigures","Fig4_EcoFunction_3.11.png"), EcoFunPlot, device = "png", width = 8, height = 5)



##########################################################
### Supplemental Table 3: Two-way Analysis of Variance table expressing independent variable and 
### interaction term effects on rates of ecosystem metabolism parameters GP, NP, Rd, and NC.
##########################################################



#############################
### CREATE KABLE TABLE
#############################

EcoFunTable <- fullmod %>%
  select(-Parameter) %>% 
  kbl(align = c("l","l","r","l","l")) %>%
  kable_classic(html_font = "Times New Roman",
                font_size = 20) %>% 
  row_spec(0, italic = TRUE, bold = TRUE) %>%  # header row
  pack_rows("GP", 1, 4) %>% 
  pack_rows("NP", 5, 8) %>% 
  pack_rows("R", 9, 12) %>% 
  pack_rows("NC", 13, 16)

EcoFunTable


# EcoFunTable %>% 
  # as_image(file = here("Output", "Thesis_Figures_Output", "EcoFunTable.png"))

