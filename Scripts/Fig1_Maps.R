#### Figure 1. Map of Moorea, French Polynesia
#### Created by Danielle Barnas


#############################
### LOAD LIBRARIES
#############################
library(here)
library(tidyverse)
library(ggmap)
library(maptools)
#devtools::install_github('oswaldosantos/ggsn')
library(ggsn) # removed from CRAN. Use install link above
library(patchwork)


##########################################################
### Figure 1A: Map of Moorea
##########################################################

#############################
### READ IN DATA
#############################
meta <- read_csv(here("Data", "Full_Metadata.csv")) 



##### MAP SITE LOCATIONS #####
meta <- meta %>% 
  # only keep platform locations
  filter(CowTagID == "VSEEP")



# mean lat and long for the maps
LocationGPS <- meta %>%
  group_by(Location) %>% 
  summarise(lon = median(lon, na.rm = TRUE),
            lat = median(lat, na.rm = TRUE))



# Moorea Base Map
MooreaMap <- get_map('Moorea', maptype = 'satellite', zoom = 12)


MooreaMapPlot <- ggmap(MooreaMap) + # base map
  
  geom_rect(data=LocationGPS, aes(xmin=lon[1] - 0.006, xmax=lon[1] + 0.006, ymin=lat[1] - 0.01, ymax=lat[1] + 0.01), color="white", alpha=0, size = 1) +
  
  labs(x = "Longitude", y = "Latitude") +
  
  geom_text(data = LocationGPS, aes(label = Location), color = "white", hjust = -0.4, size = 8) + # adds Location names to the right of the boxes
  
  geom_rect(aes(xmin = -149.842, xmax = -149.724,
                ymin = -17.638, ymax = -17.62), fill = "white") +
  
  ggsn::scalebar(x.min = -149.9, x.max = -149.74,
                 y.min = -17.625, y.max = -17.45, 
                 transform = TRUE, model = 'WGS84',
                 dist = 5, dist_unit = "km",
                 location = "bottomright",
                 st.dist = 0.04) +
    
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.title = element_text(size = 16),
        axis.text = element_text(size = 13))


MooreaMapPlot


##########################################################
### Figure 1B: Map of Platform Locations
##########################################################

##### READ IN DATA #####
meta <- read_csv(here("Data", "Full_Metadata.csv")) 



##### MAP SITE LOCATIONS #####
meta <- meta %>% 
  # only keep platform locations
  filter(CowTagID == "VSEEP" |
           CowTagID == "V17" |
           CowTagID == "V13") %>%  # low site roughly at v13
  mutate(lat = if_else(CowTagID == "V13", lat + 0.0001, lat),
         lon = if_else(CowTagID == "V13", lon + 0.0001, lon), # move point on map up to ~real in situ location
         CowTagID = if_else(CowTagID == "V17", "High\nSGD", CowTagID), # change site id names
         CowTagID = if_else(CowTagID == "V13", "Low\nSGD", CowTagID))



# isolate seep point for mapping
seeppt <- meta %>%
  filter(CowTagID == "VSEEP")

# mean lat and long for the maps
LocationGPS <- meta %>%
  summarise(lon = median(lon, na.rm = TRUE),
            lat = median(lat, na.rm = TRUE))


# Base Map
VarariBaseMap<-get_map(LocationGPS,
                       maptype = 'satellite',
                       zoom = 19)

# palette for high and low locations (match Fig1C colors)
mypal2 <- c("#0072B2", "#D55E00")

# Plot map
VmapSites <- ggmap(VarariBaseMap) +
  labs(x = "Longitude", y = "Latitude") +  #label x and y axes
  geom_label(data = meta %>% filter(CowTagID != "VSEEP"),
             aes(x = lon, y = lat,
                 label = CowTagID),
             color = "white",
             fill = mypal2,
             size = 4) +
  # add the seep point separately
  geom_point(data = seeppt,
             aes(x = lon, y = lat + 0.00003),
             fill = "white",
             color = "black",
             size = 16,
             shape = 21) +
  geom_text(data = seeppt,
            aes(x = lon, y = lat + 0.00003),
            label = "Seep",
            #fill = "white",
            size = 4) +
  
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        axis.title.y = element_blank()) +
  # add arrow outline
  geom_segment(aes(xend = -149.8999, yend = -17.5404, x = -149.8991, y = -17.5409),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 3, 
               color = "black") +
  # add flow direction arrow
  geom_segment(aes(xend = -149.8999, yend = -17.5404, x = -149.8991, y = -17.5409),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 2, 
               color = "white") +
  geom_rect(aes(xmin=-149.90005, xmax=-149.89952,
                ymin = -17.54109, ymax=-17.54095), fill = "white") +
  scalebar(x.min = -149.9000, x.max = -149.8985,
           y.min = -17.5410, y.max = -17.5395,
           location = "bottomleft", 
           dist = 20,
           dist_unit = "m",
           transform = TRUE,
           model = 'WGS84',
           st.dist = 0.03)

VmapSites 




#############################
### PATCHWORK PLOTS
#############################

Fig1plotA <- (MooreaMapPlot) +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '(',
                  tag_suffix = ')') +
  theme(plot.tag = element_text(face = 'italic'))

Fig1plotB <- (VmapSites) +
  plot_annotation(tag_levels = 'b',
                  tag_prefix = '(',
                  tag_suffix = ')') +
  theme(plot.tag = element_text(face = 'italic'))

Fig1plotAB <- Fig1plotA + Fig1plotB


ggsave(here("Output", "PaperFigures", "Figure_1ab.png"), Fig1plotAB, device = "png", height = 5, width = 9)

