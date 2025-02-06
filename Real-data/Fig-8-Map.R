library(tidyverse)
library(sf)
library(ggthemes)
library(maps)
world <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))


cost_of_living <- read.csv("/Users/gthb3/Dropbox/PhD students/Shahriari/My code/cost-of-living.csv" ,sep = ",")

Cont = unique(cost_of_living$country)
Cont[Cont %in% c("United Kingdom", "United States", "Trinidad And Tobago", 
              "Saint Vincent And The Grenadines", 
              "Antigua And Barbuda", "Vatican City", "British Virgin Islands",
              "Saint Kitts And Nevis")] = 
  c("UK", "USA", "Tobago", "Saint Vincent","Antigua", "Vatican",
    "Virgin Islands, British", "Saint Kitts")


Cont = c(Cont, "Trinidad", "Democratic Republic of the Congo", 
         "Czech Republic", "Central African Republic", "Republic of Congo")

#"Gibraltar" "Hong Kong" "Tuvalu" are removed

world1 = subset(world, toupper(world$ID) %in% toupper(Cont))
#world$ID[toupper(world$ID) %in% toupper(Cont) == F]
#Cont[toupper(Cont) %in% toupper(world$ID) == F]

Group=scan()
3 2 4 4 2 1 2 4 2 2 4 4 3 1 2 3 1 3 3 3 2 4 2 2 2 3 2 2 4 
4 4 2 5 2 3 2 3 3 5 5 4 2 4 1 4 2 5 1 4 4 1 3 2 3 3 2 2 3 
3 2 1 4 1 1 2 2 1 1 2 3 3 2 4 2 3 2 4 1 2 4 2 2 1 2 1 3 2 
2 3 3 2 3 1 2 2 3 3 2 2 4 2 4 2 2 3 2 4 3 2 2 2 1 2 2 3 2 
2 3 2 4 4 2 2 2 2 2 1 4 2 3 3 2 4 3 2 2 1 1 2 3 4 3 2 1 2 
2 4 2 3 4 2 2 4 3 2 3 3 4 2 3 1 4 3 2 2 3 2 1 2 2 2 1 2 4 
4 2 3 4 4 2 2 3 1 4 2 2 2 2 2 2 2 2 3 2 3 2 4 2 5 2 4 4 4 
3 3 1

# 244 = 81
# 245 = 238


Col = ifelse(Group == 1, "firebrick1", 
             ifelse(Group == 2, "lightsalmon2", 
                    ifelse(Group == 3, "darkgoldenrod2", 
                           ifelse(Group == 4, "darkgoldenrod4", "brown"))))

CNMFA = ifelse(Group == 1, "Group 1", 
                 ifelse(Group == 2, "Group 2", 
                        ifelse(Group == 3, "Group 3", 
                               ifelse(Group == 4, "Group 4", "NA"))))

library(countrycode)
world1$ID2 = countryname(world1$ID, destination = 'iso3c')

sf_use_s2(FALSE)
states <- cbind(world1 , st_coordinates(st_centroid(world1)))


P1 = ggplot(world1) +
  geom_sf(color = "black", aes(fill = CNMFA)) +
  theme(legend.position = c(0.1, 0.4),  
        legend.box.background = element_rect(color = "gray92", linewidth = .3),
        legend.background = element_rect(fill = "gray92", 
                                         linewidth = 0.75, linetype = "solid")) + 
  geom_text(data = states, aes(X, Y, label = ID2), size = 1.85, angle = 0, 
            fontface = 2) +
  coord_sf(expand = T) +
  theme(axis.title=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank())

ggsave(filename = "Map.eps", 
       plot = P1, 
       path = "/Users/gthb3/Dropbox/PhD students/Shahriari/My code",
       device = "eps", 
       dpi = 300, 
       width = 20,
       height = 10, 
       units = "cm")

postscript('/Users/gthb3/Dropbox/PhD students/Shahriari/My code/Map.eps', 
           width = 27, height = 27)
P1 + theme(panel.spacing = unit(c(0, 0, 0, 0), "cm"),
           plot.margin = unit(c(0, -0.04, 0, -0.03), "null"))
dev.off()







### way 2

#ggplot() + geom_sf(data = world1, fill = Col, color = "black") + 
#  geom_text(data = states, aes(X, Y, label = ID2), size = 2, angle = 0, 
#            fontface = 2, size = 2, show.legend = "polygon") +
#  coord_sf(expand = T) +
#  theme(axis.title=element_blank(), 
#        axis.text=element_blank(), 
#        axis.ticks=element_blank()) +
#  ggtitle(" clustering")

