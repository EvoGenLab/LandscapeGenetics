#
# Yellow-vented Bulbul Ne radius analysis
# 
# Jun 27, 2024
#
# Will Brooks
####

#setwd("~/Desktop/RADseq_YVBU/NeEstimator")

# Clear environment and load packages ----

rm(list = ls()) # clean environment

# Packages

library(tidyverse)
library(terra)
library(sf)
library(landscapemetrics)
library(ggplot2)
library(GGally)

# Load and explore data ----

# Ne
ne <- read.csv("../NeEstimator/NEsummary.csv") %>%
  arrange(Pop)

sites <- read.csv("../Genetic data/yvb_study_sites.csv") %>%
  subset(abbrev %in% ne$Pop) %>%
  arrange(abbrev)

# load raster data
veg5 <- rast("../Landscape data/veg5_merge78.tif") %>% # Landcover from Hansen et al.
  as.factor()
classes <- data.frame(class = c(0,3,4,5,7),
                      cover = c("water","urban",
                                "urban_green","scrub",
                                "secondary"))
levels(veg5) <- classes

# Create spatial object for sampling sites
sites.sf <- st_as_sf(sites, coords = c("Long", "Lat"), crs = 4326)
sites.utm <- st_transform(sites.sf, crs(veg5))

# Make plot
plot(veg5, 
     col=c("white","azure3","tan3","darkolivegreen1","darkgreen"))

points(sites.utm, cex =ne[,2]/250, pch = 21, bg = "red")

# add legend
legend_sizes <- c(min(ne[,2]/250),max(ne[,2]/250))
add_legend("bottom", 
       legend = round(legend_sizes*250), 
       pch = 21, 
       pt.cex = legend_sizes, 
       pt.bg = "red", 
       title = "Ne")

# # # add legend to other plot
# legend_sizes <- c(min(ne[,2]/250),max(ne[,2]/250))
# legend("bottom",
#            legend = round(legend_sizes*250),
#            pch = 21,
#        cex = 0.7,
#            pt.cex = legend_sizes,
#            pt.bg = "red",
#            title = "Ne")

# Radius land cover predict ne ----

buffer <- sample_lsm(landscape = veg5, 
                     what = "lsm_c_pland",
                     shape = "circle",
                     y = sites.utm, 
                     size = 250) %>% 
  left_join(classes,by = "class") %>% 
  pivot_wider(id_cols = plot_id, 
              names_from = cover, 
              values_from = value)
buffer[is.na(buffer)] <- 0

ne_buffer <- cbind(ne,buffer[,2:6])

# check distributions and correlation of predictors
ne_buffer %>% 
  dplyr::select(-Pop,-Ne_crit0.1) %>%
  ggpairs()

# I'll remove 0 as it has a ridiculous distribution

# Check relationships between predictors and abundance
plot(Ne_crit0.1 ~scrub, dat = ne_buffer)
plot(Ne_crit0.1 ~secondary, dat = ne_buffer)
plot(Ne_crit0.1 ~urban, dat = ne_buffer)
plot(Ne_crit0.1 ~urban_green, dat = ne_buffer)

# Model Ne against buffer covs
lm.scrub <- lm(Ne_crit0.1 ~(scrub), dat = ne_buffer)
lm.young_secondary <- lm(Ne_crit0.1 ~(secondary), dat = ne_buffer)
lm.urban <- lm(Ne_crit0.1 ~(urban), dat = ne_buffer)
lm.urban_green <- lm(Ne_crit0.1 ~(urban_green), dat = ne_buffer)

summary(lm.scrub)
summary(lm.young_secondary)
summary(lm.urban)
summary(lm.urban_green)

# I think this is nonsense

# Area predict Ne ----

# get only secondary forest
secondary <- veg5 >= 7
plot(secondary)

area <- extract_lsm(landscape = secondary, 
                     what = "lsm_p_area",
                     y = sites.utm) 


ne_area <- cbind(ne,area)

plot(Ne_crit0.1 ~ log(value), data = ne_area)

lm_area <- lm(Ne_crit0.1 ~ log(value), data = ne_area)
summary(lm_area)

# Latitude predict Ne ----

ne_lat <- cbind(ne_area,sites)

plot(Ne_crit0.1 ~ Lat, data = ne_lat)

lm_lat <- lm(Ne_crit0.1 ~ Lat, data = ne_lat)
summary(lm_lat)

(ne.lat.plot <- ggplot(data = ne_lat, 
                          mapping = aes(x = Lat, 
                                        y = Ne_crit0.1)) +
    geom_point() +
    theme_bw( ) +
    theme(plot.margin = unit(c(1,1,1,1), "cm")) +
    labs(x = "° Latitude", 
         y = bquote(N[e])) )

ggsave("ne_lat.jpeg",
       device = "jpeg",
       width = 4.5,                     # Set width (in inches)
       height = 4,                     # Set height (in inches)
       dpi = 300                       # Set DPI for raster formats
)

# ms version - MANUSCRIPT
neplot <- ne.lat.plot + 
  geom_smooth(method="lm", se=FALSE, linetype=2, color="black", linewidth=0.75) +
  theme(panel.grid=element_blank())
neplot
saveRDS(neplot, file="fig2C.RDS")


