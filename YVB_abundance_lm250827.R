#
# Yellow-vented Bulbul abundance analysis
# 
# Jun 20, 2024
#
# Will Brooks
####

setwd("~/Desktop/RADseq_YVBU/abundance_analysis")

# Clear environment and load packages ----

rm(list = ls()) # clean environment

# Packages

library(tidyverse)
library(terra)
library(sf)
library(landscapemetrics)
library(ggplot2)
library(AICcmodavg)
library(GGally)
library(stargazer)
library(patchwork)

# Load and explore data ----

all_counts <- read.csv("YVB_count.csv") %>% # Load all data
  filter(study == "after") %>%
  arrange(set) %>% 
  mutate(set_group = gsub("[0-9]", "", set)) # remove numbers from set 

# Print the count per site
site_counts <- all_counts %>% 
  group_by(site) %>%
  summarise(Total = sum(count)) 

# Load site data
sites <- read.csv("YVB_count_sites.csv")  %>%
  rename(site = X2010.names) %>%
  filter(X2004.names != "Khatib Bongsu") %>%
  arrange(site)

# load raster data
veg5 <- rast("veg5_merge78.tif") %>% # Landcover from Hansen et al.
  as.factor()
classes <- data.frame(class = c(0,3,4,5,7),
                      cover = c("Water","Urban",
                                "Urban green","Scrub",
                                "Forest"))
levels(veg5) <- classes
plot(veg5,col = c("white","gray55","khaki1","indianred1","aquamarine3")) 

# Create spatial object for sampling sites
sites.sf <- st_as_sf(sites, coords = c("Long", "Lat"), crs = 4326)
sites.utm <- st_transform(sites.sf, crs(veg5))
points(sites.utm,pch=23,bg="blue3",col="white",cex=2)

# Format data ----

# zero fill the absent count data
groups <- all_counts %>% 
  dplyr::select(site,set_group) %>%
  unique() 
groups[29,] <- c("Boon Lay Ave","boonlay")

# Get abundnace covariates

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

sites_buffer <- cbind(sites,buffer) %>% 
  left_join(site_counts,by="site")
sites_buffer[is.na(sites_buffer)] <- 0

# check distributions and correlation of predictors
buffer %>% 
  dplyr::select(-plot_id) %>%
  ggpairs()

# # I'll remove water as it has a ridiculous distribution. Secondary looks to be problematic with multicollinearity
# 
# # Set up the layout: 1 row and 3 columns (two for plots and one for the legend)
# layout(matrix(c(1, 2), nrow = 1), widths = c(3, 2))
# 
# # Plot the first graph
# par(mar = c(5, 5, 4, 3))  # Adjust margins
# plot(veg5, col=c("white","azure3","tan3","darkolivegreen1","darkolivegreen4","darkgreen"))
# points(sites.utm, cex =sites_buffer$Total/70 + 0.5, pch = 21, bg = "blue")
# 
# # add legend
# par(mar = c(5, 1, 4, 1)) 
# plot.new()
# legend_sizes <- c(min(sites_buffer$Total/70),max(sites_buffer$Total/70)) + 0.5
# legend("top", 
#        legend = round((legend_sizes - 0.5)*70), 
#        pch = 21, 
#        cex = 0.7,
#        pt.cex = legend_sizes, 
#        pt.bg = "blue", 
#        title = "Avg. count")

# Check relationships between predictors and abundance
par(mar = c(5, 5, 5, 5),mfcol=c(2,3)) 
plot(Total ~ Urban, dat = sites_buffer,
     xlab = "urban",ylab="Mean abundance")
plot(Total ~ Urban_green, dat = sites_buffer,
     xlab = "urban_green",ylab="Mean abundance")
plot(Total ~ Scrub, dat = sites_buffer,
     xlab = "scrub",ylab="Mean abundance")
plot(Total ~ secondary, dat = sites_buffer,
     xlab = "secondary",ylab="Mean abundance")
plot(Total ~ Lat, dat = sites_buffer,
     xlab = "Latitude",ylab="Mean abundance")

# Fit linear model ----

lm_global <- lm(Total ~ scale(Urban) + scale(`Urban green`) + scale(Scrub) + scale(Forest),dat = sites_buffer)

# Look at global model
summary(lm_global)

plot(lm_global)

# resids dont look good

# Fit glms ----

# check poisson
glm_global <- glm(Total ~ scale(Urban) + scale(`Urban green`) + scale(Scrub) + scale(Forest), dat = sites_buffer)

deviance(glm_global)/df.residual(glm_global) 
# Overdispersed

# Now fit negative binomial glms
library(MASS)

modList <- list()

modList[["null"]] <- m1 <- glm.nb(Total ~ 1, dat = sites_buffer)

# modList[["global"]] <- m2 <- update(m1, formula = ~ scale(urban) + scale(urban_green) + scale(scrub) + scale(secondary))

modList[["urban"]] <- m3 <- update(m1, formula = ~ scale(Urban))

modList[["urban_green"]] <- m4 <- update(m1, formula = ~ scale(`Urban green`))

modList[["scrub"]] <- m5 <- update(m1, formula = ~ scale(Scrub))

modList[["secondary"]] <- m6 <- update(m1, formula = ~ scale(Forest))

modList[["urban + urban_green"]] <- m7 <- update(m1, formula = ~ scale(Urban) + scale(`Urban green`))

modList[["urban + scrub"]] <- m8 <- update(m1, formula = ~ scale(Urban) + scale(Scrub))

modList[["urban_green + secondary"]] <- m9 <- update(m1, formula = ~ scale(`Urban green`) + scale(Forest))

modList[["urban_green + scrub"]] <- m10 <- update(m1, formula = ~ scale(`Urban green`) + scale(Scrub))

modList[["urban_green + urban + scrub"]] <- m11 <- update(m1, formula = ~ scale(`Urban green`) + scale(Urban) + scale(Scrub))

modList[["Lat"]] <- m12 <- update(m1, formula = ~ scale(Lat))

aic <- aictab(modList)
aic

stargazer(aic, 
          title = "", 
          summary = F,
          type = "html",
          digits = 2,
          font.size = "normalsize",
          out = "myAICtable_abundance.doc")

# Might as well check model averaged results
MuMIn::model.avg(modList)

# Present results best model ----

bestMod <- m3
summary(bestMod)

# check explained deviance
1 - bestMod$deviance / bestMod$null.deviance 
# This model explains 12% of dviance

# urban

newData <- expand.grid(Urban = seq(min(sites_buffer$Urban), 
                                         max (sites_buffer$Urban), 
                                         length.out = 1000))

pred.link <- as.data.frame(predict(bestMod, 
                                      newdata = newData, 
                                      type = "link", 
                                      se = T)) 

# Construct 95% CI
pred.p <- data.frame(newData, pred.link) %>%
  mutate(
    response = exp(fit),
    lcl.response = exp(fit - 1.96*se.fit),
    ucl.response = exp(fit + 1.96*se.fit)
  ) 

# Plotting

(urban.plot <- ggplot(data = pred.p, 
                            mapping = aes(x = Urban, 
                                          y = response)) +
    geom_line() +
    geom_ribbon(aes(ymin = lcl.response, 
                    ymax = ucl.response), 
                alpha = 0.2) +
    geom_point(data = sites_buffer, 
               mapping = aes(x = Urban, 
                             y = Total), 
               pch = 21,
               size = 2,
               fill="gray60") +
    theme_test( ) +
    theme(text = element_text(size = 12)) +
    labs(title = "A.",
         x = "Urban (%)", 
         y = "Count") )

# Check model assumptions best model ----

# We can first check the mean of the residuals is zero
mean(residuals(bestMod))

# Next I'm plotting the residuals along with the fitted values for our data. 
plot(residuals(bestMod) ~ fitted(bestMod), 
     xlab = "Fitted total", 
     ylab = "Deviance residuals", 
     pch = 16, 
     col = "darkgrey", 
     cex = 1.2)
abline(h = 0)
lines(lowess(residuals(bestMod) ~  fitted(bestMod)), 
      col = "red")

# I can also look at our other predictor, but no pattern is evident here, which is good. Note that the red smoothed line is simply a visual guide. It does not have to be perfectly flat for you to conclude that things look good. Don't put too much faith in the red line itself.
plot(residuals(bestMod) ~ sites_buffer$Urban, 
     xlab = "Urban", 
     ylab = "Deviance residuals", 
     pch = 16, 
     col = "darkgrey", 
     cex = 1.2)
abline(h = 0)
lines(lowess(residuals(bestMod) ~  sites_buffer$Urban), 
      col = "red")

# What about unused predictors
plot(residuals(bestMod) ~ sites_buffer$`Urban green`, 
     xlab = "urban_green", 
     ylab = "Deviance residuals", 
     pch = 16, 
     col = "darkgrey", 
     cex = 1.2)
abline(h = 0)
lines(lowess(residuals(bestMod) ~  sites_buffer$`Urban green`), 
      col = "red")

plot(residuals(bestMod) ~ sites_buffer$Scrub, 
     xlab = "scrub", 
     ylab = "Deviance residuals", 
     pch = 16, 
     col = "darkgrey", 
     cex = 1.2)
abline(h = 0)
lines(lowess(residuals(bestMod) ~  sites_buffer$Scrub), 
      col = "red")

plot(residuals(bestMod) ~ sites_buffer$Forest, 
     xlab = "forest", 
     ylab = "Deviance residuals", 
     pch = 16, 
     col = "darkgrey", 
     cex = 1.2)
abline(h = 0)
lines(lowess(residuals(bestMod) ~  sites_buffer$Forest), 
      col = "red")

# Check for spatial autocorrelation

sites_buffer %>% 
  mutate(resid = bestMod$residuals,
         resid_abs = abs(resid),
         sign = if_else(resid >= 0, "pos",
                        "neg")) %>% 
  ggplot(aes(x = Long,
             y = Lat,
             size = resid_abs,
             col = sign)) + 
  geom_point()

# no real signs of spatial autocorrelation

# Present results second best model ----

bestMod <- m9
summary(bestMod)

# check explained deviance
1 - bestMod$deviance / bestMod$null.deviance 
# This model explains 12% of dviance

# urban green

newData <- expand.grid(`Urban green` = seq(min(sites_buffer$`Urban green`), 
                                         max (sites_buffer$`Urban green`), 
                                         length.out = 1000), 
                       Forest = mean(sites_buffer$Forest))

pred.link <- as.data.frame(predict(bestMod, 
                                   newdata = newData, 
                                   type = "link", 
                                   se = T)) 

# Construct 95% CI
pred.p <- data.frame(newData, pred.link) %>%
  mutate(
    response = exp(fit),
    lcl.response = exp(fit - 1.96*se.fit),
    ucl.response = exp(fit + 1.96*se.fit)
  ) 

# Plotting

(urban_green.plot <- ggplot(data = pred.p, 
                            mapping = aes(x = `Urban.green`, 
                                          y = response)) +
    geom_line() +
    geom_ribbon(aes(ymin = lcl.response, 
                    ymax = ucl.response), 
                alpha = 0.2) +
    geom_point(data = sites_buffer, 
               mapping = aes(x = `Urban green`, 
                             y = Total), 
               pch = 21, 
               size = 2,
               fill="khaki1") +
    ylim(0,250) + 
    theme_test( ) +
    theme(text = element_text(size = 12)) +
    labs(title = "C.",
         x = "Urban green (%)", 
         y = "Count") )

# secondary

newData <- expand.grid(Forest = seq(min(sites_buffer$Forest), 
                                             max (sites_buffer$Forest), 
                                             length.out = 1000), 
                       `Urban green` = mean(sites_buffer$`Urban green`))

pred.link <- as.data.frame(predict(bestMod, 
                                   newdata = newData, 
                                   type = "link", 
                                   se = T)) 

# Construct 95% CI
pred.p <- data.frame(newData, pred.link) %>%
  mutate(
    response = exp(fit),
    lcl.response = exp(fit - 1.96*se.fit),
    ucl.response = exp(fit + 1.96*se.fit)
  ) 

# Plotting

(secondary.plot <- ggplot(data = pred.p, 
                                mapping = aes(x = Forest, 
                                              y = response)) +
    geom_line() +
    geom_ribbon(aes(ymin = lcl.response, 
                    ymax = ucl.response), 
                alpha = 0.2) +
    geom_point(data = sites_buffer, 
               mapping = aes(x = Forest, 
                             y = Total), 
               pch = 21, 
               size = 2,
               fill = "mediumaquamarine") +
    theme_test( ) +
    theme(text = element_text(size = 12)) +
    labs(title = "D.",
         x = "Forest (%)", 
         y = " ") )

# Lastly add latitude
(lat.plot <- ggplot(data = sites_buffer, 
                      mapping = aes(x = Lat, 
                                    y = Total)) +
    geom_point(pch = 21, 
               size = 2,
               fill = "black") +
    theme_test( ) +
    theme(text = element_text(size = 12)) +
    labs(title = "B.",
         x = "Latitude (°)", 
         y = "") )

(combo_plot <- (urban.plot | lat.plot) / (urban_green.plot | secondary.plot) +
  plot_annotation(title = "Abundance") & 
  theme(text = element_text(size = 12)) )

ggsave("combo_plot_abundance.jpeg",
       device = "jpeg",
       width = 5,                     # Set width (in inches)
       height = 5.25,                     # Set height (in inches)
       dpi = 300                       # Set DPI for raster formats
)

save(urban.plot,
     lat.plot,
     urban_green.plot,
     secondary.plot,
     combo_plot, 
     file = "yvb_abund_plots.RData")

# Check model assumptions second best model ----

# We can first check the mean of the residuals is zero
mean(residuals(bestMod))

# Next I'm plotting the residuals along with the fitted values for our data. 
plot(residuals(bestMod) ~ fitted(bestMod), 
     xlab = "Fitted total", 
     ylab = "Deviance residuals", 
     pch = 16, 
     col = "darkgrey", 
     cex = 1.2)
abline(h = 0)
lines(lowess(residuals(bestMod) ~  fitted(bestMod)), 
      col = "red")

# Next I'm plotting the residuals along with our predictor variables
plot(residuals(bestMod) ~ sites_buffer$Forest, 
     xlab = "Secondary", 
     ylab = "Deviance residuals", 
     pch = 16, 
     col = "darkgrey", 
     cex = 1.2)
abline(h = 0)
lines(lowess(residuals(bestMod) ~  sites_buffer$Forest), 
      col = "red")

# I can also look at our other predictor, but no pattern is evident here, which is good. Note that the red smoothed line is simply a visual guide. It does not have to be perfectly flat for you to conclude that things look good. Don't put too much faith in the red line itself.
plot(residuals(bestMod) ~ sites_buffer$`Urban green`, 
     xlab = "Urban green", 
     ylab = "Deviance residuals", 
     pch = 16, 
     col = "darkgrey", 
     cex = 1.2)
abline(h = 0)
lines(lowess(residuals(bestMod) ~  sites_buffer$`Urban green`), 
      col = "red")


# Check for spatial autocorrelation

sites_buffer %>% 
  mutate(resid = bestMod$residuals,
         resid_abs = abs(resid),
         sign = if_else(resid >= 0, "pos",
                        "neg")) %>% 
  ggplot(aes(x = Long,
             y = Lat,
             size = resid_abs,
             col = sign)) + 
  geom_point()

# no real signs of spatial autocorrelation