# process radish and rga data straight on hydra

setwd("C:/Users/hrwil/OneDrive - George Mason University - O365 Production/Studies/yellow_vented_bulbul/out_hydra/drop_KB/")
library(tidyverse)
library(radish)
library(raster)


# assess RGA results from runs n < 10 (dropping KB)
setwd("rga_drop_KB/")
rd <- str_c(rep("rga_drop_veg5_merge78_", 10), cter(1:10), rep("/Results/"))
rd_allres <- str_c(rd, rep("All_Results_Table.csv"))
rd_catres <- str_c(rd, rep("CategoricalResults.csv"))

## assess AIC scores
allres <- lapply(as.list(rd_allres), read_csv)
allres2 <- allres[[1]]
n <- 2
for (i in allres[2:10]) {
  ij <- i %>% mutate(Surface=paste0(Surface, n))
  allres2 <- add_row(allres2, ij)
  n <- n+1
}
# stats vary a bit more for this than with 11 pops.
# LL ties between IBD and surface.
# AICc favors IBD


## assess resistance values
catres <- lapply(as.list(rd_catres), read_csv)
catres2 <- catres[[1]]
for (i in catres[2:10]) {
  catres2 <- full_join(catres2, i)
}
catres2 <- catres2[,3:13]
colnames(catres2) <- c("k", "AIC", "AICc", "R2m", "R2c", "LL", "Water", "urban", "urban_green", "scrub", "forest")
catres2


## create average row
catres3 <- catres2 %>%
  add_row(as_tibble(lapply(catres2, mean))) %>%
  add_row(as_tibble(lapply(catres2, sd))) %>%
  add_row(as_tibble(lapply(catres2, function(x) {sd(x)/mean(x)}))) %>%
  add_column("run"=c(1:10, "avg", "sd", "cv"), .before=1)
catres3

## create ranked columns
catres4 <- catres3 %>%
  dplyr::select(Water, urban, urban_green, scrub, forest) %>%
  t() %>% as_tibble() %>% dplyr::select(1:10) %>%
  mutate(across(, rank)) %>%
  t() %>% as_tibble
catres4.1 <- catres4 %>%
  add_row(as_tibble(lapply(catres4, mean))) %>%
  add_row(as_tibble(lapply(catres4, sd))) %>%
  add_row(as_tibble(lapply(catres4, function(x) {sd(x)/mean(x)}))) %>%
  add_column("run"=c(1:10, "avg", "sd", "cv"), .before=1) %>%
  right_join(catres3)
colnames(catres4.1)[2:6] <- c(paste0("ranked_", c("Water", "urban", "urban_green", "scrub", "forest")))
ranked_avgs_rga <- catres4.1 %>% filter(run=="avg") %>% dplyr::select(2:6) %>% unlist %>% round(digits=5) 


## RGA supp table
catres2 <- add_column(catres2, "Surface"=rep("Landcover", nrow(catres2)))

rga_supp <- allres2 %>%
  filter((str_detect(Surface, "Distance") | str_detect(Surface, "Null"))) %>%
  arrange(Surface) %>%
  dplyr::select(Surface, AIC, AICc, R2m, R2c, LL, k) %>%
  add_column(as_tibble(matrix(ncol=5, nrow=nrow(.)))) %>%
  setNames(c("Surface", "AIC", "AICc", "R2m", "R2c", "LL", "k", "Water", "urban", "urban_green", "scrub", "forest")) %>%
  add_row(catres2) %>%
  dplyr::select(-R2m, -R2c) %>%
  mutate("deltaAICc"=AICc-min(AICc), .before=4) %>%
  arrange(deltaAICc)
write_csv(rga_supp, file="rga_supplementary_table.csv")

# get radish info
setwd("..")
models <- list.files("radish_drop_KB/") %>%
  str_subset("model_")
models <- paste0("radish_drop_KB/", models)
models <- lapply(models, readRDS)
nulls <- list.files("radish_drop_KB/") %>%
  str_subset("modelnull")
nulls <- paste0("radish_drop_KB/", nulls)
nulls <- lapply(nulls, readRDS)

# make table of coefficients
coefs <- lapply(models, function(x) {as.numeric(
  c(x$aic, x$loglik, x$df, x$mle$theta))}) %>%
  as_tibble(, .name_repair = c("minimal")) %>%
  t() %>% 
  as_tibble(.name_repair = c("unique")) %>%
  add_column("aic_null"=unlist(lapply(nulls, function(x) {x$aic})), .before=2) %>%
  add_column("ll_null"=unlist(lapply(nulls, function(x) {x$loglik})), .before=4) %>%
  add_column("k_null"=unlist(lapply(nulls, function(x) {x$df})), .before=6) %>%
  add_column("water"=rep(0, length(models)), .before=7)
colnames(coefs) <- c("aic", "aic_null", "ll", "ll_null", "k", "k_null", "water", "urban", "urban_green", "scrub", "forest")


# add summary rows
coefs1 <- coefs %>%
  mutate(delta_aic = aic - aic_null, .before=3) %>%
  add_row(as_tibble(lapply(coefs, mean))) %>%
  add_row(as_tibble(lapply(coefs, sd))) %>%
  add_row(as_tibble(lapply(coefs, function(x) {sd(x)/mean(x)}))) %>%
  add_column("runs"=c(1:nrow(coefs), "mean", "sd", "cv"), .before=1)
View(coefs1)

# add ranked average cols
coefs2 <- coefs1[1:10,] %>%
  dplyr::select(water, urban, urban_green, scrub, forest) %>%
  rowwise() %>%
  mutate(ranks = list(rank(c(water, urban, urban_green, scrub, forest)))) %>%
  pull(ranks) %>%
  as_tibble(.name_repair="minimal") %>%
  t() %>%
  as_tibble %>%
  setNames(paste0("ranked_", c("water", "urban", "urban_green", "scrub", "forest")))
# add stats rows
coefs2.1 <- coefs2 %>%
  add_row(as_tibble(lapply(coefs2, mean))) %>%
  add_row(as_tibble(lapply(coefs2, sd))) %>%
  add_row(as_tibble(lapply(coefs2, function(x) {sd(x)/mean(x)}))) %>%
  add_column("runs"=c(1:10, "mean", "sd", "cv"), .before=1)
# produce final table
coefs3 <- left_join(coefs1, coefs2.1)

rad_ranked_avgs <- coefs3 %>% filter(runs == "mean") %>% dplyr::select(ranked_water, ranked_urban, ranked_urban_green,
                                                                ranked_scrub, ranked_forest) %>%
  unlist

## radish supp table
head(coefs)
nulls <- coefs %>% 
  dplyr::select(aic_null, ll_null, k_null) %>% 
  add_column(as_tibble(matrix(nrow=nrow(.), ncol=5))) %>%
  add_column("Surface"=rep("IBD", nrow(.)), .before=1) %>%
  setNames(c("Surface", "aic", "ll", "k", "water", "urban", "urban_green", "scrub", "forest"))
rad_supp <- coefs %>%
  dplyr::select(-aic_null, -ll_null, -k_null) %>%
  add_column("Surface"=rep("Landcover", nrow(.)), .before=1) %>%
  add_row(nulls) %>%
  mutate("deltaAIC"=aic-min(aic), .before=3) %>%
  arrange(deltaAIC)
write_csv(rad_supp, file="radish_supplementary_table.csv")

## visualizations

my_ggplot <- function(x, y, dat=coefs1[1:nrow(coefs),]) {
  ggplot(dat, aes_string(x, y))+
    geom_point() +
    ggtitle(paste0(x, "_v_", y))
}

aic_plots <- list()
for (i in c("urban", "urban_green", "scrub", "forest")) {
  print(i)
  aic_plots <- append(aic_plots, list(my_ggplot("aic", i)))
}
aic_plots

ll_plots <- list()
for (i in c("urban", "urban_green", "scrub", "forest")) {
  print(i)
  ll_plots <- append(ll_plots, list(my_ggplot("ll", i)))
}
ll_plots

# Avg RGA results: order low to high resistance: 
# urban, water, scrub, urban green, forest

# Avg radish: high to low conductance: 
# urban green, water, scrub, forest, urban


# get ranked avg layers
## RGA
rga_raw_lyr <- raster("rga_drop_KB/rga_drop_veg5_merge78_1/Results/veg5_merge78.asc")
current_vals <- catres3 %>% filter(run == 1) %>% 
  dplyr::select(Water, urban, urban_green, scrub, forest) %>% unlist %>%
  round(digits=5)

rga_ranked_avg_layer <- calc(rga_raw_lyr, function(x) {
  for (i in x) {
    j <- round(i, digits=5)
    return(ifelse(j %in% current_vals,
                  ranked_avgs_rga[j == current_vals],
                  NA))
  }
})

# might need adjusting to adequately represent gaps?
cols <- c(rev(terrain.colors(255))[seq(1, 135, 25)], rev(terrain.colors(255))[seq(135, 255, 10)])
plot(rga_ranked_avg_layer, colNA="red", col=cols)


## radish layer
radish_raw_lyr <- readRDS("radish_drop_KB/fitted_veg5_merge78_drop_KB_1.RDS")
radish_raw_lyr <- log(radish_raw_lyr[["est"]])
currents <- coefs1 %>% 
  filter(runs==2) %>% 
  dplyr::select(water, urban, urban_green, scrub, forest) %>%
  unlist

radish_ranked_avg_layer <- calc(radish_raw_lyr, function(x) {
  for (i in x) {
    return(ifelse(i %in% currents, 
                  rad_ranked_avgs[currents==i],
                  NA))
  }
})

rad_cols <- c(rev(terrain.colors(255))[seq(1, 90, 30)], 
          rev(terrain.colors(255))[seq(100, 175, 10)], 
          rev(terrain.colors(255))[seq(195, 255, 20)])
raster::plot(radish_ranked_avg_layer, col=cols, colNA="red")

rm(list=ls()[!(ls() %in% c("rad_cols", "radish_ranked_avg_layer",
               "cols", "rga_ranked_avg_layer",
               "coefs3", "catres4.1", "allres2"))])
save.image(file="dropKB_visualizations.RData")





