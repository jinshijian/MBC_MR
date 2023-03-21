
library(raster)
library(tools)
library(sf)
library(tidyverse)
library(here)
library(doParallel)
library(foreach)
library(caret)
library(CAST)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(forcats)
library(magick)
library(biscale)
library(gt)
# library(glcpck)
source(here::here("code/functions.R"))

# ******************************************************************************
# load dataset 
# ******************************************************************************
cmic0 <- glc_proc_data() %>% mutate(tmean = tmean/10-273)
cmic1 <- read.csv("rawdata/mbc_meta_final2.csv") %>% 
  mutate(land_cover = factor(land_cover, levels = c("C", "FB", "FC", "FT", "G", "S"))) %>% 
  mutate(tmean = tmean/10-273) %>% 
  mutate(latitude = Latitude,
         longitude = Longitude) -> cmic1

predictors <- c("Cmic","latitude","longitude", glc_layers()) %>% unname
cmic2 <- bind_rows(
  cmic0[,predictors],
  cmic1[,predictors])  
cmic <- as.data.frame(cmic2) # needed

lm(Cmic ~ tmean, data = cmic0) %>% summary()
lm(Cmic ~ tmean, data = cmic1) %>% summary()

summary(cmic)

## create function for resample ---------------------------
sample_cmic0_lm = function(sdata){
  sdata[sample(nrow(sdata), 500, replace=TRUE), ] -> sub_cmic 
  lm_cmic = lm(Cmic ~ tmean, data = sub_cmic) %>% summary()
  lm_cmic$coefficients[2,1] -> lm_slope
  sub_cmic$Slope = lm_slope
  return(sub_cmic)
}

# ******************************************************************************
# ---------- random sample and building linear model ---------
# --scenario 1: randomly resample 100 times (n=500) from Patoine (n=762) 
# ******************************************************************************
cmic0 <- glc_proc_data() %>% mutate(tmean = tmean/10-273)
predictors <- c("Cmic","latitude","longitude", glc_layers()) %>% unname
lm(Cmic ~ tmean, data = cmic0) %>% summary()

predictors <- c("Cmic","latitude","longitude", glc_layers()) %>% unname
cmic2 <- cmic0[,predictors]  
cmic <- as.data.frame(cmic2) # needed

## resample 200 times -------------------------------------
outputs <- cmic2[0,]

# set.seed(123456)
for(i in 1:200){
  sample_cmic0_lm(cmic) -> sub_cmic
  run_number <- tibble(run_number = rep(i, 500))
  sub_cmic <- bind_cols(sub_cmic, run_number)
  outputs <- bind_rows(outputs, sub_cmic)
  print(paste0("----------",i))
}

summary(outputs)

write.csv(outputs,"output/resampled_cmic_500_patoine.csv", row.names = F)


# ******************************************************************************
# ---------- random sample and building linear model ---------
# --scenario 2: randomly resample 100 times (n=500) from combined
# data of Patoine (n=762) and control of warming experiment (n=106)
# ******************************************************************************
cmic0 <- glc_proc_data() %>% mutate(tmean = tmean/10-273)
cmic1 <- read.csv("rawdata/mbc_meta_final2.csv") %>% 
  mutate(land_cover = factor(land_cover, levels = c("C", "FB", "FC", "FT", "G", "S"))) %>% 
  mutate(tmean = tmean/10-273) %>% 
  mutate(latitude = Latitude,
         longitude = Longitude) -> cmic1

predictors <- c("Cmic","latitude","longitude", glc_layers()) %>% unname
cmic2 <- bind_rows(
  cmic0[,predictors],
  cmic1[,predictors])  
cmic <- as.data.frame(cmic2) # needed
cmic %>% summary()

# test lm between MAT and MBC
lm(cmic0$Cmic ~ cmic0$tmean) %>% summary()
lm(cmic1$Cmic ~ cmic1$tmean) %>% summary()
lm(cmic$Cmic ~ cmic$tmean) %>% summary()


## resample 200 times ---------------------------------------------------------
outputs <- cmic2[0,]

# set.seed(123456)
for(i in 1:200){
  sample_cmic0_lm(cmic) -> sub_cmic
  run_number <- tibble(run_number = rep(i, 500))
  sub_cmic <- bind_cols(sub_cmic, run_number)
  outputs <- bind_rows(outputs, sub_cmic)
  print(paste0("----------",i))
}

summary(outputs)
write.csv(outputs,"output/resampled_cmic_500_combined.csv", row.names = F)



