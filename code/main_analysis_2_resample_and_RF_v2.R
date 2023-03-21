
#*****************************************************************************************
# Use old version caret and CAST to avoid error:
# install.packages("packageHistory")
# library(packageHistory)
# packageHistory(package = "CAST", check.package = TRUE)

# remove.packages('CAST')
# remove.packages('caret')
# require(devtools)
# install_version("CAST", version = "0.5.1", repos = "http://cran.us.r-project.org")
# install_version("caret", version = "6.0-90", repos = "http://cran.us.r-project.org")

# need restart R
#*****************************************************************************************

# date: 2022-06-15
# author: jinshi jian <jinshi@vt.edu>
# adapted from Guillaume Patoine <guillaume.patoine@idiv.de>
# description: This is the script related to the publication "Drivers and trends 
# of global soil microbial carbon over two decades".

# The project uses the following folder structure
#
# project_root
# ├─ code 
# ├─ derived 
# ├─ geodata
# ├─ output
# │   ├─ figures
# │   └─ tables
# └─ rawdata

# NOTE some of the functions require additional datasets that are not provided
# in this repository, dur to space limitation. For example, `glc_get_resamp()`
# is a convenience function to load a raster layer, based on the variable name,
# year and soil depth. Code sections based on that function won't work, unless
# these layers are available.

# Load packages -----------------------------------------------------------

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

# The packages below need to be installed, but do no need to be loaded.
# The needed functions are called directly using ::
# library(readxl)
# library(fasterize)


# source functions from other script --------------------------
source(here::here("code/functions.R"))

# figure dimensions and default theme
gg_width = 11
gg_height = 5.7
ggplot2::theme_set(ggplot2::theme_bw())


# scenario 1: randomly resample 100 times (n=500) from patoine (n=762)
# scenario 2: randomly resample 100 times (n=500) from combine (n=762+106)

scenario_file <- c("resampled_cmic_500_patoine.csv", "resampled_cmic_500_combined.csv")
folder_file <- c("derived/resampe500_patoine_all", "derived/resampe500_combine")
mask_file <- c("derived/mask_patoine_500/", "derived/mask_combine_500/")

# *************************************************************************************
outputs = read.csv(paste0("output/", scenario_file[scenario_number])) %>%
  mutate(land_cover = factor(land_cover, levels = c("C", "FB", "FC", "FT", "G", "S"))) %>% 
  mutate(tmean = (tmean+273)*10)

# need update
folder_name <- paste0(folder_file[scenario_number], "/resample")
folder_fjoin <- c("derived/resampe500_patoine_fjoin", "derived/resampe500_combine_fjoin") 

# need to update the data and output folder 
scenario_number = 1 # need update based on scenario

# *************************************************************************************

outputs$run_number %>% unique() -> output_number
# non-fixed dynamic predictors
all_dp <- c("") # changed by jian
fclim_dp <- c("tmean", "prec") 
fLC_dp <- c("ndvi", "land_cover") # original condition
# fLC_temp <- c("ndvi", "land_cover", "prec") # changed by jian


# par setup 
cl <- makeCluster(6, outfile = "")
registerDoParallel(cl)

# for (i in 1){
for (i in 1:output_number[1:2]){
  predictors <- glc_layers() %>% unname
  cmic = outputs %>% filter(run_number == i)
  set.seed(202)
  
  model <- train(x = cmic[,predictors], 
                 y = cmic$Cmic,
                 method = "rf",
                 importance = TRUE,
                 tuneGrid = expand.grid(mtry = c(2:4)), # length(predictors) or 2:6
                 trControl = trainControl(method = "cv", 
                                          number = 20,
                                          p = 0.75,
                                          savePredictions = TRUE))
  
  model # most often mtry = 2
  
  # RMSE + R2
  model$results %>% as_tibble %>% filter(mtry == model$bestTune %>% unlist) %>% select(RMSE, Rsquared) %>% print
  # variable importance
  # varImp(model) %>% plot
  # varImp(model, scale = FALSE) %>% plot
  
  # ************************************************************
  # ------- Fixed global change drivers (09-3 + 09-4) -----------
  # ************************************************************
  
  # run predictions with variables fixed one by one
  # based on 06-0-RF_and_predictions_v2
  
  # Predictions all years ---------------------------------------------------
  
  lc_levs <- levels(model$trainingData$land_cover)
  
  # predict, save rds, (raster), png
  # takes ~ 50min
  for (iyear in 1992:2013) {
    # iyear <- 1992
    
    fixvars <- all_dp
    
    # either use iyear or starting year 1992
    st <- stack(purrr::map(predictors, ~glc_get_resamp(.x, iyear, "5-15")))
    
    # rename
    names(st) <- predictors
    
    # as.data.frame (~ 30 sec, ~ 2min with na.rm = T)
    gridyear <- as.data.frame(st, xy = TRUE, na.rm = T) %>%
      mutate(pid = glc_pIDfromXY(x,y))
    
    # land_cover to text
    gridyear <- gridyear %>% 
      mutate(land_cover = glc_LC_num2chr(land_cover)) %>% 
      filter(land_cover %in% lc_levs) %>% 
      mutate(land_cover = factor(land_cover))
    
    # mask: update for dynamic mask ---------------------------------------updte
    fullmask <- readRDS(paste0(mask_file[scenario_number], "fullmask_df", i, ".rds")) %>% select(pid, mask)
    
    gridyear1 <- gridyear %>% left_join(fullmask, by = "pid") %>% 
      filter(is.na(mask))  # %>% select(-mask) # include mask does not influence the result
    
    # predict (100 sec) 
    gridyear1$pred <- predict(model, gridyear1)
    
    
    if(!dir.exists((folder_file[scenario_number]))){
      dir.create(folder_file[scenario_number])
      }
    
    if (!dir.exists(here(paste0(folder_name, i ,'_pred_all_vars')))) {
      dir.create(here(paste0(folder_name, i ,'_pred_all_vars')))
      }
    
    saveRDS(gridyear1, here(paste0(folder_name, i ,'_pred_all_vars'),
                            paste0("resample", i ,"_all_", iyear, ".rds")))
    
    
    # using same mask (fulljoin of all 200 runs) -------------------------update
    fullmask_fjoin <- readRDS(paste0(mask_file[scenario_number], "fullmask_fjoin", ".rds")) %>% select(pid, mask)
    
    gridyear_fjoin <- gridyear %>% left_join(fullmask_fjoin, by = "pid") %>% 
      filter(is.na(mask))  # %>% select(-mask) # include mask column does not influence the results
    
    # predict for fjoin (100 sec)
    gridyear_fjoin$pred <- predict(model, gridyear_fjoin)
    if(!dir.exists((folder_fjoin[scenario_number]))){
      dir.create(folder_fjoin[scenario_number])
      }
    
    if (!dir.exists(here(paste0(folder_fjoin[scenario_number], "/resample", i ,'_pred_all_vars')))) {
      dir.create(here(paste0(folder_fjoin[scenario_number], "/resample", i ,'_pred_all_vars')))}
    
    saveRDS(gridyear_fjoin, here(paste0(folder_fjoin[scenario_number], "/resample", i ,'_pred_all_vars'),
                                 paste0("resample", i ,"_all_", iyear, ".rds")))
    
    print(paste0("*************", iyear, "----------", Sys.time()))
    
  }
  print(paste0("resample*************", i, "----------", Sys.time()))
}

#end----------------------------------------------------------------------------
stopCluster(cl)



# ******************************************************************************
# for scenario 2, go to line 84, and change scenario_number to 2
# ******************************************************************************


