
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










# The code below are used to perform the spatial coverage and generate masks
# ******************************************************************************
# ------- function for Environmental coverage analysis (06-2,3) -----------
# ******************************************************************************
# notes: need copy (or prepare) geodata before run codes below

get_mask <- function(sdata, scenario){
  predictors <- glc_layers() %>% unname
  sdata <- as.data.frame(sdata) # needed
  # random forest model
  set.seed(202)
  print(paste0("**********", i, "th random forest"))
  Sys.time() %>% print()
  mod <- train(x = sdata[,predictors], 
               y = sdata$Cmic,
               method = "rf",
               importance = TRUE,
               tuneGrid = expand.grid(mtry = c(2:4)), # length(predictors) or 2:6
               trControl = trainControl(method = "cv", 
                                        number = 20,
                                        p = 0.75,
                                        savePredictions = TRUE))
  
  mod # most often mtry = 2
  
  # RMSE + R2
  mod$results %>% as_tibble %>% filter(mtry == mod$bestTune %>% unlist) %>% select(RMSE, Rsquared) %>% print()
  
  # no LC
  vars <- glc_layers() %>% .[!. %in% c("land_cover")] %>% unname
  
  # stack all layers
  st <- stack(map(vars, ~glc_get_resamp(.x, 2013, "5-15")))
  names(st) <- vars
  
  # ~2 min
  grid <- as.data.frame(st, xy = TRUE, na.rm = TRUE)
  
  # maha sdata
  maha_cmic <- mahadist(dat = sdata, world = grid, vars = vars, xy = c("x", "y"))
  mahamask <- maha_cmic %>% filter(mahatype == "chisq > 0.975") %>% select(x, y)
  sdata <- mahadist(dat = sdata, world = grid, vars = vars, xy = c("x", "y"))
  
  ### RF model ----------------------------------------------------------------
  predictors <- mod$trainingData %>% names %>% .[-length(.)]
  
  # prediction layers as stack
  raspred <- stack(map(predictors, ~glc_get_resamp(.x, 2013, "5-15")))
  # cbind(predictors, names(raspred) %>% word(sep = "_"))
  names(raspred) <- predictors
  
  # ~2min
  print(paste0("**********", i, "th grid"))
  preddf <- as.data.frame(raspred, xy = TRUE, na.rm = TRUE)
  all(names(preddf %>% select(-c(x, y))) %in% names(mod$trainingData))
  
  # fix land_cover col
  preddf <- preddf %>% 
    mutate(land_cover = glc_LC_num2chr(land_cover)) %>% 
    filter(land_cover %in% levels(mod$trainingData$land_cover)) %>% 
    mutate(land_cover = factor(land_cover)) #remove unused levels
  
  # predict from stack.as.df (~1min)
  prediction <- predict(mod, preddf)
  
  ### aoa calculation ---------------------------------------------------------
  # with variable weighting: (~5-10 min)
  AOA <- aoa(preddf, model = mod) # previously with argument returnTrainDI, cl
  # AOA$AOA %>% table
  
  # saveRDS(AOA, here("derived", "06-2-AOA_object.rds"))
  
  ## put together ------------------------------------------------------------
  
  preds <- mod$pred[mod$pred$mtry==mod$bestTune$mtry,]
  
  absError <- abs(preds$pred[order(preds$rowIndex)]-preds$obs[order(preds$rowIndex)])
  
  preddf <- preddf %>% mutate(DI = AOA$DI,
                              AOA = AOA$AOA)
  
  ## Create mask
  aoamask <- preddf %>% filter(AOA == 0) %>% select(x, y)
  aoamask$aoam <- 1
  aoamask <- aoamask %>% 
    mutate(pid = glc_pIDfromXY(x, y)) %>% 
    select(-c(x,y))
  
  # mahalanobis mask
  mahamask
  mahamask$maham <- 1
  mahamask <- mahamask %>% 
    mutate(pid = glc_pIDfromXY(x, y)) %>% 
    select(-c(x,y))
  
  # combine masks
  fullmask <- full_join(aoamask, mahamask, by = "pid")
  
  fullmask <- fullmask %>% select(pid, maha_mask = maham, aoa_mask = aoam) %>% 
    mutate(mask = 1) # as.numeric(maha_mask | aoa_mask)
  
  print(paste0("**********", i, "th save full mask"))
  Sys.time() %>% print()
  
  # need update based on scenario
  if(!dir.exists("derived")){
    dir.create("derived")
    }
  if (!dir.exists(paste0("derived/", scenario))){
    dir.create(paste0("derived/", scenario))
  }
  saveRDS(fullmask, here(paste0("derived/", scenario, "/", "fullmask_df", i, ".rds"))) 
}


# ******************************************************************************
# ---------------- get mask for patoine_500 ----------------
# ******************************************************************************

# get mask for scenario mask_patoine_500 --------------------------------------
# get masks for every 200 resamples

outputs = read.csv("output/resampled_cmic_500_patoine.csv") %>% # need update
  mutate(land_cover = factor(land_cover, levels = c("C", "FB", "FC", "FT", "G", "S"))) %>% 
  mutate(tmean = (tmean+273)*10)

outputs$run_number %>% unique() -> sample_number
cl <- makeCluster(6) #8-10
registerDoParallel(cl)

for(i in sample_number){
  sdata <- outputs %>% filter(run_number == i)
  get_mask(sdata, "mask_patoine_500")
}

stopCluster(cl)


# ---------------- get Intersection of 100 masks ----------------

mask_file <- list.files(path = "derived/mask_patoine_500/")

fullmask <- readRDS(paste0("derived/mask_patoine_500/", mask_file[1])) %>% 
  select(pid, mask)

for(i in 2:length(mask_file)){
  fullmask_i <- readRDS(paste0("derived/mask_patoine_500/", mask_file[i])) %>% 
    select(pid, mask)
  
  fullmask %>% full_join(fullmask_i, by = c("pid", "mask")) -> fullmask
  print(paste0(i, "th run ****************"))
}

# need update
saveRDS(fullmask, here(paste0("derived/mask_patoine_500/", "fullmask_fjoin", ".rds"))) 





# ******************************************************************************
# ---------------- get mask for combine_500 ----------------
# ******************************************************************************
# get mask for scenario mask_combine_500 --------------------------------------
# get masks for every 200 resamples

outputs = read.csv("output/resampled_cmic_500_combined.csv") %>% # need update
  mutate(land_cover = factor(land_cover, levels = c("C", "FB", "FC", "FT", "G", "S"))) %>% 
  mutate(tmean = (tmean+273)*10)

outputs$run_number %>% unique() -> sample_number
cl <- makeCluster(6) #8-10
registerDoParallel(cl)

for(i in sample_number){
  sdata <- outputs %>% filter(run_number == i)
  get_mask(sdata, "mask_combine_500")
}

stopCluster(cl)


# ---------------- get Intersection of 100 masks ----------------

mask_file <- list.files(path = "derived/mask_combine_500/")

fullmask <- readRDS(paste0("derived/mask_combine_500/", mask_file[1])) %>% 
  select(pid, mask)

for(i in 2:length(mask_file)){
  fullmask_i <- readRDS(paste0("derived/mask_combine_500/", mask_file[i])) %>% 
    select(pid, mask)
  
  fullmask %>% full_join(fullmask_i, by = c("pid", "mask")) -> fullmask
  print(paste0(i, "th run ****************"))
}

# need update
saveRDS(fullmask, here(paste0("derived/mask_combine_500/", "fullmask_fjoin", ".rds"))) 




# ******************************************************************************
# ---------------- plot mask of fjoin ----------------
# ******************************************************************************
### plot patoine --------------------------------------------------------------------
fullmask <- readRDS("derived/mask_patoine_500/fullmask_fjoin.rds")
preddf <- readRDS("derived/resampe500_combine/resample1_pred_all_vars/resample1_all_1992.rds")

preddf %>% dplyr::select(-mask) -> preddf


preddf <- preddf %>% mutate(pid = glc_pIDfromXY(x, y)) %>% 
  left_join(fullmask, by = "pid")

preddf %>% names

sum(preddf$mask == 1, na.rm = TRUE)

# Compare both
preddf <- preddf %>% mutate(Category = case_when(
  mask == 1 ~ "Area of outlier",
  TRUE ~ "Predicted") %>% factor(levels = c("Area of outlier", "Predicted", "Deficient data")))

preddf %>% count(Category)

cmic <- glc_proc_data()


# need to separate legend
(p_mask <- ggplot(preddf, aes(x, y, fill = Category))+
    borders(size = 0.3, fill = "grey90", ylim = c(-60, 90), colour = NA,
            show.legend = TRUE)+
    geom_raster(alpha = 0.7)+
    scale_fill_manual(values = c("Area of outlier" = "#cc78df",
                                 "Predicted" = "#8ef284",
                                 "Deficient data" = "grey90"),
                      drop = FALSE)+
    theme_void()+
    theme(
      legend.position = c(0.2, 0.25),
      legend.title = element_text(face = "bold"))+
    coord_fixed())

gg_width = 11
gg_height = 5.7

ggsave(plot = p_mask,
       here("output/figures", "Figure_X4_predicted_patoine_fjoin.png"),
       width = gg_width, height = gg_height)


### plot combine --------------------------------------------------------------------
fullmask <- readRDS("derived/mask_combine_500/fullmask_fjoin.rds")
preddf <- readRDS("derived/resampe500_combine/resample1_pred_all_vars/resample1_all_1992.rds")

preddf %>% select(-mask) -> preddf


preddf <- preddf %>% mutate(pid = glc_pIDfromXY(x, y)) %>% 
  left_join(fullmask, by = "pid")

preddf %>% names

sum(preddf$mask == 1, na.rm = TRUE)

# Compare both
preddf <- preddf %>% mutate(Category = case_when(
  mask == 1 ~ "Area of outlier",
  TRUE ~ "Predicted") %>%
    factor(levels = c("Area of outlier", "Predicted", "Deficient data")))

preddf %>% count(Category)

cmic <- glc_proc_data()


# need to separate legend
(p_mask <- ggplot(preddf, aes(x, y, fill = Category))+
    borders(size = 0.3, fill = "grey90", ylim = c(-60, 90), colour = NA,
            show.legend = TRUE)+
    geom_raster(alpha = 0.7)+
    scale_fill_manual(values = c("Area of outlier" = "#cc78df",
                                 "Predicted" = "#8ef284",
                                 "Deficient data" = "grey90"),
                      drop = FALSE)+
    theme_void()+
    theme(
      legend.position = c(0.2, 0.25),
      legend.title = element_text(face = "bold"))+
    coord_fixed())

gg_width = 11
gg_height = 5.7

ggsave(plot = p_mask,
       here("output/figures", "Figure_X4_predicted_combine_fjoin.png"),
       width = gg_width, height = gg_height)




### plot mask ------------------------------------------------------------------
## creat function
### plot mask ------------------------------------------------------------------

preddf <- readRDS("derived/resampe500_combine/resample1_pred_all_vars/resample1_all_1992.rds") %>% 
  mutate(Category = "Predicted")

preddf2 <- readRDS("derived/resampe500_combine/resample2_pred_all_vars/resample2_all_1992.rds") %>% 
  mutate(Category = "Predicted")

plot_mask <- function (sdata) {
  p_mask <- ggplot(sdata, aes(x, y, fill = Category))+
    borders(size = 0.3, fill = "grey90", ylim = c(-60, 90), colour = NA,
            show.legend = TRUE)+
    geom_raster(alpha = 0.7)+
    scale_fill_manual(values = c("Predicted" = "#8ef284"),
                      drop = FALSE)+
    theme_void()+
    theme(legend.position = c(0.2, 0.25),
          legend.title = element_text(face = "bold"))+
    coord_fixed()
  p_mask
}

plot_mask(preddf)
plot_mask(preddf2)






