
#*******************************************************************************
# Calculate and compare temporal trends of all predictions -----------------
#*******************************************************************************
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
# install.packages('qdapRegex')
library(qdapRegex)
source(here::here("code/functions.R"))

# par setup 
cl <- makeCluster(6, outfile = "")
registerDoParallel(cl)

# calculate cstocks
bulk <- raster(here("geodata/resampled_0p05deg/static/bdod", "bdod_5-15cm_mean_resamp_c2021-02-14_153856.tif"))
cofr <- raster(here("geodata/resampled_0p05deg/static/cfvo", "cfvo_5-15cm_mean_resamp_c2021-02-14_163601.tif"))
rarea <- raster::area(cofr)

st <- stack(bulk, cofr, rarea)
names(st) <- c("bulk", "cofr", "area")

# Global

temporal_trend <- function(sdata, var_scenario) {
  # load all predictions (~ 2min)
  preds <- map_dfr(sdata, readRDS, .id = "year")
  
  # fix year
  preds$year <- as.numeric(preds$year) + 1991
  
  # ~ 2 min
  bulk_df <- as.data.frame(st, xy = TRUE, na.rm = TRUE)
  bulk_df <- bulk_df %>% mutate(pid = glc_pIDfromXY(x, y)) %>% 
    select(-c(x, y))
  
  preds <- preds %>% left_join(bulk_df, by = "pid")
  
  # stock is in tonnes/ha (weight/area)
  # cell_stock is the total weight for the cell in tonnes
  preds <- preds %>% mutate(stock = pred * bulk * (1000-cofr) * 12.01 / 10^8,
                            cell_stock = stock * area * 100) 
  
  # remove missing stocks
  preds <- preds %>% drop_na(stock)
  
  preds %>% names %>% writeLines
  
  # nest per cell --------------------------------------------------------
  prn <- preds %>% select(pid, year, cell_stock, prec, tmean, ndvi, land_cover) %>%
    group_by(pid,) %>%
    nest %>% ungroup
  
  # remove if missing temporal datapoints
  nrow(prn) #2658500
  prn <- prn %>% filter(map_lgl(data, ~ nrow(.x) == 22))
  nrow(prn) #2596897
  
  jreg <- prn %>% unnest(cols = c(data))
  
  # global temporal trend ---------------------------------------------------
  
  # NOTE cell_stock is in tonnes
  glob_sum <- jreg %>% group_by(year) %>% 
    summarise(cstock = sum(cell_stock))
  
  # mean value cstock
  mcs <- glob_sum$cstock %>% mean #4.34 Pg Cmic
  
  glob_sum <- glob_sum %>% mutate(change = cstock - cstock[1],
                                  perc_ch = change / cstock[1] * 100,
                                  ysca = year - 1992)
  
  glob_sum <- glob_sum %>% mutate(mod_fix = "full")
  
  # Calculation global cmstock changes ---------------------------------------
  lm(perc_ch ~ ysca, glob_sum) %>% summary() %>% print()
  
  lm(perc_ch ~ ysca, glob_sum) %>% summary() -> lm_model_sum
  lm_slope = lm_model_sum$coefficients[2,1]
  lm_slope_se = lm_model_sum$coefficients[2,2]
  lm_slope_p = lm_model_sum$coefficients[2,4]
  
  outpt_slope = tibble(scenario = var_scenario,
                       slope = lm_slope,
                       slope_se = lm_slope_se,
                       slope_p = lm_slope_p,
                       folder = sdata[1])
  
  return(outpt_slope)
}


# need update baseed on scenario
# note that vmask_output are outputs for runs with different mask
vmask_output <- c("derived/resampe500_patoine", "derived/resampe500_combine")
# fjoin_output are outputs for runs with uniform full joined mask
fjoin_output <- c("derived/resampe500_patoine_fjoin", "derived/resampe500_combine_fjoin")

scenario_file <- c("resampled_cmic_500_patoine.csv", "resampled_cmic_500_combined.csv")
folder_file <- c("derived/resampe500_patoine", "derived/resampe500_combine")






#*******************************************************************************
# scenario 1: patoine ---------------------------
#*******************************************************************************
scenario_number <- 1

# each resample run with diffirent mask ****************************************
list.files(here(vmask_output[scenario_number]),  
           pattern = c("resample"), full.names = TRUE) -> folder_list 

resample_lm_output = tibble(scenario = NA,
                            slope = NA,
                            slope_se = NA,
                            slope_p = NA,
                            folder = NA)

for (i in 1:length(folder_list)){
  list.files(folder_list[i], pattern = c("_all_.+rds"), full.names = TRUE) -> pred_files
  temporal_trend(pred_files, "fixed_no_var") -> resample_lm_output_i
  resample_lm_output = bind_rows(resample_lm_output, resample_lm_output_i)
  print(paste0("resample", "**************", i))
}

#end----------------------------------------------------------
stopCluster(cl)


# join and get run_number ---------------------------------------------

resample_lm_output %>% 
  filter(!is.na(slope)) %>% 
  ggplot(aes(x = slope)) +
  geom_histogram(bins = 30, color="black", fill="white")

mean(resample_lm_output$slope, na.rm = T)

resample_lm_output2 <- resample_lm_output %>% filter(!is.na(slope))
resample_lm_output2$folder 


outputs = read.csv(paste0("output/", scenario_file[scenario_number])) %>%  # need update
  mutate(land_cover = factor(land_cover, levels = c("C", "FB", "FC", "FT", "G", "S"))) %>% 
  mutate(tmean = (tmean+273)*10)

# tibble(run_number = c(1:100), outslope = outputs$Slope %>% unique()) -> optput_lm_slope

outputs %>% 
  mutate(outslope = Slope) %>% 
  select(run_number, outslope) %>% 
  unique() -> optput_lm_slope

# need update

rm_between(resample_lm_output2$folder, here(paste0(folder_file[scenario_number]),"resample"),
           '_pred_all_vars', extract=TRUE)

resample_lm_output2$run_number = unlist(rm_between(resample_lm_output2$folder,
                                                   here(paste0(folder_file[scenario_number]),"resample"),
                                                   '_pred_all_vars', extract=TRUE))

resample_lm_output2$run_number <- as.numeric(resample_lm_output2$run_number)

resample_lm_output2 %>%  
  left_join(optput_lm_slope, by = "run_number") -> resample_lm_output2

# need update based on scenario
# write.csv(resample_lm_output2, "output/resample_lm_output_500_patoine.csv", row.names = FALSE) 



# all resample runs with a same mask ******************************************
list.files(here(fjoin_output[scenario_number]),  
           pattern = c("resample"), full.names = TRUE) -> folder_list 

resample_lm_output = tibble(scenario = NA,
                            slope = NA,
                            slope_se = NA,
                            slope_p = NA,
                            folder = NA)

for (i in 1:length(folder_list)){
  list.files(folder_list[i], pattern = c("_all_.+rds"), full.names = TRUE) -> pred_files
  temporal_trend(pred_files, "fixed_no_var") -> resample_lm_output_i
  resample_lm_output = bind_rows(resample_lm_output, resample_lm_output_i)
  print(paste0("resample", "**************", i))
}

#end----------------------------------------------------------
stopCluster(cl)


# join and get run_number ---------------------------------------------

resample_lm_output %>% 
  filter(!is.na(slope)) %>% 
  ggplot(aes(x = slope)) +
  geom_histogram(bins = 30, color="black", fill="white")

mean(resample_lm_output$slope, na.rm = T)

resample_lm_output2 <- resample_lm_output %>% filter(!is.na(slope))
resample_lm_output2$folder 


outputs = read.csv(paste0("output/", scenario_file[scenario_number])) %>%  # need update
  mutate(land_cover = factor(land_cover, levels = c("C", "FB", "FC", "FT", "G", "S"))) %>% 
  mutate(tmean = (tmean+273)*10)

# tibble(run_number = c(1:100), outslope = outputs$Slope %>% unique()) -> optput_lm_slope

outputs %>% 
  mutate(outslope = Slope) %>% 
  select(run_number, outslope) %>% 
  unique() -> optput_lm_slope

# need update

rm_between(resample_lm_output2$folder, here(paste0(folder_file[scenario_number], "_fjoin"),"resample"),
           '_pred_all_vars', extract=TRUE)

resample_lm_output2$run_number = 
  unlist(rm_between(resample_lm_output2$folder,
                    here(paste0(folder_file[scenario_number], "_fjoin"),"resample"),
                    '_pred_all_vars', extract=TRUE))

resample_lm_output2$run_number <- as.numeric(resample_lm_output2$run_number)

resample_lm_output2 %>%  
  left_join(optput_lm_slope, by = "run_number") -> resample_lm_output2

# need update based on scenario
# write.csv(resample_lm_output2, "output/resample_lm_output_500_patoine_fjoin.csv", row.names = FALSE)







#*******************************************************************************
# scenario 2: combine ---------------------------
#*******************************************************************************
scenario_number <- 2

# each resample run with diffirent mask ****************************************
list.files(here(vmask_output[scenario_number]),  
           pattern = c("resample"), full.names = TRUE) -> folder_list 

resample_lm_output = tibble(scenario = NA,
                            slope = NA,
                            slope_se = NA,
                            slope_p = NA,
                            folder = NA)

for (i in 1:length(folder_list)){
  list.files(folder_list[i], pattern = c("_all_.+rds"), full.names = TRUE) -> pred_files
  temporal_trend(pred_files, "fixed_no_var") -> resample_lm_output_i
  resample_lm_output = bind_rows(resample_lm_output, resample_lm_output_i)
  print(paste0("resample", "**************", i))
}

#end----------------------------------------------------------
stopCluster(cl)


# join and get run_number ---------------------------------------------

resample_lm_output %>% 
  filter(!is.na(slope)) %>% 
  ggplot(aes(x = slope)) +
  geom_histogram(bins = 30, color="black", fill="white")

mean(resample_lm_output$slope, na.rm = T)

resample_lm_output2 <- resample_lm_output %>% filter(!is.na(slope))
resample_lm_output2$folder 


outputs = read.csv(paste0("output/", scenario_file[scenario_number])) %>%  # need update
  mutate(land_cover = factor(land_cover, levels = c("C", "FB", "FC", "FT", "G", "S"))) %>% 
  mutate(tmean = (tmean+273)*10)

# tibble(run_number = c(1:100), outslope = outputs$Slope %>% unique()) -> optput_lm_slope

outputs %>% 
  mutate(outslope = Slope) %>% 
  select(run_number, outslope) %>% 
  unique() -> optput_lm_slope

# need update

rm_between(resample_lm_output2$folder, here(paste0(folder_file[scenario_number]),"resample"),
           '_pred_all_vars', extract=TRUE)

resample_lm_output2$run_number = unlist(rm_between(resample_lm_output2$folder,
                                                   here(paste0(folder_file[scenario_number]),"resample"),
                                                   '_pred_all_vars', extract=TRUE))

resample_lm_output2$run_number <- as.numeric(resample_lm_output2$run_number)

resample_lm_output2 %>%  
  left_join(optput_lm_slope, by = "run_number") -> resample_lm_output2

# need update based on scenario
# write.csv(resample_lm_output2, "output/resample_lm_output_500_combine.csv", row.names = FALSE) 



# all resample runs with a same mask ******************************************
list.files(here(fjoin_output[scenario_number]),  
           pattern = c("resample"), full.names = TRUE) -> folder_list 

resample_lm_output = tibble(scenario = NA,
                            slope = NA,
                            slope_se = NA,
                            slope_p = NA,
                            folder = NA)

for (i in 1:length(folder_list)){
  list.files(folder_list[i], pattern = c("_all_.+rds"), full.names = TRUE) -> pred_files
  temporal_trend(pred_files, "fixed_no_var") -> resample_lm_output_i
  resample_lm_output = bind_rows(resample_lm_output, resample_lm_output_i)
  print(paste0("resample", "**************", i))
}

#end----------------------------------------------------------
stopCluster(cl)


# join and get run_number ---------------------------------------------

resample_lm_output %>% 
  filter(!is.na(slope)) %>% 
  ggplot(aes(x = slope)) +
  geom_histogram(bins = 30, color="black", fill="white")

mean(resample_lm_output$slope, na.rm = T)

resample_lm_output2 <- resample_lm_output %>% filter(!is.na(slope))
resample_lm_output2$folder 


outputs = read.csv(paste0("output/", scenario_file[scenario_number])) %>%  # need update
  mutate(land_cover = factor(land_cover, levels = c("C", "FB", "FC", "FT", "G", "S"))) %>% 
  mutate(tmean = (tmean+273)*10)

# tibble(run_number = c(1:100), outslope = outputs$Slope %>% unique()) -> optput_lm_slope

outputs %>% 
  mutate(outslope = Slope) %>% 
  select(run_number, outslope) %>% 
  unique() -> optput_lm_slope

# need update
rm_between(resample_lm_output2$folder, here(paste0(folder_file[scenario_number], "_fjoin"),"resample"),
           '_pred_all_vars', extract=TRUE)

resample_lm_output2$run_number = 
  unlist(rm_between(resample_lm_output2$folder,
                    here(paste0(folder_file[scenario_number], "_fjoin"),"resample"),
                    '_pred_all_vars', extract=TRUE))

resample_lm_output2$run_number <- as.numeric(resample_lm_output2$run_number)

resample_lm_output2 %>%  
  left_join(optput_lm_slope, by = "run_number") -> resample_lm_output2

# need update based on scenario
# write.csv(resample_lm_output2, "output/resample_lm_output_500_combine_fjoin.csv", row.names = FALSE)


