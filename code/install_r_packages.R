
#*****************************************************************************************
# Use old version caret and CAST to avoid error:
packageHistory(package = "CAST", check.package = TRUE)

remove.packages('CAST')
remove.packages('caret')

install.packages('devtools')
require(devtools)
install_version("CAST", version = "0.5.1", repos = "http://cran.us.r-project.org")
install_version("caret", version = "6.0-90", repos = "http://cran.us.r-project.org")

# then need to restart R
#*****************************************************************************************

# date: 2022-06-15

# Load packages -----------------------------------------------------------
# skip this if all packages needed have been installed
install.packages('raster')
install.packages('tools')
install.packages('sf')
install.packages('tidyverse')
install.packages('here')
install.packages('doParallel')
install.packages('foreach')
install.packages('caret')
install.packages('CAST')
install.packages('cowplot')
install.packages('grid')
install.packages('gridExtra')
install.packages('ggExtra')
install.packages('forcats')
install.packages('magick')
install.packages('biscale')
install.packages('gt')
install.packages('readxl')
install.packages('fasterize')
install.packages('randomForest')
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
