




## ------------------------------------------------------------------------
## HMSC script for South region
## ------------------------------------------------------------------------

# 1. LOAD REQUIRED LIBRARIES
library(readxl)
library(Hmsc)
library(dplyr)
library(tidyverse)


# 2. LOAD DATASETS (ensure these .rds files are in your working directory)

setwd("D:/Datasets and code/datasets/Datasets_south_model")#set working directory 

XData_south <- readRDS("XData_south.R")
YData_south <- readRDS("YData_south.R")
xy_south <- readRDS("xy_south.R")
studyDesign_south <- readRDS("studyDesign_south.R")

# 3. DEFINE RANDOM EFFECT STRUCTURE

# Spatial random effect
rL1 <- HmscRandomLevel(sData = xy_south, sMethod = 'NNGP', nNeighbours = 10)

# Temporal random effect (by year)
rL2 <- HmscRandomLevel(units = unique(studyDesign_south$Year))

# 4. SET MCMC PARAMETERS

set.seed(42)

thinning <- 100          # Thinning interval
samples <- 1000          # Number of posterior samples
burn_in <- 0.5 * samples * thinning  # Transient iterations (burn-in)
cores <- 10              # Number of parallel chains

# 5. DEFINE THE HMSC MODEL

m_south <- Hmsc(
  Y = as.matrix(YData_south),
  XData = XData_south,
  XFormula = ~ dead_wood_index + forest_cover + peatbogs +
    streams_distance + streams_length + terrain_ruggedness +
    snow_depth + roads_distance + residential_distance + log_effort,
  studyDesign = data.frame(studyDesign_south),
  ranLevels = list("Triangles" = rL1, "Year" = rL2),
  distr = "lognormal poisson"
)

# 6. RUN MCMC SAMPLING

south_model <- sampleMcmc(
  m_south,
  thin = thinning,
  samples = samples,
  transient = burn_in,
  nChains = cores,
  nParallel = cores,
  verbose = 1
)


# 7. SAVE THE OUTPUT

save(south_model, file = paste0("south_model", thinning, ".RData"))# running time may still be long due to intensive MCMC sampling. 

## ------------------------------------------------------------------------
## End of South region script
## ------------------------------------------------------------------------

