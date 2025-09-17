

## ------------------------------------------------------------------------
## HMSC script for Centre region
## ------------------------------------------------------------------------

# 1. LOAD REQUIRED LIBRARIES
library(readxl)
library(Hmsc)
library(dplyr)
library(tidyverse)
setwd("D:/Datasets and code/datasets/Datasets_centre_model")#set working directory 
# 2. LOAD DATASETS 
XData_centre <- readRDS("XData_centre.R")
YData_centre <- readRDS("YData_centre.R")
xy_centre <- readRDS("xy_centre.R")
studyDesign_centre <- readRDS("studyDesign_centre.R")

# 3. DEFINE RANDOM EFFECT STRUCTURE

# Spatial random effect
rL1 <- HmscRandomLevel(sData = xy_centre, sMethod = 'NNGP', nNeighbours = 10)

# Temporal random effect (by year)
rL2 <- HmscRandomLevel(units = unique(studyDesign_centre$Year))

# 4. SET MCMC PARAMETERS

set.seed(42)

# 6. RUN MCMC SAMPLING

thinning <- 100          # Thinning interval
samples <- 1000          # Number of posterior samples
burn_in <- 0.5 * samples * thinning  # Transient iterations (burn-in)
cores <- 10              # Number of parallel chains



# 5. DEFINE THE HMSC MODEL

m_centre <- Hmsc(
  Y = as.matrix(YData_centre),
  XData = XData_centre,
  XFormula = ~ dead_wood_index + forest_cover + peatbogs +
    streams_distance + streams_length + terrain_ruggedness +
    snow_depth + roads_distance + residential_distance + log_effort,
  studyDesign = data.frame(studyDesign_centre),
  ranLevels = list("Triangles" = rL1, "Year" = rL2),
  distr = "lognormal poisson"
)

# 6. RUN MCMC SAMPLING

centre_model <- sampleMcmc(
  m_centre,
  thin = thinning,
  samples = samples,
  transient = burn_in,
  nChains = cores,
  nParallel = cores,
  verbose = 1
)

# 7. SAVE THE OUTPUT

save(centre_model, file = paste0("centre_model", thinning, ".RData")) # running time may still be long due to intensive MCMC sampling. 

## ------------------------------------------------------------------------
## End of Centre region script
## ------------------------------------------------------------------------

