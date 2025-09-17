## ------------------------------------------------------------------------
## 'Habitat complexity and prey composition shape an apex predator’s habitat use across contrasting landscapes'
## ------------------------------------------------------------------------



## ------------------------------------------------------------------------
## R script to reproduce HMSC analyses.
## This code closely follows the protocol from the HMSC book:
## Ovaskainen, Otso, and Nerea Abrego. Joint species distribution modelling: 
## with applications in R. Cambridge University Press, 2020.
## ------------------------------------------------------------------------

# 1. LOAD REQUIRED LIBRARIES
library(readxl)       # For reading Excel files (if needed)
library(Hmsc)         # For Hierarchical Modelling of Species Communities
library(dplyr)        # For data manipulation
library(tidyverse)    # For general data wrangling and plotting

# 2. LOAD DATASETS (please ensure these .rds files are in your working directory)
setwd("D:/Datasets and code/datasets/Datasets_fullmodel")# example of working directory

# Load environmental covariates (predictor variables)
XData <- readRDS("XData.R")

# Load species response matrix (e.g., snow track counts)
YData <- readRDS("YData.R")

# Load spatial coordinates
xy <- readRDS("xy.R")

# Load study design with hierarchical grouping factors (e.g., Region, Year, Spatial Unit of the triangles)
studyDesign <- readRDS("studyDesign.R")

# 3. DEFINE RANDOM EFFECT STRUCTURE

# Spatial random effect using Nearest Neighbor Gaussian Process (NNGP)
rL1 <- HmscRandomLevel(sData = xy, sMethod = 'NNGP', nNeighbours = 10)

# Temporal random effect (by year)
rL2 <- HmscRandomLevel(units = unique(studyDesign$Year))

# Regional random effect
rL3 <- HmscRandomLevel(units = unique(studyDesign$Region))

# 4. SET MCMC PARAMETERS

set.seed(42)  # For reproducibility

thinning <- 100          # Thinning interval
samples <- 1000          # Number of posterior samples
burn_in <- 0.5 * samples * thinning  # Transient iterations (burn-in)
cores <- 10              # Number of parallel chains

# 5. DEFINE THE HMSC MODEL

m <- Hmsc(
  Y = as.matrix(YData), 
  XData = XData, 
  XFormula = ~ dead_wood_index + forest_cover + peatbogs +
    streams_distance + terrain_ruggedness +
    snow_depth + roads_distance + residential_distance + log_effort,
  studyDesign = data.frame(studyDesign),
  ranLevels = list("Triangles" = rL1, "Year" = rL2, "Region" = rL3),
  distr = "lognormal poisson"  # for count data
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

# Save the fitted model object for downstream analysis
save(full_model, file = paste0("full_model", thinning, ".RData"))# running time may still be long due to intensive MCMC sampling. 

## ------------------------------------------------------------------------
## End of script
## ------------------------------------------------------------------------
