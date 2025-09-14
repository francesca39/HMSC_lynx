## ------------------------------------------------------------------------
## 'Habitat complexity and prey composition shape an apex predator’s habitat use across contrasting landscapes'
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
# 'R script to reproduce the prediction maps from HMSC in Supplementary material S1

# 1. LOAD LIBRARIES
# ========================


library(Hmsc)# Hmsc is the core package used for Hierarchical Modelling of Species Communities
library(raster)# Used for handling and manipulating spatial raster data
library(ggplot2)# ggplot2 is a popular library for making plots and maps
library(dplyr)# dplyr is used for data wrangling and transformation

library(abind)# abind allows us to merge multidimensional arrays
library(readxl)# readxl is used to import Excel files (.xlsx)

# ========================
# 2. LOAD MODEL AND  DATA
# ========================
# Load the dataset  used to fit the model
setwd("D:/Datasets and code/datasets/Datasets_fullmodel")
XData<- readRDS("XData.R")

# Load a pre-trained Hmsc model from file
fullmodel<- load(".../fullmodel.R")

# Extract the actual model object from memory
m <- (fullmodel)




# ========================
# 3. PREPARE GRID FOR PREDICTION
# ========================

# Load focal raster data for which we want to make prediction map for example  "dead wood index" 
setwd("D:/Datasets and code/datasets/rasters for prediction maps_2km resolution")

dead_wood_index_raster <- raster("dead_wood_index_raster.tif")

# Convert raster into a dataframe format suitable for prediction
grid <- as.data.frame(dead_wood_index_raster, xy = TRUE)
colnames(grid)[3] <- "dead_wood_index"
grid <- na.omit(grid)  # remove missing values

# keep all other variables to the grid constant(0 = mean)
grid <- data.frame(grid,
                   snow_depth = 0, 
                   roads_distance = 0,
                   streams_distance = 0,
                   terrain_ruggedness = 0, 
                   covered_distance = 0,
                   forest_cover = 0,  # this will be overwritten
                   days_since_snow = 0,
                   residential_distance = 0,
                   streams_length = 0,
                   peatbogs = 0,
                   log_effort = 0)




# ========================
# 4. set null the random factors
# ========================

m$studyDesign <- NULL
m$ranLevels <- NULL
m$rL <- NULL
m$rLNames <- NULL
m$ranLevelsUsed <- NULL
m$nr <- 0



# ========================
# 5. BLOCK PREDICTIONS
# ========================

# To reduce memory use, we clear some internal model components not needed for prediction

# Split the grid into 300 chunks for block-wise predictions
grid_list <- split(grid, cut(1:nrow(grid), 300, labels = FALSE))
predictions_list <- vector("list", length = length(grid_list))

# Loop through each grid chunk and make predictions
for (i in seq_along(grid_list)) {
  cat("Processing chunk", i, "of", length(grid_list), "\n")
  foo <- grid_list[[i]]
  
  pred <- predict(m,
                  studyDesign = NULL,
                  XData = foo,
                  ranLevels = NULL,      # marginal prediction (no random effects)
                  expected = TRUE,
                  predictEtaMean = TRUE)
  
  # Average across MCMC samples and store
  predictions_list[[i]] <- apply(abind(pred, along = 3), c(1, 2), mean)
  gc()  # free memory
}

# ========================
# 6. COMBINE PREDICTIONS INTO FINAL GRID
# ========================

# Reconstruct full prediction and grid by binding all chunks
preds <- do.call(rbind, predictions_list)
result_grid <- do.call(rbind, grid_list)

# ========================
# 7. SAVE PREDICTIONS
# ========================

# Save the predictions list to disk for reuse
saveRDS(predictions_list,
        "predictions_dead_wood_index.rds")

# ========================
# 8. LOAD SAVED PREDICTIONS
# ========================

# Reload predictions if needed (e.g., in another session)
prediction.fullmodel <- readRDS("predictions_lists_dead_wood_index.rds")

# ========================
# 10. RECREATE GRID FOR COORDINATES
# ========================

# Use any raster to extract spatial coordinates
forest_raster <- raster(".../forest_cover_raster.tif")
grid <- as.data.frame(forest_raster, xy = TRUE)
grid <- na.omit(grid)

# ========================
# 9. VISUALIZE PREDICTIONS FOR A SELECTED SPECIES
# ========================

# Extract species names from the model
species_names <- colnames(m$Y)
print(species_names)

# Choose which species to plot according the order of columns in the dataset(e.g., index 3 for lynx )
lynx_tracks <- 3
selected_species_name <- species_names[lynx_tracks]
cat("Selected species:", selected_species_name, "\n")

# Extract and combine predictions for selected species
combined_predictions_lynx <- unlist(lapply(prediction.fullmodel, function(x) x[, lynx_tracks]))

# Add predictions to grid for plotting
grid$combined_predictions_lynx <- combined_predictions_lynx

# Plot raster map using ggplot2
ggplot(grid, aes(x = x, y = y, fill = (combined_predictions_lynx))) +
  geom_raster() +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    name = "Lynx track count (log)"
  ) +
  labs(title = "Prediction map") +
  coord_equal() +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )


###################################################################################

#BIVARIATE MAPS TO SEE OVERLAP OF HABITAT USE BETWEEN LYNX AND PREY, SUPPLEMENTARY MATERIAL S6

# ========================
# 1. LOAD LIBRARIES
# ========================

library(raster)         # For raster data handling
library(ggplot2)        # For data visualization
library(dplyr)          # For data manipulation
library(abind)          # For combining arrays
library(Hmsc)           # For ecological modeling (Hierarchical Modeling of Species Communities)
library(gridGraphics)   # For advanced graphics manipulation
library(grid)           # For layout and graphical objects

# ========================
# 2. LOAD THE TRAINED HMSC MODEL
# ========================

# Load a trained model from disk
fullmodel_100 <- load(".../fullmodel.R")

# Extract the actual model object into 'm'
m <- get(fullmodel_100)

# ========================
# 3. LOAD MULTIPLE PREDICTION FILES
# ========================

# Set working directory where prediction files are stored
setwd("...")

# Load prediction outputs for different environmental variables
prediction.nat2               <- readRDS("predictions_lists_dead_wood_index.rds")
prediction.forestcover        <- readRDS("predictions_lists_forestcover.rds")
prediction.roadsdistance      <- readRDS("predictions_lists_roads_distance.rds")
prediction.residentialdistance<- readRDS("predictions_lists_residential_distance.rds")
prediction.streamslength      <- readRDS("predictions_lists_streams_length.rds")
prediction.streamsdistance    <- readRDS("predictions_lists_streams_distance.rds")
prediction.snowdepth          <- readRDS("predictions_lists_snowdepth.rds")
prediction.tri                <- readRDS("predictions_lists_terrain_ruggedness.rds")
prediction.peatbogs           <- readRDS("predictions_lists_peatbogs.rds")


# ========================
# 4. COMBINE PREDICTIONS INTO DATAFRAMES
# ========================

# Each list of predictions is converted into a single dataframe
combined_snowdepth           <- bind_rows(lapply(prediction.snowdepth, as.data.frame))
combined_nat2                <- bind_rows(lapply(prediction.nat2, as.data.frame))
combined_ruggedness          <- bind_rows(lapply(prediction.ruggedness, as.data.frame))
combined_forestcover         <- bind_rows(lapply(prediction.forestcover, as.data.frame))
combined_tri                 <- bind_rows(lapply(prediction.tri, as.data.frame))
combined_roadsdistance       <- bind_rows(lapply(prediction.roadsdistance, as.data.frame))
combined_residentialdistance <- bind_rows(lapply(prediction.residentialdistance, as.data.frame))
combined_streamsdistance     <- bind_rows(lapply(prediction.streamsdistance, as.data.frame))
combined_streamslength       <- bind_rows(lapply(prediction.streamslength, as.data.frame))
combined_peatbogs            <- bind_rows(lapply(prediction.peatbogs, as.data.frame))

# ========================
# 5. SELECT SPECIES AND AVERAGE PREDICTIONS
# ========================

# Extract species names
species_names <- colnames(m$Y)

# Select the species index (e.g., Lynx = 3)
species_tracks <- 3
selected_species_name <- species_names[species_tracks]
print(paste("Selected species:", selected_species_name))

# Combine selected predictions into a list
pred_list <- list(combined_snowdepth,
                  combined_streamsdistance,
                  combined_nat2,
                  combined_forestcover,
                  combined_residentialdistance,
                  combined_peatbogs,
                  combined_roadsdistance,
                  combined_streamslength,
                  combined_tri)

# Average all predictions across variables
mean_pred <- Reduce("+", pred_list) / length(pred_list)

# ========================
# 6. LOAD GRID FOR SPATIAL COORDINATES
# ========================

#use any raster for making the grid
dead_wood_index_raster <- raster("dead_wood_index_raster.tif")

# Extract coordinates from raster
grid <- as.data.frame(dead_wood_index_raster, xy = TRUE)
grid <- na.omit(grid)

# Merge spatial grid with predicted values
final_pred <- cbind(grid, mean_pred)

# ========================
# 7. BIVARIATE MAPPING (LYNX VS PREY)
# ========================

library(biscale)  # For bivariate classification
library(cowplot)  # For plotting composite maps

# --- Example: Lynx vs White-tailed Deer ---
data <- bi_class(final_pred, x = white_tailed_deer_ntracks, y = lynx_ntracks, style = "quantile", dim = 3)

map <- data %>%
  ggplot(aes(x = x, y = y, fill = bi_class)) +
  geom_raster(show.legend = FALSE) +
  bi_scale_fill(pal = "Brown2", dim = 3) +
  theme_map() +
  coord_fixed()

legend <- bi_legend(pal = "Brown2", dim = 3,
                    xlab = "White-tailed deer", 
                    ylab = "Lynx", size = 15)

finalPlot <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, -0.02, .65, .3, .3) +
  draw_figure_label("Distribution of Lynx and White-tailed Deer", fontface = "bold", size = 16)

print(finalPlot)

# Repeat similar bivariate plots for:
# Roe deer
# Brown hare
# Mountain hare

# Just replace the X variable (e.g., roe_deer_ntracks) and label accordingly

# ========================
# 8. EXPORT BIVARIATE MAP AS RASTER
# ========================

# Convert the bivariate class into numeric codes for raster export
data <- data %>%
  mutate(bivariate_code = as.numeric(factor(bi_class)))

# Create raster from classified map
bivariate_raster <- rasterFromXYZ(data[, c("x", "y", "bivariate_code")])

# Export as GeoTIFF
writeRaster(bivariate_raster, 
            filename = ".../bivariate_map_lynx_white_tailed_deer.tif", 
            format = "GTiff", 
            overwrite = TRUE)


####################END OF THE CODE################################################
