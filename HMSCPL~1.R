## ------------------------------------------------------------------------
## ## Habitat complexity and prey composition shape an apex predator’s habitat use across contrasting landscapes
## ------------------------------------------------------------------------



## ------------------------------------------------------------------------
# 'R script to reproduce the PLOTS from HMSC output 

# 1. LOAD LIBRARIES
# ========================


library(janitor)     # Data cleaning (e.g., clean column names)
library(dplyr)       # Data wrangling
library(ggplot2)     # Data visualization
library(corrplot)    # Correlation matrix plotting
library(coda)        # MCMC diagnostics
library(Hmsc)        # Hierarchical Modelling of Species Communities (HMSC)
library(tidyr)       # Tidy data transformation
library(tidyverse)   # Meta-package that includes dplyr, tidyr, ggplot2, etc.

# Load the fitted HMSC models for different regions
centre <- load(".../centre_model.R")     # Central region model
south <- load(".../south_model.R")       # Southern region model
fullmodel <- load(".../fullmodel.R")          # Full study area model

# Assign loaded models to specific variables
m <- get(fullmodel)
m1 <- get(centre_model)
m2 <- get(south_model)

############################################################################################


#SUPPLEMENTARY MATERIAL S1 

# Convert model to a coda object to extract MCMC diagnostics
mpost <- convertToCodaObject(m)

# Calculate effective sample size (ESS) for Beta parameters
effectiveSizes <- effectiveSize(mpost$Beta)

# Plot ESS and Gelman-Rubin diagnostic (PSRF) to assess convergence
par(mfrow = c(2, 2))
hist(effectiveSizes, main = "Effective Sample Size (ESS) for Beta")

gelmanDiagnostics <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
hist(gelmanDiagnostics, main = "Potential Scale Reduction Factor (PSRF) for Beta")


#################################################################################################
#1.  VARIANCE PARTITIONING PLOT 

# Compute variance partitioning per predictor group
VP <- computeVariancePartitioning(m,  
                                  group = attr(m$X[[1]], "assign"),
                                  groupnames = attr(terms(m$XFormula, data = m$XData), "term.labels"))


# Rename predictors for plotting
predictor_names <- c(
  "log_effort" = "Sampling effort",
  "dead_wood_index" = "Dead Wood Potential",
  "forest_cover" = "Forest cover proportion",
  "residential_distance" = "Distance to residential areas (m)",
  "roads_distance" = "Distance to roads (m)",
  "snow_depth" = "Snow depth (cm)",
  "peatbogs" = "Peatbogs cover proportion",
  "streams_distance" = "Distance to streams (m)",
  "streams_length" = "Streams length (m)",
  "terrain_ruggedness" = "Terrain ruggedness",
  "Random: Triangles" = "Spatial random effect",
  "Random: Region" = "Region random effect",
  "Random: Year" = "Year random effect"
)

# Rename species for plotting (ALL SPECIES SUPPLEMENTARY MATERIAL S1, FIG.2)
species_names <- c(
  "wolf_ntracks" = "Wolf",
  "white_tailed_deer_ntracks" = "White-tailed deer",
  "roe_deer_ntracks" = "Roe deer",
  "red_fox_ntracks" = "Red fox",
  "mountain_hare_ntracks" = "Mountain hare",
  "moose_ntracks" = "Moose",
  "lynx_ntracks" = "Lynx",
  "forest_reindeer_ntracks" = "Forest reindeer",
  "brown_hare_ntracks" = "Brown hare"
)

colors <- c(
  "Sampling effort" = "#A89B8C",
  "Dead Wood Potential" = "#C05A50",
  "Forest cover proportion" = "#2d6a4f",
  "Distance to residential areas (m)" = "#fcbf49",
  "Distance to roads (m)" = "darkslateblue",
  "Snow depth (cm)" = "#B0C4DE",
  "Peatbogs cover proportion" = "mediumseagreen",
  "Distance to streams (m)" = "royalblue1",
  "Streams length (m)" = "powderblue",
  "Terrain ruggedness" = "chocolate4",
  "Region random effect" = "lightsalmon3",
  "Spatial random effect" = "#e9d8a6",
  "Year random effect" = "#939F88"
)
########################## SELECTED SPECIES MANUSCRIPT FIG. 2#######################################################
selected_species <- c("Lynx", "Brown hare", "Mountain hare", "Roe deer", "White-tailed deer")

VP_all <- VP |> 
  pluck("vals") |> 
  data.frame() |> 
  rownames_to_column("Predictors") |>
  janitor::clean_names() |> 
  pivot_longer(!predictors) |> 
  mutate(
    predictors = recode(predictors, !!!predictor_names),
    name = recode(name, !!!species_names),
    predictors = factor(
      predictors,
      levels = c(
        "Dead Wood Potential", "Forest cover proportion", "Snow depth (cm)",
        "Peatbogs cover proportion", "Distance to streams (m)", "Streams length (m)",
        "Terrain ruggedness", "Distance to residential areas (m)", "Distance to roads (m)",
        "Sampling effort", "Region random effect", "Spatial random effect", "Year random effect"
      )
    )
  )

# Full model plot
ggplot(VP_all, aes(x = name, y = value, fill = predictors)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", angle = 90, hjust = 1),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(
    y = "Variance Explained",
    x = "Species",
    fill = "Predictors",
    title = "Full model variance partitioning - All species"
  ) + 
  coord_flip()
VP_selected <- VP_all |> 
  filter(name %in% selected_species) |> 
  mutate(name = factor(name, levels = selected_species))

ggplot(VP_selected, aes(x = name, y = value, fill = predictors)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16, face = "bold", angle = 90, hjust = 1),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(
    y = "Variance Explained",
    x = "Species",
    fill = "Predictors",
    title = "Variance partitioning full model"
  ) +
  coord_flip()
#######################################################################################################################################

#2. POSTERIOR BETA ESTIMATES PLOT, MANUSCRIPT FIG.3#####################################



# Transform posterior Beta estimates into tidy format and annotate support
beta_val <- postBeta |> 
  imap(~ .x |>  
         data.frame() |> 
         slice(-1) |> 
         magrittr::set_rownames(c(
           "Deadwood Potential", "Forest Cover proportion", "Peatbogs proportion", 
           "Distance to streams (m)", "Stream length (m)", "Terrain ruggedness",
           "Snow depth (cm)", "Distance to roads (m)", "Distance to residencial area (m)", 
           "Log effort")) |>
         magrittr::set_colnames(c(
           "Brown Hare", "Forest Reindeer", "Lynx", "Moose", "Mountain Hare",
           "Red fox", "Roe Deer", "White-tailed Deer", "Wolf")) |> 
         rownames_to_column("Predictor") |>
         pivot_longer(!Predictor, names_to = "Species", values_to = .y)
  ) |> 
  reduce(left_join, by = c("Predictor", "Species")) |>
  mutate(support_final = ifelse(mean > 0, support, supportNeg)) |> 
  select(Predictor, Species, mean, support_final) |> 
  mutate(
    support_final = case_when(
      support_final < 0.8 ~ "Unsupported",
      support_final >= 0.8 & support_final < 0.9 ~ "Weak support",
      support_final >= 0.9 & support_final < 0.95 ~ "Moderate support",
      support_final >= 0.95 ~ "Strong support",
      TRUE ~ NA_character_
    ),
    support_final = factor(
      support_final,
      levels = c("Unsupported", "Weak support", "Moderate support", "Strong support"),
      ordered = TRUE
    )
  )
#Define Color Gradient for Effects
colors <- colorRampPalette(c('#370f71', 'orangered3'))  # From dark purple to orange-red
colorLevels <- 2
cols_to_use <- colors(colorLevels)



# Filter out species not of main interest
fig_x <- beta_val |> 
  filter(!Species %in% c("Wolf", "Red fox", "Forest Reindeer", "Moose")) |> 
  ggplot(aes(x = Predictor, y = Species, size = support_final, fill = factor(sign(mean)))) +
  labs(x = "Predictors", y = "Species", fill = "Effect:", size = "Support:") +
  geom_point(color = 'gray60', shape = 21) +
  scale_fill_manual(values = cols_to_use, labels = c("Negative", "Positive")) +
  scale_size_manual(values = c(1, 4, 7, 10)) +
  scale_x_discrete(expand = c(0.1, 0.1),
                   labels = c(
                     "Deadwood index" = "Deadwood potential index",
                     "Forest cover proportion" = "Forest cover proportion",
                     "Distance to residential areas (m)" = "Distance to residential areas (m)",
                     "Distance to roads (m)" = "Distance to roads (m)",
                     "Snow depth (cm)" = "Snow depth (cm)",
                     "Peatbogs cover proportion" = "Peatbogs cover proportion",
                     "Distance to streams (m)" = "Distance to streams (m)",
                     "Streams length (m)" = "Streams length (m)",
                     "Terrain ruggedness" = "Terrain ruggedness",
                     "Log effort" = "Sampling effort")) +
  scale_y_discrete(expand = c(0.1, 0.1)) +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white"),
    legend.margin = margin(l = 3, unit = 'cm'),
    legend.title = element_text(hjust = 0.5, size = 15),
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(face = "bold", size = 16),
    text = element_text(size = 14)
  )

print(fig_x)

##################################################################################################################

#3. RESIDUAL SPATIAL LYNX-PREY ASSOCIATIONS PLOT MANUSCRIPT FIG. 5

# Convert regional models to coda format for association analysis
mpost <- convertToCodaObject(m1)  # Central region
mpost <- convertToCodaObject(m2)  # Southern region
#Compute Associations – FULL MODEL 
correlation_full <- computeAssociations(m) |> 
  set_names(c("Triangles", "Year", "Region")) |> 
  map(~ .x |> 
        imap(~ .x |>  
               data.frame() |> 
               magrittr::set_rownames(c("Brown Hare", "Forest Reindeer", "Lynx", "Moose", "Mountain Hare", 
                                        "Red fox", "Roe Deer", "White-tailed Deer", "Wolf")) |>
               magrittr::set_colnames(c("Brown Hare", "Forest Reindeer", "Lynx", "Moose", "Mountain Hare", 
                                        "Red fox", "Roe Deer", "White-tailed Deer", "Wolf")) |> 
               rownames_to_column("Pair_1") |>
               pivot_longer(!Pair_1, names_to = "Pair_2", values_to = .y)
        ) |> 
        reduce(left_join, by = c("Pair_1", "Pair_2")) |>
        mutate(
          support = case_when(
            support > 0.2 & support < 0.8 ~ "Unsupported",
            (support >= 0.8 & support < 0.9) | (support > 0.1 & support <= 0.2) ~ "Weak support",
            (support >= 0.9 & support < 0.95) | (support > 0.05 & support <= 0.1) ~ "Moderate support",
            (support >= 0.95) | (support <= 0.05) ~ "Strong support",
            TRUE ~ NA_character_
          ),
          support = factor(support, levels = c("Unsupported", "Weak support", "Moderate support", "Strong support"), ordered = TRUE)
        )
  ) |> 
  bind_rows(.id = "Random_effect") |> 
  filter(Random_effect == "Triangles") |> 
  mutate(model = "Full model")

#Compute Associations – Southern Region

correlation_southwest <- computeAssociations(m2) |> 
  set_names(c("Triangles", "Year")) |> 
  map(~ .x |> 
        imap(~ .x |>  
               data.frame() |> 
               magrittr::set_rownames(c("Brown Hare", "Forest Reindeer", "Lynx", "Moose", "Mountain Hare", 
                                        "Red fox", "Roe Deer", "White-tailed Deer", "Wolf")) |>
               magrittr::set_colnames(c("Brown Hare", "Forest Reindeer", "Lynx", "Moose", "Mountain Hare", 
                                        "Red fox", "Roe Deer", "White-tailed Deer", "Wolf")) |> 
               rownames_to_column("Pair_1") |>
               pivot_longer(!Pair_1, names_to = "Pair_2", values_to = .y)
        ) |> 
        reduce(left_join, by = c("Pair_1", "Pair_2")) |>
        mutate(
          support = case_when(
            support > 0.2 & support < 0.8 ~ "Unsupported",
            (support >= 0.8 & support < 0.9) | (support > 0.1 & support <= 0.2) ~ "Weak support",
            (support >= 0.9 & support < 0.95) | (support > 0.05 & support <= 0.1) ~ "Moderate support",
            (support >= 0.95) | (support <= 0.05) ~ "Strong support",
            TRUE ~ NA_character_
          ),
          support = factor(support, levels = c("Unsupported", "Weak support", "Moderate support", "Strong support"), ordered = TRUE)
        )
  ) |> 
  bind_rows(.id = "Random_effect") |> 
  filter(Random_effect == "Triangles") |> 
  mutate(model = "Southern region")
#Compute Associations – Central Region

correlation_northeast <- computeAssociations(m1) |> 
  set_names(c("Triangles", "Year")) |> 
  map(~ .x |> 
        imap(~ .x |>  
               data.frame() |> 
               magrittr::set_rownames(c("Brown Hare", "Forest Reindeer", "Lynx", "Moose", "Mountain Hare", 
                                        "Red fox", "Roe Deer", "White-tailed Deer", "Wolf")) |>
               magrittr::set_colnames(c("Brown Hare", "Forest Reindeer", "Lynx", "Moose", "Mountain Hare", 
                                        "Red fox", "Roe Deer", "White-tailed Deer", "Wolf")) |> 
               rownames_to_column("Pair_1") |>
               pivot_longer(!Pair_1, names_to = "Pair_2", values_to = .y)
        ) |> 
        reduce(left_join, by = c("Pair_1", "Pair_2")) |>
        mutate(
          support = case_when(
            support > 0.2 & support < 0.8 ~ "Unsupported",
            (support >= 0.8 & support < 0.9) | (support > 0.1 & support <= 0.2) ~ "Weak support",
            (support >= 0.9 & support < 0.95) | (support > 0.05 & support <= 0.1) ~ "Moderate support",
            (support >= 0.95) | (support <= 0.05) ~ "Strong support",
            TRUE ~ NA_character_
          ),
          support = factor(support, levels = c("Unsupported", "Weak support", "Moderate support", "Strong support"), ordered = TRUE)
        )
  ) |> 
  bind_rows(.id = "Random_effect") |> 
  filter(Random_effect == "Triangles") |> 
  mutate(model = "Central region")
#Combine All Correlation Results into One Dataset

correlation_species <- 
  list(correlation_full, correlation_northeast, correlation_southwest) |> 
  bind_rows()
# Plot Species Associations Across Models

correlation_species |> 
  filter(Pair_1 != Pair_2) |> 
  ggplot(aes(x = Pair_1, y = Pair_2, size = support, fill = mean)) +
  geom_point(shape = 21) +
  facet_grid(~model) +
  scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white"),
    legend.margin = margin(l = 1, unit = 'cm'),
    legend.title = element_text(hjust = 0.1),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    axis.title = element_text(face = "bold")
  ) +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_discrete(expand = c(0.1, 0.1)) +
  scale_size_manual(values = c(0, 3, 6, 9)) +
  labs(
    x = "Species",
    y = "Species",
    fill = "Direction"
  )

#######################################################################################

#4. PREDICTIVE VALUES SR2 HISTOGRAM, MANUSCRIPT FIG.4

## Predict species occurrences using the fitted models
preds_central <- computePredictedValues(m1)  # Central Region
preds_south <- computePredictedValues(m2)    # Southern Region

# Evaluate model fit using predictive R² (SR²) values
MF_central <- evaluateModelFit(hM = m1, predY = preds_central)
MF_south <- evaluateModelFit(hM = m2, predY = preds_south)

# Extract and Convert SR² Values

# Ensure SR² values are numeric
MF_central$SR2 <- as.numeric(MF_central$SR2)
MF_south$SR2 <- as.numeric(MF_south$SR2)

# Extract the SR² for Lynx (assuming 2nd species in the dataset is Lynx)
lynxR2_north <- MF_central$SR2[2]
lynxR2_south <- MF_south$SR2[2]



#Histogram of SR² for Lynx in Each Region

par(mfrow = c(1, 2))  # Split the plotting window in 2 columns

# Central Region
hist(lynxR2_north, xlim = c(0, 1), breaks = 10,
     main = paste0("SR² Central Region: ", round(lynxR2_north, 2)),
     xlab = "SR²", col = "lightblue", border = "black")

# Southern Region
hist(lynxR2_south, xlim = c(0, 1), breaks = 10,
     main = paste0("SR² Southern Region: ", round(lynxR2_south, 2)),
     xlab = "SR²", col = "darkgrey", border = "black")

# Reset layout to single plot
par(mfrow = c(1, 1))
# Barplot Comparison of Lynx SR² Between Regions

# Create data for the barplot
regions <- c("Central Region", "Southern Region")
values <- c(lynxR2_north, lynxR2_south)

# Plot comparative barplot
barplot(values, names.arg = regions, col = c("lightblue", "khaki"),
        main = "Explained Variance (SR²) of environmental factors for Lynx",
        ylab = "SR²", ylim = c(0, 0.3),
        cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2)



###############END OF THE CODE #####################

