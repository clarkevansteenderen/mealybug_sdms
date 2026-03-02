library(biomod2)
library(magrittr)
library(tidyverse)
library(terra)
library(ggthemes)
library(earth)

PRESPOINTS = readRDS(file = "data/PRESPOINTS.ALLSPP.rds")
BACKPOINTS = readRDS(file = "data/BACKPOINTS.ALLSPP.rds")

sapply(PRESPOINTS, nrow)

#############################################
# N. viridis
#############################################

# NB!! long first , then lat
# it took me forever to find that this was the problem

n.viridis.coords.pres = PRESPOINTS$n.viridis %>%
  dplyr::select(longitude, latitude)

n.viridis.coords.back = BACKPOINTS$n.viridis %>%
  dplyr::select(longitude, latitude)

n.viridis.coords = rbind(n.viridis.coords.pres, n.viridis.coords.back)
n.viridis.all_resp = c(rep(1, nrow(n.viridis.coords.pres)), rep(0, nrow(n.viridis.coords.back)))

n.viridis.preds = terra::rast("data/n.viridis_current.tif")
#terra::plot(n.viridis.preds)

n.viridis.PRES = PRESPOINTS$n.viridis %>%
  dplyr::select(-c(species, latitude, longitude, source)) %>%
  janitor::clean_names()

n.viridis.BACK = BACKPOINTS$n.viridis %>%
  dplyr::select(-c(latitude, longitude)) %>%
  janitor::clean_names()

# bind the presence and absence data together into one data frame
n.viridis.data = dplyr::bind_rows(n.viridis.PRES, n.viridis.BACK)

rownames(n.viridis.data) = NULL

# Create a vector containing 0 (indicating background points) and 1 (indicating presence points)
n.viridis.vector = c(
  replicate(nrow(n.viridis.PRES), "1"),
  replicate(nrow(n.viridis.BACK), "0")
) 

n.viridis.vector

length(n.viridis.vector)

# Tuned params = rm 3 and fc of Hinge

myBiomodData <- BIOMOD_FormatingData(
  # Ensure this is a simple 0/1 numeric vector
  resp.var = as.numeric(n.viridis.vector),
  # Ensure this is a SpatRaster (from terra)
  expl.var = n.viridis.preds,
  # Force to a basic data.frame
  resp.xy = as.data.frame(n.viridis.coords), 
  resp.name = "N_viridis",
  PA.nb.rep = 0,
  filter.raster = TRUE
)

myBiomodData


avail_opts <- names(bm_ModelingOptions(
  data.type = "binary",
  models = "MAXNET",
  strategy = "default",
  bm.format = myBiomodData
)@options)

maxnet_key <- grep("MAXNET", avail_opts, value = TRUE)

print(maxnet_key)

maxnet_user_options <- list(
  classes = "lqh",    
  regmult = 3   
)

myBiomodOptions <- bm_ModelingOptions(
  data.type = "binary",
  models = c('MAXNET','RF','RFd','GLM',"GBM","ANN","MARS","FDA","CTA"),
  strategy = "user.defined",
  user.val = setNames(
    list(
      list('_allData_allRun' = maxnet_user_options)
    ),
    maxnet_key
  ),
  bm.format = myBiomodData
)

print(myBiomodOptions)


# 1. Generate the standard default options (this avoids the structure error)
myBiomodOptions <- bm_ModelingOptions(
  data.type = 'binary',
  models = c('MAXNET', 'RF', 'RFd', 'GLM', "GBM", "ANN", "MARS", "FDA", "CTA"),
  strategy = 'default',
  bm.format = myBiomodData
)

# specify hinge
myBiomodOptions@options$MAXNET.binary.maxnet.maxnet@args.values[['_allData_allRun']]$classes <- "h"
# specify rm = 3
myBiomodOptions@options$MAXNET.binary.maxnet.maxnet@args.values[['_allData_allRun']]$regmult <- 3.0

# 3. Verify the change
print(myBiomodOptions@options$MAXNET.binary.maxnet.maxnet@args.values[['_allData_allRun']])

myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  modeling.id = "N_viridis",
  models = c('MAXNET', 'RF', 'RFd', 'GLM', "GBM", "ANN", "MARS", "FDA", "CTA"),
  bm.options = myBiomodOptions,
  CV.strategy = 'random',
  CV.nb.rep = 15,
  CV.perc = 0.8,
  metric.eval = c('TSS', 'ROC', 'KAPPA'),
  var.import = 3
)

myBiomodModelOut

# 1. Get raw evaluations
myEval <- get_evaluations(myBiomodModelOut)
head(myEval)

# 2. Calculate Mean and SD using dplyr
eval_summary <- myEval %>%
  group_by(algo, metric.eval) %>%
  summarise(
    Mean = mean(validation, na.rm = TRUE),
    SD = sd(validation, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Round for readability
  mutate(across(c(Mean, SD), ~ round(.x, 3))) %>%
  # Combine Mean and SD into a single string for a clean table: "0.850 (0.04)"
  dplyr::mutate(Display = paste0(Mean, " (+/- ", SD, ")"))

# 3. Pivot to wide format (Algorithms as rows, Metrics as columns)
eval_final <- eval_summary %>%
  select(algo, metric.eval, Display) %>%
  pivot_wider(names_from = metric.eval, values_from = Display) %>%
  arrange(desc(ROC))

print(eval_final)


write_csv(
  eval_final,
  file = paste0("figures/", "TESTER", "_model_metrics.csv")
)

writexl::write_xlsx(eval_final, "figures/n.viridis_model_metrics.xlsx")

biomod_performance = bm_PlotEvalBoxplot(
  bm.out = myBiomodModelOut,
  group.by = c('algo', 'algo'), # Group by Algorithm
  dataset = 'calibration',      # Look at Validation (Testing) scores
  scales = 'free'
)


perf_data <- get_evaluations(myBiomodModelOut)

head(perf_data)

biomod_performance_box <- ggplot(perf_data, aes(x = algo, y = calibration, fill = algo)) +
  geom_boxplot() +
  facet_wrap(~ metric.eval, scales = "free_y") +
  labs(
    x = "Model",
    y = "Calibration",
    subtitle = bquote(italic(.(species)))   # italic species name
  ) +
  scale_fill_tableau() +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold"),  # facet labels bold
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate x labels
  )

biomod_performance_box

biomod_performance_box = biomod_performance$plot
biomod_performance_box = biomod_performance_box +
  labs(x = "Model", y = "Calibration", subtitle = expression(italic("Nipaecoccus viridis"))) +
  scale_fill_tableau()  
biomod_performance_box

ggsave("figures/n.viridis_biomod_performance_box.png", plot = biomod_performance_box, 
       width = 8, height = 5, dpi = 450)

# Get variable importance scores
myVarImp <- get_variables_importance(myBiomodModelOut)

# Plot it
var_importance = bm_PlotVarImpBoxplot(
  bm.out = myBiomodModelOut,
  group.by = c('expl.var', 'algo', 'algo') 
)

myVarImp <- get_variables_importance(myBiomodModelOut)
head(myVarImp)

# Combine importance across algorithms
varimp_df <- myVarImp %>%
  select(algo, expl.var, var.imp)
head(varimp_df)

# Make boxplot
varimp_box <- ggplot(varimp_df, aes(x = expl.var, y = var.imp, fill = algo)) +
  geom_boxplot() +
  facet_wrap(~algo) +
  labs(
    x = "Predictor",
    y = "Variable Importance",
    subtitle = bquote(italic(.("sp")))
  ) +
  scale_fill_tableau() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  );varimp_box

# Aggregate the data to get the average importance per variable per algorithm
var_summary <- var_importance$tab %>%
  group_by(algo, expl.var) %>%
  summarize(mean_imp = mean(var.imp, na.rm = TRUE),
            sd_imp = sd(var.imp, na.rm = TRUE)) %>%
  arrange(algo, desc(mean_imp))

#var_summary
write.csv(var_summary, "figures/n.viridis_var_importance.csv", row.names = F) 

# var_importance.plot = var_importance$plot
# var_importance.plot = var_importance.plot +
#   labs(x = "Model", y = "Contribution", subtitle = expression(italic("Nipaecoccus viridis"))) +
#   scale_fill_tableau() 


var_importance.tab = var_importance$tab

var_importance.boxplot <- ggplot(var_importance.tab, aes(x = expl.var, y = var.imp, fill = expl.var)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  # Create a separate panel for each algorithm
  facet_wrap(~algo, scales = "free_y") + 
  labs(
    x = "Environmental Predictors", 
    y = "Contribution (Variable Importance)", 
    subtitle = expression(italic("D. aberiae")),
    fill = "Predictors"
  ) +
  scale_fill_tableau() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",        # Hide legend since x-axis labels identify the variables
    strip.background = element_rect(fill = "gray95"), # Style the facet headers
    strip.text = element_text(face = "bold")
  )
var_importance.boxplot 

ggsave("figures/n.viridis_var_contrib.png", plot = var_importance.plot, 
       width = 8, height = 8, dpi = 450)

# Plot the response curves for the top variables

predictions <- get_predictions(myBiomodModelOut)
head(predictions)
head(n.viridis.PRES)

env_data <- n.viridis.PRES %>%
  dplyr::mutate(points = row_number())  # add a 'points' column to match predictions

# Merge predictions with environmental values
resp_data <- predictions %>%
  left_join(env_data, by = "points")

head(resp_data)

resp_data_long <- resp_data %>%
  pivot_longer(
    cols = starts_with("bio_"),   # select your environmental variables
    names_to = "expl.var",
    values_to = "expl.val"
  )

resp.plot <- ggplot(resp_data_long, aes(x = expl.val, y = pred, color = algo, fill = algo)) +
  geom_smooth(method = "loess", span = 0.5, se = FALSE) +
  scale_color_tableau() +
  scale_fill_tableau() +
  facet_wrap(~ expl.var, scales = "free_x") +
  theme_classic() +
  labs(
    title = "Model Response Curves",
    subtitle = bquote(italic(.("N.viridis"))),
    y = "Predicted Probability",
    x = "Environmental Variable",
    color = "Model",
    fill = "Model"
  ) +
  theme(legend.position = "right");resp.plot



myMessyPlot <- bm_PlotResponseCurves(
  bm.out = myBiomodModelOut,
  models.chosen = 'all',
  fixed.var = 'median',
  do.plot = FALSE          # <-- Don't draw the plot
)

class(myMessyPlot)


# 2. Extract the underlying data table
resp_data <- myMessyPlot$tab
resp_data$algo <- sapply(strsplit(as.character(resp_data$pred.name), "_"), tail, 1)
head(resp_data)

resp.plot = ggplot(resp_data, aes(x = expl.val, y = pred.val, color = algo, fill = algo)) +
  # geom_smooth automatically calculates the average consensus line and a confidence band
  geom_smooth(method = "loess", span = 0.5, se = FALSE) + 
  scale_color_tableau() +
  scale_fill_tableau() +
  facet_wrap(~ expl.name, scales = "free_x") +
  theme_classic() +
  labs(title = "Model Response Curves",
       subtitle = expression(italic("Nipaecoccus viridis")),
       y = "Probability of Occurrence",
       x = "Environmental Variable", color = "Model", fill = "Model") +
  theme(legend.position = "right")

resp.plot

ggsave("figures/n.viridis_resp_plot.png", plot = resp.plot, 
       width = 8, height = 7, dpi = 450)



myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut,
  models.chosen = 'all',
  em.by = 'all',
  
  # COMBINATION METHOD:
  # We use Weighted Mean to trust the better models (MAXNET) more
  em.algo = c('EMwmean'),
  
  # FILTER -> keep ROC > a threshold
  metric.select = 'ROC',
  metric.select.thresh = 0.85,
  
  # Outputs
  var.import = 3,
  prob.mean = TRUE,
  prob.ci = TRUE
)

# Check if it worked
myBiomodEM

myKeptModels <- get_kept_models(myBiomodEM)

# Print the list
print(myKeptModels)
length(myKeptModels)


####################################################################################
# PROJECT THIS ENSEMBLE MODEL ONTO AFRICA
####################################################################################

africa_ext = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(continent == "Africa")
africa_vect = terra::vect(africa_ext)


years = c("current", "2050", "2070", "2100")
mycols = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
diffcols = rev(colorRampPalette(brewer.pal(9, "PuOr"))(100))

predict_list = list()

for (yr in years) {
  
  print(paste("Processing year:", yr))
  
  # 1. Read raster
  raster_file = paste0("data/n.viridis_", yr, ".tif")
  r = terra::rast(raster_file)
  
  # 2. Crop to Africa extent
  r = terra::crop(r, terra::ext(africa_vect)) # using terra::ext instead of raster::extent
  r = terra::mask(r, africa_vect)             # optional: clips out the ocean
  
  # Set CRS and Names
  crs(r) = "EPSG:4326"
  names(r) = names(n.viridis.preds)
  
  # -------------------------------------------------------------
  # THE BIOMOD2 MAGIC HAPPENS HERE
  # -------------------------------------------------------------
  
  # 3a. Project the individual models onto the new climate raster
  myProj <- BIOMOD_Projection(
    bm.mod = myBiomodModelOut,       # Note: Use the base models here!
    proj.name = paste0("Proj_", yr), # Gives each year a unique folder
    new.env = r,
    models.chosen = 'all',
    build.clamping.mask = FALSE
  )
  
  # 3b. Combine them into the Ensemble forecast
  myEnsProj <- BIOMOD_EnsembleForecasting(
    bm.em = myBiomodEM,              # Now use the Ensemble weights
    bm.proj = myProj,                # Feed it the projections from 3a
    models.chosen = 'all'
  )
  
  # 4. Extract the final spatial raster to save in your list
  # By default, biomod2 spits out prob.mean or prob.wmean based on your EM setup
  final_raster <- get_predictions(myEnsProj)
  
  # Save in list
  predict_list[[yr]] <- final_raster
  
  # Plotting (Optional)
  # plot(final_raster, main = paste("Ensemble Suitability:", yr), col = mycols)
}


current_prob <- predict_list$current / 1000
current_prob_masked <- mask(current_prob, africa_vect)
current_prob_masked <- clamp(current_prob_masked, lower=0, upper=1)
current_prob_masked 

# calculate a threshold for each map
pts <- terra::vect(n.viridis.coords.pres, geom=c("longitude", "latitude"), crs="EPSG:4326")

# extract values specifically from the raster you are plotting
vals_at_presences <- terra::extract(current_prob_masked, pts)
# calculate the NEW threshold from these exact values
THRESH <- quantile(vals_at_presences[,2], probs = 0.1, na.rm = TRUE)
THRESH

map.curr = PLOT.MYMAP(current_prob_masked, "Current", "Nipaecoccus viridis", 
                      MAXP = 1, THRESH = THRESH)
map.curr

map.curr +
  geom_point(data = n.viridis.coords.pres, aes(x = longitude, y = latitude),
             color = "black", fill = "black", shape = 21, size = 1.5)

hist(current_prob_masked, 
     main="Distribution of Suitability Scores", 
     xlab="Probability", 
     col="skyblue", 
     breaks=50)



prob_2050 = predict_list$`2050` / 1000
prob_masked_2050 <- mask(prob_2050, africa_vect)
prob_masked_2050 <- clamp(prob_masked_2050, lower=0, upper=1)
map.2050 = PLOT.MYMAP(prob_masked_2050, "2050", "Nipaecoccus viridis",
                      MAXP = 1, THRESH = THRESH)
map.2050



prob_2070 = predict_list$`2070` / 1000
prob_masked_2070 <- mask(prob_2070, africa_vect)
prob_masked_2070 <- clamp(prob_masked_2070, lower=0, upper=1)
map.2070 = PLOT.MYMAP(prob_masked_2070, "2070", "Nipaecoccus viridis",
                      MAXP = 1, THRESH = THRESH)
map.2070


prob_2100 = predict_list$`2100` / 1000
prob_masked_2100 <- mask(prob_2100, africa_vect)
prob_masked_2100 <- clamp(prob_masked_2100, lower=0, upper=1)
map.2100 = PLOT.MYMAP(prob_masked_2100, "2100", "Nipaecoccus viridis",
                      MAXP = 1, THRESH = THRESH)
map.2100

combomaps = gridExtra::grid.arrange(map.curr, map.2050, map.2070, map.2100)

ggsave("figures/n.viridis_combomaps.svg", plot = combomaps, 
       width = 10, height = 10, dpi = 450)
ggsave("figures/n.viridis_combomaps.png", plot = combomaps, 
       width = 10, height = 10, dpi = 450)

#diff = prob_masked_2100 - current_prob_masked
#PLOT.MYMAP(diff, "Diff", "Nipaecoccus viridis", MAXP = 0.2)


ggplot() +
  geom_spatraster(data = current_prob_masked) +
  scale_fill_gradientn(
    colors = mycols,       # your custom RdYlBu palette
    limits = c(0, 1),      # ensure 0–1 range
    na.value = "white",
    name = "Occurrence Probability"
  ) +
  geom_spatvector(
    data = africa_vect,
    fill = NA,
    color = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    xlim = c(-25, 55),
    ylim = c(-38, 38),
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    title = "Modeled Climatic Suitability for N. viridis",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )




PLOT.MYMAP = function(IN.DATA, YEAR, SPECIES, MAXP = 1, THRESH){
  
  bin2 = THRESH + (MAXP - THRESH) * 0.33
  bin3 = THRESH + (MAXP - THRESH) * 0.66
  
  mymap = ggplot() +
    # 1. Plot the raster (automatically handles terra objects)
    tidyterra::geom_spatraster(data = IN.DATA) +
    
    # 2. Add the country borders over the top
    tidyterra::geom_spatvector(data = africa_vect, fill = NA, color = "black", linewidth = 0.3) +
    
    # 3. Apply your custom colors and handle the ocean (NA)
    # scale_fill_gradientn(
    #   colors = mycols, 
    #   na.value = "white",       # Custom ocean color (or "transparent")
    #   name = "Suitability", 
    #   limits = c(0, MAXP)           # Locks the scale so it's consistent across all years
    # ) +
    scale_fill_stepsn(
      # unsuitable -> low suitability -> moderate suitability -> high suitability
      colors = c("transparent", "grey80", "#fdae61", "#d73027"), # Grey to Red
      breaks = c(0, THRESH, bin2, bin3), 
      limits = c(0, MAXP),
      na.value = "white",
      name = "Suitability",
      #labels = c("< 10p Threshold", "Low", "Medium", "High"),
      labels = c(
        paste0("10p (<", round(THRESH, 2), ")"),
        paste0("Low (", round(bin2, 2), ")"),
        paste0("Medium (", round(bin3, 2), ")"),
        paste0("High (", round(bin3, 2), "-", MAXP, ")")
      ),
      # This ensures the labels are centered in the color blocks
      guide = guide_colorsteps(
        barheight = unit(6, "cm"),
        ticks = FALSE,
        even.steps = TRUE, # This makes the 4 boxes equal height!
        show.limits = FALSE
      )
    ) +
    coord_sf(xlim = c(-25, 55),
             ylim = c(-38, 38),
             expand = FALSE) +
    # 4. Clean up the theme
    theme_classic() +
    labs(
      title = paste0(YEAR, " Climatic Suitability for ", SPECIES),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme(
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA) 
    )
  return(mymap)
}




# change over time

mean_curr <- terra::global(current_prob_masked, fun="mean", na.rm=TRUE)$mean
mean_2050 <- terra::global(prob_masked_2050, fun="mean", na.rm=TRUE)$mean
mean_2070 <- terra::global(prob_masked_2070, fun="mean", na.rm=TRUE)$mean
mean_2100 <- terra::global(prob_masked_2100, fun="mean", na.rm=TRUE)$mean

# Put it in a small table
trend_data <- data.frame(
  Year = c("Current", 2050, 2070, 2100),
  Mean_Suitability = c(mean_curr, mean_2050, mean_2070, mean_2100)
)

print(trend_data)

df_curr <- data.frame(Year = "Current", Value = values(current_prob_masked, na.rm=TRUE))
df_2100 <- data.frame(Year = "2100", Value = values(prob_masked_2100, na.rm=TRUE))
plot_df <- rbind(df_curr, df_2100)

colnames(plot_df) = c("Year", "Prob")
head(plot_df)

density.plot = ggplot(plot_df, aes(x = Prob, fill = Year)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(title = "Shift in Climatic Suitability: Current vs 2100",
       subtitle = expression(italic("Nipaecoccus viridis")),
       x = "Suitability Probability (0-1)",
       y = "Density") +
  scale_fill_manual(values = c("Current" = "yellow", "2100" = "firebrick")) +
  scale_y_continuous(limits = c(0, 10),
                     breaks = seq(0, 10, 2))
density.plot 

ggsave("figures/n.viridis_density_plot.png", plot = density.plot, 
       width = 8, height = 5, dpi = 450)

# A function to get percentage of land in each suitability bracket
get_suitability_breakdown <- function(raster_layer, name) {
  # Total non-NA pixels (landmass)
  total_land_pixels <- terra::freq(raster_layer > -1)$count 
  
  # Count pixels in brackets
  low    <- terra::freq(raster_layer < 0.2)$count
  med    <- terra::freq(raster_layer >= 0.2 & raster_layer < 0.4)$count
  high   <- terra::freq(raster_layer >= 0.4 & raster_layer < 0.6)$count
  v_high <- terra::freq(raster_layer >= 0.6)$count
  
  # Convert to percentage
  data.frame(
    Period = name,
    Low_0.2 = (low / total_land_pixels) * 100,
    Med_0.4 = (med / total_land_pixels) * 100,
    High_0.6 = (high / total_land_pixels) * 100,
    VHigh_Above = (v_high / total_land_pixels) * 100
  )
}

# Run for Current and 2100
breakdown_curr <- get_suitability_breakdown(current_prob_masked, "Current")
breakdown_2100 <- get_suitability_breakdown(prob_masked_2100, "2100")

# Combine and print
final_table <- rbind(breakdown_curr, breakdown_2100)
print(final_table)

stats_summary <- data.frame(
  Metric = c("Mean Suitability", "Max Suitability", "Standard Deviation"),
  Current = c(
    terra::global(current_prob_masked, mean, na.rm=TRUE)$mean,
    terra::global(current_prob_masked, max, na.rm=TRUE)$max,
    terra::global(current_prob_masked, sd, na.rm=TRUE)$sd
  ),
  Future_2100 = c(
    terra::global(prob_masked_2100, mean, na.rm=TRUE)$mean,
    terra::global(prob_masked_2100, max, na.rm=TRUE)$max,
    terra::global(prob_masked_2100, sd, na.rm=TRUE)$sd
  )
)

print(stats_summary)


# Calculate "Suitability Volume"
# expanse() returns area, then we multiply by the cell values
vol_curr <- sum(values(current_prob_masked * terra::cellSize(current_prob_masked, unit="km")), na.rm=TRUE)
vol_2100 <- sum(values(prob_masked_2100 * terra::cellSize(prob_masked_2100, unit="km")), na.rm=TRUE)

cat("Total Habitat Capacity (Current):", round(vol_curr, 0), "Units\n")
cat("Total Habitat Capacity (2100):", round(vol_2100, 0), "Units\n")
cat("Percentage Loss of Capacity:", round((1 - (vol_2100 / vol_curr)) * 100, 2), "%\n")

write.csv(stats_summary, "figures/n.viridis_area_change.csv", row.names = F)


# another way

# current_prob.nvir = terra::rast("N.viridis/proj_Proj_current/proj_Proj_current_N.viridis_ensemble.tif")
# current_prob.nvir <- current_prob.nvir / 1000
# current_prob.nvir_masked <- mask(current_prob.nvir, africa_vect)
# 
# prob.nvir.2100 = terra::rast("N.viridis/proj_Proj_2100/proj_Proj_2100_N.viridis_ensemble.tif")
# prob.nvir.2100 <- prob.nvir.2100 / 1000
# prob.nvir.2100_masked <- mask(prob.nvir.2100, africa_vect)

# 1. Create binary maps (1 = Suitable, 0 = Unsuitable)
# use THRESH here as the 10p threshold
curr_binary <- current_prob_masked >= THRESH
futu_binary <- prob_masked_2100 >= THRESH

plot(curr_binary)
plot(futu_binary)

# 2. Calculate the area of only the "1" cells (Suitable habitat)
# expanse() with transform=TRUE handles the earth's curvature correctly
area_curr <- terra::expanse(curr_binary, unit="km", transform=TRUE, byValue=TRUE)
area_2100 <- terra::expanse(futu_binary, unit="km", transform=TRUE, byValue=TRUE)

# 3. Extract the area for the '1' (Suitable) class
# (expanse returns a table; we want the value where value == 1)
km2_curr <- area_curr$area[2]
km2_2100 <- area_2100$area[2]

cat("Suitable Area (Current):", round(km2_curr, 0), "km2\n")
cat("Suitable Area (2100):", round(km2_2100, 0), "km2\n")
cat("Percentage Change in Range:", round(((km2_2100 - km2_curr) / km2_curr) * 100, 3), "%\n")
