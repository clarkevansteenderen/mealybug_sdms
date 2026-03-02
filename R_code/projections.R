# Project MaxEnt models

options(java.parameters = "-Xmx8g")

if(!dir.exists("raster_projections")){
  dir.create("raster_projections")
}

PRESPOINTS = readRDS(file = "data/PRESPOINTS.ALLSPP.rds")
BACKPOINTS = readRDS(file = "data/BACKPOINTS.ALLSPP.rds")

# get the raster for the world
world_ext = rnaturalearth::ne_countries(
  scale = "medium", 
  returnclass = "sf"
) 

# Convert to SpatVector (terra format)
world_vect = vect(world_ext)

sapply(PRESPOINTS, nrow)

#############################################
# D. aberiae
#############################################

# read in our predictors again
d.aberiae.preds = terra::rast("data/d.aberiae_current.tif")
names(d.aberiae.preds)
# 
# # we removed bio6 and bio8 based on response plots, so remove here
d.aberiae.preds = d.aberiae.preds[[-c(2, 3)]]

# # crop rasters
world_env_layers_d.aberiae = raster::crop(d.aberiae.preds, raster::extent(world_ext))
#terra::plot(world_env_layers_d.aberiae)

# # set the CRS
crs(world_env_layers_d.aberiae) = "EPSG:4326"
names(world_env_layers_d.aberiae) = names(d.aberiae.data)
names(world_env_layers_d.aberiae)

# # use the model to predict climate suitability 
predict_maxent_d.aberiae = terra::predict(model = d.aberiae.model,
                                          object = world_env_layers_d.aberiae,
                                          na.rm = TRUE)

terra::plot(predict_maxent_d.aberiae)

#############################################
# N. viridis
#############################################

# Define years and file suffixes
#years = c("current", "2100")
years = c("current", "2050", "2070", "2100")
mycols = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
diffcols = rev(colorRampPalette(brewer.pal(9, "PuOr"))(100))

# Set extent excluding Antarctica
extent_no_antarctica = ext(-180, 180, -60, 90)

# Store predicted rasters in a list
predict_list = list()

for (yr in years) {
  
  # Read raster
  raster_file = paste0("data/n.viridis_", yr, ".tif")
  r = terra::rast(raster_file)
  
  # Remove layer 3 if needed (adjust indexing as in your original code)
  r = r[[-c(3)]]
  
  # Crop to world extent
  r = raster::crop(r, raster::extent(world_ext))
  
  # Set CRS
  crs(r) = "EPSG:4326"
  
  # Rename layers to match model predictors
  names(r) = names(n.viridis.data)
  
  # Predict with MaxEnt model
  pred = terra::predict(model = n.viridis.model, object = r, na.rm = TRUE)
  
  # Crop to exclude Antarctica
  pred = crop(pred, extent_no_antarctica)
  
  # Save in list
  predict_list[[yr]] = pred
  
  # Plot and save PNG
  # png_filename = paste0("figures/raster_projections/n_viridis_projection.", yr, ".png")
  # png(png_filename, width = 2200, height = 1400, res = 350)
  # plot(pred, xlab = "Longitude", ylab = "Latitude", main = "", col = mycols)
  # lines(world_vect, col = "black", lwd = 0.5)
  # dev.off()
}

# get diff between 2100 and current

diff_r = predict_list[["2100"]] - predict_list[["current"]]

png_filename = paste0("figures/raster_projections/n_viridis_diff_", yr, "-current.png")
png(png_filename, width = 2200, height = 1400, res = 350)
plot(diff_r, xlab = "Longitude", ylab = "Latitude", main = "", col = diffcols)
lines(world_vect, col = "black", lwd = 0.5)
dev.off()

# plot(diff_r > 0.1, col = c("white", "forestgreen"))
# lines(world_vect, col = "black", lwd = 0.5)
# 
# plot(diff_r < -0.1, col = c("NA", "darkred"))
# lines(world_vect, col = "black", lwd = 0.5)
# 
# plot(diff_r == 0)

# 10th-percentile threshold (suitability)
gps.pts = PRESPOINTS$n.viridis %>%
  dplyr::select(latitude, longitude)

presence_vals = terra::extract(predict_list[["current"]], gps.pts)[,2] 
threshold_10 = quantile(presence_vals, probs = 0.10, na.rm = TRUE)
threshold_10

# youden threshold
gps.back = BACKPOINTS$n.viridis %>%
  dplyr::select(latitude, longitude)

back_vals = terra::extract(predict_list[["current"]], gps.back)[,2] 

vals = c(presence_vals, back_vals)
labs = c(rep(1, length(presence_vals)), rep(0, length(back_vals)))  # 1 = presence, 0 = background

# compute ROC and best threshold (Youden)
roc_obj = pROC::roc(labs, vals, quiet = TRUE)
best    = pROC::coords(roc_obj, "best", best.method = "youden", transpose = FALSE)
threshold_youden = as.numeric(best["threshold"])
threshold_youden
# you can also get sensitivity/specificity:
best[c("sensitivity","specificity")]

# binary_map = predict_list[["2100"]] >= threshold_10
# plot(binary_map)
# lines(world_vect, col = "black", lwd = 0.5)
# 
# 
# suitable = predict_list[["current"]] > 0.5
# cell_areas = terra::cellSize(suitable, unit = "m")
# # Mask to only suitable cells
# suitable_area_cells = mask(cell_areas, suitable, maskvalues = 0, updatevalue = NA)
# 
# # Sum the cell areas
# total_area_m2 = global(suitable_area_cells, "sum", na.rm = TRUE)
# 
# # Convert to km²
# total_area_km2 = total_area_m2 / 1e6
# total_area_km2

# current: area suitable > 0.5 = 368140108 km2
#2100: 363356530 km2

##############################################################################
# some loops

# Compute cell areas once (they have the same extent/resolution)
cell_areas = terra::cellSize(predict_list[["current"]], unit = "m")

# Define thresholds
thresholds = seq(0, 1, by = 0.1)

# Initialize an empty list to store results for both periods
area_list = list()

# Loop through each period
for (period in c("current", "2100")) {
  
  pred = predict_list[[period]]
  results = data.frame(
    threshold = thresholds,
    period = period,
    area_km2 = NA_real_
  )
  
  # Loop through thresholds
  for (i in seq_along(thresholds)) {
    thr = thresholds[i]
    suitable = pred > thr
    
    # Mask and sum suitable area
    suitable_area_cells = mask(cell_areas, suitable, maskvalues = 0, updatevalue = NA)
    total_area_m2 = global(suitable_area_cells, "sum", na.rm = TRUE)
    
    # Store converted area (m² → km²)
    results$area_km2[i] = as.numeric(total_area_m2) / 1e6
  }
  
  area_list[[period]] = results
}

# Combine both into one data frame
area_table = do.call(rbind, area_list)

area_wide = tidyr::pivot_wider(
  area_table,
  names_from = period,
  values_from = area_km2
)

area_wide$diff = area_wide$current - area_wide$'2100'

head(area_wide)

ggplot(area_wide, aes(x = threshold, y = diff)) +
  geom_line(color = "firebrick", linewidth = 1) +
  theme_classic() +
  labs(x = "Threshold", y = "Difference in suitable area (km²)")

# Convert wide data to long format for plotting
area_long <- area_wide %>%
  pivot_longer(cols = c(current, `2100`), names_to = "period", values_to = "area_km2")

# Plot both periods
n.viridis.area = ggplot(area_long, aes(x = threshold, y = area_km2 / 1e8, 
                      color = period, shape = period)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("black", "royalblue"),
    labels = c("2100", "Current")
  ) +
  scale_shape_manual(
    values = c(16, 15),
    labels = c("2100", "Current")
  ) +
  scale_y_continuous(
    breaks = seq(3.3, 5.0, by = 0.2),
    limits = c(3.3, 5.0)
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    x = "MaxEnt suitability",
    y = "Suitable area (×10⁸ km²)",
    title = expression("a) " * italic("N. viridis")),
    color = "Period",  # same name for both color and shape
    shape = "Period"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right"
    #legend.title = element_text(face = "bold")
  ) +
  geom_vline(xintercept = threshold_10, linetype = "dashed", color = "grey25") +
  geom_vline(xintercept = threshold_youden, linetype = "dotted", color = "grey25") +
  annotate("text", x = threshold_10 + 0.02, y = Inf, label = "10-percentile", vjust = 2, color = "grey25") +
  annotate("text", x = threshold_youden + 0.02, y = Inf, label = "Youden/TSS", vjust = 3.5, color = "grey25")

n.viridis.area

thresholds = c("10pct" = 0.1, "Youden" = 0.72)

# Pre-compute cell areas (m²)
cell_areas <- terra::cellSize(predict_list[["current"]], unit = "m")

area_at_threshold <- function(pred_raster, thr, cell_areas_m2) {
  suitable <- pred_raster > thr
  suitable_cells <- mask(cell_areas_m2, suitable, maskvalues = 0, updatevalue = NA)
  total_m2 <- terra::global(suitable_cells, "sum", na.rm = TRUE)
  as.numeric(total_m2) / 1e6  # convert to km²
}

results <- data.frame(
  threshold_name = names(thresholds),
  threshold_val  = thresholds,
  area_current   = NA_real_,
  area_2100      = NA_real_
)

for (i in seq_along(thresholds)) {
  thr <- thresholds[i]
  results$area_current[i] <- area_at_threshold(predict_list[["current"]], thr, cell_areas)
  results$area_2100[i]    <- area_at_threshold(predict_list[["2100"]], thr, cell_areas)
}

# Compute difference as 2100 minus current
results$diff_km2   <- results$area_2100 - results$area_current
results$pct_change <- 100 * results$diff_km2 / results$area_current

results.n.viridis = results
results.n.viridis

# key info
n.viridis.area
results.n.viridis

#############################################
# P. burnerae
#############################################

# read in our predictors again
# p.burnerae.preds = terra::rast("data/p.burnerae_2100.tif")
# p.burnerae.preds = p.burnerae.preds[[-c(1, 2, 3, 4)]]
# names(p.burnerae.preds)
# 
# # crop rasters
# world_env_layers_p.burnerae = raster::crop(p.burnerae.preds, raster::extent(world_ext))
# #terra::plot(world_env_layers_p.burnerae)
# 
# # set the CRS
# crs(world_env_layers_p.burnerae) = "EPSG:4326"
# names(world_env_layers_p.burnerae) = names(p.burnerae.data)
# names(world_env_layers_p.burnerae)
# 
# # use the model to predict climate suitability 
# predict_maxent_p.burnerae = terra::predict(model = p.burnerae.model, 
#                                           object = world_env_layers_p.burnerae,
#                                           na.rm = TRUE)
# 
# terra::plot(predict_maxent_p.burnerae)


#############################################
# P. citri
#############################################

# Define years and file suffixes
#years = c("current", "2050", "2070", "2100")
years = c("current", "2100")
mycols = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
diffcols = rev(colorRampPalette(brewer.pal(9, "PuOr"))(100))

# Set extent excluding Antarctica
extent_no_antarctica = ext(-180, 180, -60, 90)

# Store predicted rasters in a list
predict_list = list()

options(java.parameters = "-Xmx16g")

for (yr in years) {
  
  # Read raster
  raster_file = paste0("data/p.citri_", yr, ".tif")
  r = terra::rast(raster_file)
  
  # Remove layer 3 if needed (adjust indexing as in your original code)
  r = r[[-c(3,5)]]
  
  # Crop to world extent
  r = raster::crop(r, raster::extent(world_ext))
  
  # Set CRS
  crs(r) = "EPSG:4326"
  
  # Rename layers to match model predictors
  names(r) = names(p.citri.data)
  
  # Predict with MaxEnt model
  pred = terra::predict(model = p.citri.model, object = r, na.rm = TRUE)
  
  # Crop to exclude Antarctica
  pred = crop(pred, extent_no_antarctica)
  
  # Save in list
  predict_list[[yr]] = pred
  
  # Plot and save PNG
  # png_filename = paste0("figures/raster_projections/p.citri_projection.", yr, ".png")
  # png(png_filename, width = 2200, height = 1400, res = 350)
  # plot(pred, xlab = "Longitude", ylab = "Latitude", main = "", col = mycols)
  # lines(world_vect, col = "black", lwd = 0.5)
  # dev.off()
}

# get diff between 2100 and current

# diff_r = predict_list[["2100"]] - predict_list[["current"]]
# 
# png_filename = paste0("figures/raster_projections/p.citri_diff_", yr, "-current.png")
# png(png_filename, width = 2200, height = 1400, res = 350)
# plot(diff_r, xlab = "Longitude", ylab = "Latitude", main = "", col = diffcols)
# lines(world_vect, col = "black", lwd = 0.5)
# dev.off()

# 10th-percentile threshold (suitability)
gps.pts = PRESPOINTS$p.citri %>%
  dplyr::select(latitude, longitude)

presence_vals = terra::extract(predict_list[["current"]], gps.pts)[,2] 
threshold_10 = quantile(presence_vals, probs = 0.10, na.rm = TRUE)
threshold_10

# youden threshold
gps.back = BACKPOINTS$p.citri %>%
  dplyr::select(latitude, longitude)

back_vals = terra::extract(predict_list[["current"]], gps.back)[,2] 

vals = c(presence_vals, back_vals)
labs = c(rep(1, length(presence_vals)), rep(0, length(back_vals)))  # 1 = presence, 0 = background

# compute ROC and best threshold (Youden)
roc_obj = pROC::roc(labs, vals, quiet = TRUE)
best    = pROC::coords(roc_obj, "best", best.method = "youden", transpose = FALSE)
threshold_youden = as.numeric(best["threshold"])
threshold_youden
# you can also get sensitivity/specificity:
best[c("sensitivity","specificity")]

# binary_map = predict_list[["2100"]] >= threshold_10
# plot(binary_map)
# lines(world_vect, col = "black", lwd = 0.5)
# 
# 
# suitable = predict_list[["current"]] > 0.5
# cell_areas = terra::cellSize(suitable, unit = "m")
# # Mask to only suitable cells
# suitable_area_cells = mask(cell_areas, suitable, maskvalues = 0, updatevalue = NA)
# 
# # Sum the cell areas
# total_area_m2 = global(suitable_area_cells, "sum", na.rm = TRUE)
# 
# # Convert to km²
# total_area_km2 = total_area_m2 / 1e6
# total_area_km2

# current: area suitable > 0.5 = 368140108 km2
#2100: 363356530 km2

##############################################################################
# some loops

# Compute cell areas once (they have the same extent/resolution)
cell_areas = terra::cellSize(predict_list[["current"]], unit = "m")

# Define thresholds
thresholds = seq(0, 1, by = 0.1)

# Initialize an empty list to store results for both periods
area_list = list()

# Loop through each period
for (period in c("current", "2100")) {
  
  pred = predict_list[[period]]
  results = data.frame(
    threshold = thresholds,
    period = period,
    area_km2 = NA_real_
  )
  
  # Loop through thresholds
  for (i in seq_along(thresholds)) {
    thr = thresholds[i]
    suitable = pred > thr
    
    # Mask and sum suitable area
    suitable_area_cells = mask(cell_areas, suitable, maskvalues = 0, updatevalue = NA)
    total_area_m2 = global(suitable_area_cells, "sum", na.rm = TRUE)
    
    # Store converted area (m² → km²)
    results$area_km2[i] = as.numeric(total_area_m2) / 1e6
  }
  
  area_list[[period]] = results
}

# Combine both into one data frame
area_table = do.call(rbind, area_list)

area_wide = tidyr::pivot_wider(
  area_table,
  names_from = period,
  values_from = area_km2
)

area_wide$diff = area_wide$current - area_wide$'2100'

head(area_wide)

ggplot(area_wide, aes(x = threshold, y = diff)) +
  geom_line(color = "firebrick", linewidth = 1) +
  theme_classic() +
  labs(x = "Threshold", y = "Difference in suitable area (km²)")

# Convert wide data to long format for plotting
area_long <- area_wide %>%
  pivot_longer(cols = c(current, `2100`), names_to = "period", values_to = "area_km2")

# Plot both periods
p.citri.area = ggplot(area_long, aes(x = threshold, y = area_km2 / 1e8, 
                                       color = period, shape = period)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("black", "royalblue"),
    labels = c("2100", "Current")
  ) +
  scale_shape_manual(
    values = c(16, 15),
    labels = c("2100", "Current")
  ) +
  scale_y_continuous(
    breaks = seq(3.3, 5.0, by = 0.2),
    limits = c(3.3, 5.0)
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    x = "MaxEnt suitability",
    y = "Suitable area (×10⁸ km²)",
    title = expression("b) " * italic("P. citri")),
    color = "Period",  # same name for both color and shape
    shape = "Period"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right"
    #legend.title = element_text(face = "bold")
  ) +
  geom_vline(xintercept = threshold_10, linetype = "dashed", color = "grey25") +
  geom_vline(xintercept = threshold_youden, linetype = "dotted", color = "grey25") +
  annotate("text", x = threshold_10 + 0.02, y = Inf, label = "10-percentile", vjust = 2, color = "grey25") +
  annotate("text", x = threshold_youden + 0.15, y = Inf, label = "Youden/TSS", vjust = 3.5, color = "grey25")

p.citri.area

thresholds = c("10pct" = threshold_10, "Youden" = threshold_youden)

# Pre-compute cell areas (m²)
cell_areas <- terra::cellSize(predict_list[["current"]], unit = "m")

area_at_threshold <- function(pred_raster, thr, cell_areas_m2) {
  suitable <- pred_raster > thr
  suitable_cells <- mask(cell_areas_m2, suitable, maskvalues = 0, updatevalue = NA)
  total_m2 <- terra::global(suitable_cells, "sum", na.rm = TRUE)
  as.numeric(total_m2) / 1e6  # convert to km²
}

results <- data.frame(
  threshold_name = names(thresholds),
  threshold_val  = thresholds,
  area_current   = NA_real_,
  area_2100      = NA_real_
)

for (i in seq_along(thresholds)) {
  thr <- thresholds[i]
  results$area_current[i] <- area_at_threshold(predict_list[["current"]], thr, cell_areas)
  results$area_2100[i]    <- area_at_threshold(predict_list[["2100"]], thr, cell_areas)
}

# Compute difference as 2100 minus current
results$diff_km2   <- results$area_2100 - results$area_current
results$pct_change <- 100 * results$diff_km2 / results$area_current

results.p.citri = results
results.p.citri

# key info
p.citri.area
results.p.citri


p.citri.area
n.viridis.area

ggsave("figures/n.viridis.area.png", n.viridis.area, width = 8, height = 6, dpi = 350)
ggsave("figures/p.citri.area.png", p.citri.area, width = 8, height = 6, dpi = 350)

results.n.viridis$sp = "n.viridis"
results.p.citri$sp = "p.citri"
results = rbind(results.n.viridis, results.p.citri)
write.csv(results, "data/area.changes.csv")
writexl::write_xlsx(results, "data/area.changes.xlsx")

#############################################
# P. calceolariae
#############################################

options(java.parameters = "-Xmx8g")

# Define years and file suffixes
years = c("current", "2100")
mycols = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
diffcols = rev(colorRampPalette(brewer.pal(9, "PuOr"))(100))

# Set extent excluding Antarctica
extent_no_antarctica = ext(-180, 180, -60, 90)

# Store predicted rasters in a list
predict_list = list()

for (yr in years) {
  
  # Read raster
  raster_file = paste0("data/p.calceolariae_", yr, ".tif")
  r = terra::rast(raster_file)
  
  # Remove layer 3 if needed (adjust indexing as in your original code)
  r = r[[-c(5)]]
  
  # Crop to world extent
  r = raster::crop(r, raster::extent(world_ext))
  
  # Set CRS
  crs(r) = "EPSG:4326"
  
  # Rename layers to match model predictors
  names(r) = names(p.calceolariae.data)
  
  # Predict with MaxEnt model
  pred = terra::predict(model = p.calceolariae.model, object = r, na.rm = TRUE)
  
  # Crop to exclude Antarctica
  pred = crop(pred, extent_no_antarctica)
  
  # Save in list
  predict_list[[yr]] = pred
  
  # Plot and save PNG
  # png_filename = paste0("figures/raster_projections/p.calceolariae_projection.", yr, ".png")
  # png(png_filename, width = 2200, height = 1400, res = 350)
  # plot(pred, xlab = "Longitude", ylab = "Latitude", main = "", col = mycols)
  # lines(world_vect, col = "black", lwd = 0.5)
  # dev.off()
}

# get diff between 2100 and current

diff_r = predict_list[["2100"]] - predict_list[["current"]]

png_filename = paste0("figures/raster_projections/p.calceolariae_diff_", yr, "-current.png")
png(png_filename, width = 2200, height = 1400, res = 350)
plot(diff_r, xlab = "Longitude", ylab = "Latitude", main = "", col = diffcols)
lines(world_vect, col = "black", lwd = 0.5)
dev.off()

##############################################################################
# some loops

# Compute cell areas once (they have the same extent/resolution)
cell_areas = terra::cellSize(predict_list[["current"]], unit = "m")

# Define thresholds
thresholds = seq(0, 1, by = 0.1)

# Initialize an empty list to store results for both periods
area_list = list()

# Loop through each period
for (period in c("current", "2100")) {
  
  pred = predict_list[[period]]
  results = data.frame(
    threshold = thresholds,
    period = period,
    area_km2 = NA_real_
  )
  
  # Loop through thresholds
  for (i in seq_along(thresholds)) {
    thr = thresholds[i]
    suitable = pred > thr
    
    # Mask and sum suitable area
    suitable_area_cells = mask(cell_areas, suitable, maskvalues = 0, updatevalue = NA)
    total_area_m2 = global(suitable_area_cells, "sum", na.rm = TRUE)
    
    # Store converted area (m² → km²)
    results$area_km2[i] = as.numeric(total_area_m2) / 1e6
  }
  
  area_list[[period]] = results
}

# Combine both into one data frame
area_table = do.call(rbind, area_list)

area_wide = tidyr::pivot_wider(
  area_table,
  names_from = period,
  values_from = area_km2
)

area_wide$diff = area_wide$current - area_wide$'2100'

head(area_wide)

ggplot(area_wide, aes(x = threshold, y = diff)) +
  geom_line(color = "firebrick", linewidth = 1) +
  theme_classic() +
  labs(x = "Threshold", y = "Difference in suitable area (km²)")

# Convert wide data to long format for plotting
area_long <- area_wide %>%
  pivot_longer(cols = c(current, `2100`), names_to = "period", values_to = "area_km2")

# 10th-percentile threshold (suitability)
gps.pts = PRESPOINTS$p.calceolariae %>%
  dplyr::select(latitude, longitude)

presence_vals = terra::extract(predict_list[["current"]], gps.pts)[,2] 
threshold_10 = quantile(presence_vals, probs = 0.10, na.rm = TRUE)
threshold_10

# youden threshold
gps.back = BACKPOINTS$p.calceolariae %>%
  dplyr::select(latitude, longitude)

back_vals = terra::extract(predict_list[["current"]], gps.back)[,2] 

vals = c(presence_vals, back_vals)
labs = c(rep(1, length(presence_vals)), rep(0, length(back_vals)))  # 1 = presence, 0 = background

# compute ROC and best threshold (Youden)
roc_obj = pROC::roc(labs, vals, quiet = TRUE)
best    = pROC::coords(roc_obj, "best", best.method = "youden", transpose = FALSE)
threshold_youden = as.numeric(best["threshold"])
threshold_youden
# you can also get sensitivity/specificity:
best[c("sensitivity","specificity")]

# Plot both periods
p.calceolariae.area = ggplot(area_long, aes(x = threshold, y = area_km2 / 1e8, 
                                     color = period, shape = period)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("black", "royalblue"),
    labels = c("2100", "Current")
  ) +
  scale_shape_manual(
    values = c(16, 15),
    labels = c("2100", "Current")
  ) +
  scale_y_continuous(
    breaks = seq(3.3, 5.0, by = 0.2),
    limits = c(3.3, 5.0)
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    x = "MaxEnt suitability",
    y = "Suitable area (×10⁸ km²)",
    title = expression("c) " * italic("P. calceolariae")),
    color = "Period",  # same name for both color and shape
    shape = "Period"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right"
    #legend.title = element_text(face = "bold")
  ) +
  geom_vline(xintercept = threshold_10, linetype = "dashed", color = "grey25") +
  geom_vline(xintercept = threshold_youden, linetype = "dotted", color = "grey25") +
  annotate("text", x = threshold_10 + 0.02, y = Inf, label = "10-percentile", vjust = 2, color = "grey25") +
  annotate("text", x = threshold_youden + 0.15, y = Inf, label = "Youden/TSS", vjust = 3.5, color = "grey25")

p.calceolariae.area

thresholds = c("10pct" = threshold_10, "Youden" = threshold_youden)

# Pre-compute cell areas (m²)
cell_areas <- terra::cellSize(predict_list[["current"]], unit = "m")

area_at_threshold <- function(pred_raster, thr, cell_areas_m2) {
  suitable <- pred_raster > thr
  suitable_cells <- mask(cell_areas_m2, suitable, maskvalues = 0, updatevalue = NA)
  total_m2 <- terra::global(suitable_cells, "sum", na.rm = TRUE)
  as.numeric(total_m2) / 1e6  # convert to km²
}

results <- data.frame(
  threshold_name = names(thresholds),
  threshold_val  = thresholds,
  area_current   = NA_real_,
  area_2100      = NA_real_
)

for (i in seq_along(thresholds)) {
  thr <- thresholds[i]
  results$area_current[i] <- area_at_threshold(predict_list[["current"]], thr, cell_areas)
  results$area_2100[i]    <- area_at_threshold(predict_list[["2100"]], thr, cell_areas)
}

# Compute difference as 2100 minus current
results$diff_km2   <- results$area_2100 - results$area_current
results$pct_change <- 100 * results$diff_km2 / results$area_current

results.p.calceolariae = results
results.p.calceolariae

# key info
p.calceolariae.area
results.p.calceolariae

ggsave("figures/p.calceolariae.area.png", p.calceolariae.area, width = 8, height = 6, dpi = 350)

#############################################
# P. longispinus
#############################################

# Define years and file suffixes
years = c("current", "2100")
mycols = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
diffcols = rev(colorRampPalette(brewer.pal(9, "PuOr"))(100))

# Set extent excluding Antarctica
extent_no_antarctica = ext(-180, 180, -60, 90)

# Store predicted rasters in a list
predict_list = list()

for (yr in years) {
  
  # Read raster
  raster_file = paste0("data/p.longispinus_", yr, ".tif")
  r = terra::rast(raster_file)
  
  # Remove layer 3 if needed (adjust indexing as in your original code)
  r = r[[-c(6)]]
  
  # Crop to world extent
  r = raster::crop(r, raster::extent(world_ext))
  
  # Set CRS
  crs(r) = "EPSG:4326"
  
  # Rename layers to match model predictors
  names(r) = names(p.longispinus.data)
  
  # Predict with MaxEnt model
  pred = terra::predict(model = p.longispinus.model, object = r, na.rm = TRUE)
  
  # Crop to exclude Antarctica
  pred = crop(pred, extent_no_antarctica)
  
  # Save in list
  predict_list[[yr]] = pred
  
  # Plot and save PNG
  # png_filename = paste0("figures/raster_projections/p.longispinus_projection.", yr, ".png")
  # png(png_filename, width = 2200, height = 1400, res = 350)
  # plot(pred, xlab = "Longitude", ylab = "Latitude", main = "", col = mycols)
  # lines(world_vect, col = "black", lwd = 0.5)
  # dev.off()
}

# get diff between 2100 and current

# diff_r = predict_list[["2100"]] - predict_list[["current"]]
# 
# png_filename = paste0("figures/raster_projections/p.longispinus_diff_", yr, "-current.png")
# png(png_filename, width = 2200, height = 1400, res = 350)
# plot(diff_r, xlab = "Longitude", ylab = "Latitude", main = "", col = diffcols)
# lines(world_vect, col = "black", lwd = 0.5)
# dev.off()

##############################################################################
# some loops

# Compute cell areas once (they have the same extent/resolution)
cell_areas = terra::cellSize(predict_list[["current"]], unit = "m")

# Define thresholds
thresholds = seq(0, 1, by = 0.1)

# Initialize an empty list to store results for both periods
area_list = list()

# Loop through each period
for (period in c("current", "2100")) {
  
  pred = predict_list[[period]]
  results = data.frame(
    threshold = thresholds,
    period = period,
    area_km2 = NA_real_
  )
  
  # Loop through thresholds
  for (i in seq_along(thresholds)) {
    thr = thresholds[i]
    suitable = pred > thr
    
    # Mask and sum suitable area
    suitable_area_cells = mask(cell_areas, suitable, maskvalues = 0, updatevalue = NA)
    total_area_m2 = global(suitable_area_cells, "sum", na.rm = TRUE)
    
    # Store converted area (m² → km²)
    results$area_km2[i] = as.numeric(total_area_m2) / 1e6
  }
  
  area_list[[period]] = results
}

# Combine both into one data frame
area_table = do.call(rbind, area_list)

area_wide = tidyr::pivot_wider(
  area_table,
  names_from = period,
  values_from = area_km2
)

area_wide$diff = area_wide$current - area_wide$'2100'

head(area_wide)

ggplot(area_wide, aes(x = threshold, y = diff)) +
  geom_line(color = "firebrick", linewidth = 1) +
  theme_classic() +
  labs(x = "Threshold", y = "Difference in suitable area (km²)")

# Convert wide data to long format for plotting
area_long <- area_wide %>%
  pivot_longer(cols = c(current, `2100`), names_to = "period", values_to = "area_km2")

# 10th-percentile threshold (suitability)
gps.pts = PRESPOINTS$p.longispinus %>%
  dplyr::select(latitude, longitude)

presence_vals = terra::extract(predict_list[["current"]], gps.pts)[,2] 
threshold_10 = quantile(presence_vals, probs = 0.10, na.rm = TRUE)
threshold_10

# youden threshold
gps.back = BACKPOINTS$p.longispinus %>%
  dplyr::select(latitude, longitude)

back_vals = terra::extract(predict_list[["current"]], gps.back)[,2] 

vals = c(presence_vals, back_vals)
labs = c(rep(1, length(presence_vals)), rep(0, length(back_vals)))  # 1 = presence, 0 = background

# compute ROC and best threshold (Youden)
roc_obj = pROC::roc(labs, vals, quiet = TRUE)
best    = pROC::coords(roc_obj, "best", best.method = "youden", transpose = FALSE)
threshold_youden = as.numeric(best["threshold"])
threshold_youden
# you can also get sensitivity/specificity:
best[c("sensitivity","specificity")]

# Plot both periods
p.longispinus.area = ggplot(area_long, aes(x = threshold, y = area_km2 / 1e8, 
                                            color = period, shape = period)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("black", "royalblue"),
    labels = c("2100", "Current")
  ) +
  scale_shape_manual(
    values = c(16, 15),
    labels = c("2100", "Current")
  ) +
  scale_y_continuous(
    breaks = seq(3.3, 5.0, by = 0.2),
    limits = c(3.3, 5.0)
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    x = "MaxEnt suitability",
    y = "Suitable area (×10⁸ km²)",
    title = expression("d) " * italic("P. longispinus")),
    color = "Period",  # same name for both color and shape
    shape = "Period"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right"
    #legend.title = element_text(face = "bold")
  ) +
  geom_vline(xintercept = threshold_10, linetype = "dashed", color = "grey25") +
  geom_vline(xintercept = threshold_youden, linetype = "dotted", color = "grey25") +
  annotate("text", x = threshold_10 + 0.05, y = Inf, label = "10-percentile", vjust = 2, color = "grey25") +
  annotate("text", x = threshold_youden + 0.02, y = Inf, label = "Youden/TSS", vjust = 3.5, color = "grey25")

p.longispinus.area

thresholds = c("10pct" = threshold_10, "Youden" = threshold_youden)

# Pre-compute cell areas (m²)
cell_areas <- terra::cellSize(predict_list[["current"]], unit = "m")

area_at_threshold <- function(pred_raster, thr, cell_areas_m2) {
  suitable <- pred_raster > thr
  suitable_cells <- mask(cell_areas_m2, suitable, maskvalues = 0, updatevalue = NA)
  total_m2 <- terra::global(suitable_cells, "sum", na.rm = TRUE)
  as.numeric(total_m2) / 1e6  # convert to km²
}

results <- data.frame(
  threshold_name = names(thresholds),
  threshold_val  = thresholds,
  area_current   = NA_real_,
  area_2100      = NA_real_
)

for (i in seq_along(thresholds)) {
  thr <- thresholds[i]
  results$area_current[i] <- area_at_threshold(predict_list[["current"]], thr, cell_areas)
  results$area_2100[i]    <- area_at_threshold(predict_list[["2100"]], thr, cell_areas)
}

# Compute difference as 2100 minus current
results$diff_km2   <- results$area_2100 - results$area_current
results$pct_change <- 100 * results$diff_km2 / results$area_current

results.p.longispinus = results
results.p.longispinus

# key info
p.longispinus.area
results.p.longispinus

ggsave("figures/p.longispinus.area.png", p.longispinus.area, width = 8, height = 6, dpi = 350)


results.p.calceolariae$sp = "p.calceolariae"
results.p.longispinus$sp = "p.longispinus"

results = rbind(results.p.calceolariae, results.p.longispinus)

writexl::write_xlsx(results, "data/results.xlsx")
