# downstream plotting after you have the raster outputs from the ensemble modelling

######################################################## 
# MAP PLOTTING FUNCTIONS
######################################################## 

# Continuous map
PLOT.CONT.MAP = function(IN.DATA, YEAR){
  
  mycols = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
  
  map = ggplot() +
    tidyterra::geom_spatraster(data = IN.DATA) +
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
      title = YEAR,
      x = "Longitude",
      y = "Latitude"
    ) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right"
    )
  
  return(map)
}

# binary thresholded map

PLOT.MYMAP = function(IN.DATA, YEAR, MAXP = 1, THRESH){
  
  # Define bins
  bin2 = THRESH + (MAXP - THRESH) * 0.33
  bin3 = THRESH + (MAXP - THRESH) * 0.66
  
  # Reclassify into 4 discrete categories
  IN.CLASS <- terra::classify(
    IN.DATA,
    rbind(
      c(0, THRESH, 1),
      c(THRESH, bin2, 2),
      c(bin2, bin3, 3),
      c(bin3, MAXP, 4)
    ),
    include.lowest = TRUE
  )
  
  IN.CLASS <- terra::as.factor(IN.CLASS)
  
  mymap = ggplot() +
    
    # Plot categorical raster
    tidyterra::geom_spatraster(data = IN.CLASS) +
    
    # Country borders
    tidyterra::geom_spatvector(
      data = africa_vect,
      fill = NA,
      color = "black",
      linewidth = 0.3
    ) +
    
    # Manual discrete colours
    scale_fill_manual(
      values = c(
        "1" = "grey80",
        "2" = "darkmagenta",
        "3" = "orange",
        "4" = "#d73027"
      ),
      name = "Suitability",
      labels = c(
        paste0("< ", round(THRESH, 2)),
        paste0("Low (", round(THRESH, 2), "–", round(bin2, 2), ")"),
        paste0("Moderate (", round(bin2, 2), "–", round(bin3, 2), ")"),
        paste0("High (", round(bin3, 2), "–", round(MAXP, 2), ")")
      ),
      na.value = "white",
      na.translate = FALSE   # ← this removes NA from legend
    ) +
    
    coord_sf(
      xlim = c(-25, 55),
      ylim = c(-38, 38),
      expand = FALSE
    ) +
    
    theme_classic() +
    
    labs(
      title = paste0(YEAR),
      x = "Longitude",
      y = "Latitude"
    ) +
    
    theme(
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(mymap)
}

########################################################
########################################################

africa_ext = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(continent == "Africa")
africa_vect = terra::vect(africa_ext)

PRESPOINTS = readRDS(file = "data/PRESPOINTS.ALLSPP.rds")
BACKPOINTS = readRDS(file = "data/BACKPOINTS.ALLSPP.rds")

species_key <- c(
  "Delottococcus aberiae"     = "d.aberiae",
  "Nippaecoccus viridis"      = "n.viridis",
  "Planococcus citri"         = "p.citri",
  "Pseudococcus calceolariae" = "p.calceolariae",
  "Pseudococcus longispinus"  = "p.longispinus"
)

# set this each time
species <- "Pseudococcus longispinus"
species_code <- species_key[species]


coords.pres = PRESPOINTS[[species_code]] %>%
  dplyr::select(longitude, latitude)

coords.back = BACKPOINTS[[species_code]] %>%
  dplyr::select(longitude, latitude)

current_prob = terra::rast(paste0("projections/", species_code, "_ensemble_current.tif"))
current_prob <- current_prob / 1000
current_prob.masked <- mask(current_prob, africa_vect)
current_prob.masked <- clamp(current_prob.masked, lower=0, upper=1)

prob_2050 = terra::rast(paste0("projections/", species_code, "_ensemble_2050.tif"))
prob_2050 <- prob_2050 / 1000
prob_2050.masked <- mask(prob_2050, africa_vect)
prob_2050.masked <- clamp(prob_2050.masked, lower=0, upper=1)

prob_2070 = terra::rast(paste0("projections/", species_code, "_ensemble_2070.tif"))
prob_2070 <- prob_2070 / 1000
prob_2070.masked <- mask(prob_2070, africa_vect)
prob_2070.masked <- clamp(prob_2070.masked, lower=0, upper=1)

prob_2100 = terra::rast(paste0("projections/", species_code, "_ensemble_2100.tif"))
prob_2100 <- prob_2100 / 1000
prob_2100.masked <- mask(prob_2100, africa_vect)
prob_2100.masked <- clamp(prob_2100.masked, lower=0, upper=1)


##########################################
# 10p threshold

pres_vect <- vect(
  coords.pres,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)
pres_vals <- extract(current_prob.masked, pres_vect)
# Extract returns a data frame with ID + value column
head(pres_vals)
# Remove ID column
pres_suitability <- pres_vals[,2]
ten.p_thresh <- quantile(pres_suitability, 0.10, na.rm = TRUE)

print(paste("The 10p threshold is:", ten.p_thresh))


##########################################

# Get the maxSS
# optimises the balance between true + and true -
abs_vals <- extract(current_prob.masked, coords.back)
abs_suitability <- abs_vals[,2]

# Use the evaluate function
eval <- dismo::evaluate(p = pres_suitability, a = abs_suitability)

all_thresholds <- dismo::threshold(eval)

# Access MaxSS directly
max_ss_val <- all_thresholds$spec_sens
all_thresholds$sensitivity

print(paste("The MaxSS threshold is:", max_ss_val))

##########################################

THRESH = ten.p_thresh

# binary thresholded maps

map.curr = PLOT.MYMAP(current_prob.masked, "Current", species, 
                      MAXP = 1, THRESH = THRESH)

# add presence points to check
# map.curr +
#   geom_point(data = coords.pres, aes(x = longitude, y = latitude),
#              color = "black", fill = "black", shape = 21, size = 1.5)
# 
# hist(current_prob.masked, 
#      main="Distribution of Suitability Scores", 
#      xlab="Probability", 
#      col="skyblue", 
#      breaks=50)

map.2050 = PLOT.MYMAP(prob_2050.masked, "2050", species,
                      MAXP = 1, THRESH = THRESH)

map.2070 = PLOT.MYMAP(prob_2070.masked, "2070", species,
                      MAXP = 1, THRESH = THRESH)

map.2100 = PLOT.MYMAP(prob_2100.masked, "2100", species,
                      MAXP = 1, THRESH = THRESH)

# Continuous maps

map.curr.cont = PLOT.CONT.MAP(current_prob.masked, "Current")

map.2050.cont = PLOT.CONT.MAP(prob_2050.masked, "2050")

map.2070.cont = PLOT.CONT.MAP(prob_2070.masked, "2070")

map.2100.cont = PLOT.CONT.MAP(prob_2100.masked, "2100")

# binary maps

map_list <- list(map.curr, map.2050, map.2070, map.2100)
map_list <- lapply(map_list, function(x) x + theme(plot.margin = margin(0,0,0,0, "pt")))

# create a version where only the last map has a legend
map_list_noleg <- map_list
map_list_noleg[-length(map_list_noleg)] <- 
  lapply(map_list_noleg[-length(map_list_noleg)], function(p) {
    p + theme(legend.position = "none")
  })

# Create side strip
species_strip <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, 
           label = species,
           angle = 90, size = 5, fontface = "italic") +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "grey90", colour = NA)
  )

combomaps <- species_strip + 
  wrap_plots(map_list_noleg, nrow = 1) +
  plot_layout(widths = c(0.05, 1))

ggsave(paste0("figures/", species, "_combomaps.svg"), plot = combomaps, 
       width = 10, height = 10, dpi = 450)
ggsave(paste0("figures/", species, "_combomaps.png"), plot = combomaps, 
       width = 14, height = 3, dpi = 450)


# continuous maps

map_list_cont <- list(map.curr.cont, map.2050.cont, map.2070.cont, map.2100.cont)
map_list_cont <- lapply(map_list_cont, function(x) x + theme(plot.margin = margin(0,0,0,0, "pt")))

# create a version where only the last map has a legend
map_list_cont_noleg <- map_list_cont
map_list_cont_noleg[-length(map_list_cont_noleg)] <- 
  lapply(map_list_cont_noleg[-length(map_list_cont_noleg)], function(p) {
    p + theme(legend.position = "none")
  })

# Create side strip
species_strip <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, 
           label = species,
           angle = 90, size = 5, fontface = "italic") +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "grey90", colour = NA)
  )

combomaps.cont <- species_strip + 
  wrap_plots(map_list_cont_noleg, nrow = 1) +
  plot_layout(widths = c(0.05, 1))


ggsave(paste0("figures/", species, "_combomaps.cont.png"), plot = combomaps.cont, 
       width = 14, height = 3, dpi = 450)
ggsave(paste0("figures/", species, "_combomaps.cont.svg"), plot = combomaps.cont, 
       width = 10, height = 10, dpi = 450)

# Change over time

# 1. Create binary maps (1 = Suitable, 0 = Unsuitable)
# use THRESH here as the 10p threshold
curr_binary <- current_prob.masked >= THRESH
futu_binary <- prob_2100.masked >= THRESH

ext_crop <- ext(-25, 55, -38, 38)

curr_binary <- crop(curr_binary, ext_crop)
futu_binary <- crop(futu_binary, ext_crop)

plot(curr_binary, col = c("khaki", "#d73027"), legend = FALSE)
plot(futu_binary, col = c("khaki", "#d73027"), legend = FALSE)

curr.plot = ggplot() +
  geom_spatraster(data = curr_binary) +
  scale_fill_manual(values = c("khaki", "#d73027"), na.value = "white") +
  coord_sf(xlim = c(-25, 55),
           ylim = c(-38, 38),
           expand = FALSE) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", 
       title = "Current",
       #subtitle = expression(italic(species))
       ) +
  theme(legend.position = "none")

future.plot = ggplot() +
  geom_spatraster(data = futu_binary) +
  scale_fill_manual(values = c("khaki", "#d73027"), na.value = "white") +
  coord_sf(xlim = c(-25, 55),
           ylim = c(-38, 38),
           expand = FALSE) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", 
       title = "2100",
       #subtitle = expression(italic(species))
       ) +
  theme(legend.position = "none")

comparison.bin.plots = gridExtra::grid.arrange(curr.plot, future.plot, ncol = 2)

ggsave(paste0("figures/", species, "_binary_comp.plots.png"), plot = comparison.bin.plots, 
       width = 10, height = 10, dpi = 450)

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


#################################################
# variable importance
#################################################

# read in all csv files that contain the name _var_importance.csv

file_list <- list.files(path = "figures", 
                        pattern = "var_importance", # Simplified pattern
                        full.names = TRUE,
                        ignore.case = TRUE)

# CHECK: If this returns 0, the path "figures" is wrong or empty
if(length(file_list) == 0) stop("No files found! Check if 'figures/' is the correct subfolder.")

# 3. Read and Name in one clean pipe (requires 'purrr' or just base R)
var_imp_list <- lapply(file_list, read.csv)

# 4. Extract clean species names
# This regex removes everything from the underscore onwards
names(var_imp_list) <- gsub("_.*", "", basename(file_list))

head(var_imp_list[[species_code]])

master_var_imp <- bind_rows(var_imp_list, .id = "species")

head(master_var_imp)

# Calculate both Mean and SD
avg_var_imp <- master_var_imp %>%
  group_by(species, algo, expl.var) %>%
  summarise(
    mean_imp = mean(var.imp, na.rm = TRUE),
    sd_imp   = sd(var.imp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(species, algo, desc(mean_imp))

# Preview the result
head(avg_var_imp)

ggplot(avg_var_imp, aes(x = reorder(expl.var, mean_imp), y = mean_imp, fill = algo)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_imp - sd_imp, ymax = mean_imp + sd_imp), 
                width = 0.2, position = position_dodge(0.9)) +
  coord_flip() +
  theme_classic() +
  facet_wrap(~species) +
  labs(title = "Variable Importance with Uncertainty (SD)", 
       x = "Variable", y = "Mean Importance")

var_importance_heatmap = ggplot(avg_var_imp, aes(x = algo, y = expl.var, fill = mean_imp)) +
  geom_tile(color = "white") +
  facet_wrap(~species, scales = "free_y") +
  scale_fill_viridis_c(option = "magma", name = "Importance") +
  labs(title = "Model Consensus on Variable Importance",
       x = "Model", y = "Bioclimatic Variable") +
  theme_minimal() +
  theme(
    strip.text = ggtext::element_markdown(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) ;var_importance_heatmap

ggsave("figures/var_importance_heatmap.png", plot = var_importance_heatmap, 
       width = 10, height = 6, dpi = 450)
write.csv(x = avg_var_imp, file = "figures/full_var_importances.txt", row.names = F)

# 1. Clean and format the species names
avg_var_imp <- avg_var_imp %>%
  mutate(species = species %>% 
           # 1. Replace the dot with a dot + space
           str_replace("\\.", ". ") %>% 
           # 2. Ensure the first letter is uppercase
           # Use sub() to just hit the very first character
           sub("(^.)", "\\U\\1", ., perl = TRUE) %>% 
           # 3. Wrap in markdown italics
           paste0("*", ., "*"))

# 2. Get the top 3 with the new formatted names
top_3_summary <- avg_var_imp %>%
  group_by(species, algo) %>%
  slice_max(mean_imp, n = 3) %>%
  mutate(importance_label = paste0(round(mean_imp, 3), " (± ", round(sd_imp, 3), ")"))

# Preview the result
# Now 'species' will look like "*D. aberiae*"
head(top_3_summary)

var_importance.plot = ggplot(top_3_summary, aes(x = reorder(expl.var, mean_imp), y = mean_imp, fill = expl.var)) +
  geom_col(show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_imp - sd_imp, ymax = mean_imp + sd_imp), width = 0.2) +
  coord_flip() +
  facet_grid(species ~ algo, scales = "free_y") + # 'free_y' allows different top vars per panel
  labs(title = "Top 3 bioclimatic predictors by mealybug species and model",
       x = "Bioclimatic Variable",
       y = "Mean Importance (± SD)") +
  scale_fill_viridis_d(option = "mako", begin = 0.2, end = 0.8) + # 'mako' or 'plasma' work well
  theme_bw() +
  theme(strip.text.y = ggtext::element_markdown())

ggsave("figures/var_importance_full.png", plot = var_importance.plot, 
       width = 14, height = 8, dpi = 450)
