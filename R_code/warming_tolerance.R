##############################################

# warming tolerance (WT) = CTmax − Tmax
# CTmax = critial thermal maximum
# Tmax = mean temperature of warmest quarter (BIO10)

##############################################

# e.g. species A's CTmax is 25 degrees, and Tmax is 20 degrees
# WT = 5 degrees (some margin to adapt to a warmer climate)

# spB's CTmax is 20 degrees, Tmax is 20 degrees
# WT = 0 degrees (no margin -> a warming climate will be bad)

# spC's CTmax is 40 degrees, Tmax is -40 degrees
# WT = 40 - (-40) = 80 degrees

# module load chpc/BIOMODULES R/4.2.0
# export LANG=en_US.UTF-8
# export LC_ALL=en_US.UTF-8

##############################################

source("setup.R")
source("create_wt_map.R")
source("get_climate.R")
library(emmeans)

CT.VALS = readxl::read_xlsx("data/Modelling data adult mealybugs CTL.xlsx",
                            sheet = 2) %>%
  janitor::clean_names()

head(CT.VALS)

CT.VALS$species = as.factor(CT.VALS$species)

levels(CT.VALS$species)

ct.glm = glm(data = CT.VALS, formula = ctmax ~ species)
summary(ct.glm)
anova(ct.glm, test="F")

DHARMa::plotResiduals(ct.glm)
DHARMa::plotQQunif(ct.glm)

em = emmeans(ct.glm, pairwise ~ species, adjust = "tukey")

cld = multcomp::cld(em$emmeans, Letters = letters)

cld

cld_letters <- cld %>%
  mutate(.group = gsub(" ", "", .group))

max_vals <- CT.VALS %>%
  group_by(species) %>%
  summarise(y_pos = max(ctmax) + 0.4)

cld_letters <- left_join(cld_letters, max_vals, by = "species")

CT.MEANS = CT.VALS %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_ct = mean(ctmax, na.rm = TRUE),
    sd_ct = sd(ctmax, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(desc(mean_ct))

CT.MEANS

ctmeans.plot = ggplot(CT.VALS, aes(x = species, y = ctmax)) +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Species", y = "Critical Thermal Maxima (°C)", title = "") +
  scale_y_continuous(breaks = seq(45, 50, 1)) +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1));ctmeans.plot

ggsave("figures/CTmeans.svg", plot = ctmeans.plot, 
       width = 6, height = 4, dpi = 350)

# boxplot


ctmeans.boxplot = ggplot(CT.VALS, aes(x = species, y = ctmax)) +
  geom_boxplot(fill = "gray85", color = "black", width = 0.6, outlier.shape = NA) +
  #geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +  # Optional: show individual points
  labs(
    x = "Species", 
    y = "Critical Thermal Maxima (°C)", 
    title = ""
  ) +
  scale_y_continuous(breaks = seq(45, 55, 1)) +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  stat_summary(fun = mean, geom = "point", shape = 23, fill = "black", size = 3, color = "black") +
  geom_text(data = cld_letters,
            aes(x = species, y = y_pos, label = .group),
            size = 4) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

ctmeans.boxplot 

ggsave("figures/CTmeans_new.png", plot = ctmeans.boxplot, 
       width = 5, height = 5, dpi = 350)

####################################################
# CLIMATE DATA
####################################################

africa_ext = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(continent == "Africa")

southernafrica_ext = rnaturalearth::ne_countries(scale = "medium",
                                                 returnclass = "sf") %>%
  dplyr::filter(name %in% c("South Africa", "Namibia", "Botswana", 
                            "Zimbabwe", "Mozambique", "Angola",
                            "Zambia", "Malawi", "Lesotho", "eSwatini"))

africa_vect = vect(africa_ext)
southern_africa_vect = vect(southernafrica_ext)

world_map = rnaturalearth::ne_countries(
  scale = "medium", 
  returnclass = "sf"
) 

world_vect = vect(world_map)

############################
# current climate
############################

#########
# GLOBAL
#########

WT.global.current.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_current$bio_10,
  ctmax_df = CT.MEANS,
  vect = world_vect,
  EXT = world_map,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 2,
  xlim = c(-180, 180),
  ylim = c(-60, 90),
  lower_bound = 10,
  upper_bound = 32,
  bin_width = 4
) ;WT.global.current.plot$plot

ggsave("figures/WT.global.current.png", plot = WT.global.current.plot$plot, 
       width = 8, height = 8, dpi = 450)
ggsave("figures/WT.global.current.svg", plot = WT.global.current.plot$plot, 
       width = 8, height = 8, dpi = 450)

#########
# AFRICA
#########

WT.Africa.current.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_current$bio_10,
  ctmax_df = CT.MEANS,
  vect = africa_vect,
  EXT = africa_ext,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 3,
  xlim = c(-25, 55),
  ylim = c(-38, 38),
  lower_bound = 10,
  upper_bound = 40,
  bin_width = 5
) ;WT.Africa.current.plot$plot

ggsave("figures/WT.Africa.current.png", plot = WT.Africa.current.plot$plot, 
       width = 10, height = 8, dpi = 450)
ggsave("figures/WT.Africa.current.svg", plot = WT.Africa.current.plot$plot, 
       width = 10, height = 8, dpi = 450)

###################
# SOUTHERN AFRICA
###################

WT.SAfrica.current.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_current$bio_10,
  ctmax_df = CT.MEANS,
  vect = southern_africa_vect,
  EXT = southernafrica_ext,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 2,
  xlim = c(10, 43),
  ylim = c(-35, -3),
  lower_bound = 18,
  upper_bound = 35,
  bin_width = 3
) ;WT.SAfrica.current.plot$plot

ggsave("figures/WT.SAfrica.current.png", plot = WT.SAfrica.current.plot$plot, 
       width = 8, height = 8, dpi = 450)
ggsave("figures/WT.SAfrica.current.svg", plot = WT.SAfrica.current.plot$plot, 
       width = 8, height = 8, dpi = 450)

############################
# 2030 climate
############################

#########
# GLOBAL
#########

WT.global.2030.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_2030$bio_10,
  ctmax_df = CT.MEANS,
  vect = world_vect,
  EXT = world_map,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 2,
  xlim = c(-180, 180),
  ylim = c(-60, 90),
  lower_bound = 10,
  upper_bound = 32,
  bin_width = 4
) ;WT.global.2030.plot$plot

ggsave("figures/WT.global.2030.png", plot = WT.global.2030.plot$plot, 
       width = 8, height = 8, dpi = 450)
ggsave("figures/WT.global.2030.svg", plot = WT.global.2030.plot$plot, 
       width = 8, height = 8, dpi = 450)

#########
# AFRICA
#########

range(pred_clim_2030$`wc2.1_2.5m_bioc_AWI-CM-1-1-MR_ssp245_2021-2040_10`)

WT.Africa.2030.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_2030$bio_10,
  ctmax_df = CT.MEANS,
  vect = africa_vect,
  EXT = africa_ext,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 2,
  xlim = c(-25, 55),
  ylim = c(-38, 38),
  lower_bound = 10,
  upper_bound = 32,
  bin_width = 4
) ;WT.Africa.2030.plot$plot

ggsave("figures/WT.Africa.2030.png", plot = WT.Africa.2030.plot$plot, 
       width = 8, height = 8, dpi = 450)
ggsave("figures/WT.Africa.2030.svg", plot = WT.Africa.2030.plot$plot, 
       width = 8, height = 8, dpi = 450)

###################
# SOUTHERN AFRICA
###################

WT.SAfrica.2030.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_2030$bio_10,
  ctmax_df = CT.MEANS,
  vect = southern_africa_vect,
  EXT = southernafrica_ext,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 2,
  xlim = c(10, 43),
  ylim = c(-35, -3),
  lower_bound = 18,
  upper_bound = 35,
  bin_width = 3
) ;WT.SAfrica.2030.plot$plot

ggsave("figures/WT.SAfrica.2030.png", plot = WT.SAfrica.2030.plot$plot, 
       width = 8, height = 8, dpi = 450)
ggsave("figures/WT.SAfrica.2030.svg", plot = WT.SAfrica.2030.plot$plot, 
       width = 8, height = 8, dpi = 450)

############################
# 2050 climate
############################

#########
# GLOBAL
#########

WT.global.2050.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_2050$bio_10,
  ctmax_df = CT.MEANS,
  vect = world_vect,
  EXT = world_map,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 2,
  xlim = c(-180, 180),
  ylim = c(-60, 90),
  lower_bound = 10,
  upper_bound = 32,
  bin_width = 4
) ;WT.global.2050.plot$plot

ggsave("figures/WT.global.2050.png", plot = WT.global.2050.plot$plot, 
       width = 8, height = 8, dpi = 450)
ggsave("figures/WT.global.2050.svg", plot = WT.global.2050.plot$plot, 
       width = 8, height = 8, dpi = 450)

#########
# AFRICA
#########

WT.Africa.2050.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_2050$bio_10,
  ctmax_df = CT.MEANS,
  vect = africa_vect,
  EXT = africa_ext,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 3,
  xlim = c(-25, 55),
  ylim = c(-38, 38),
  lower_bound = 10,
  upper_bound = 32,
  bin_width = 4
) ;WT.Africa.2050.plot$plot

ggsave("figures/WT.Africa.2050.png", plot = WT.Africa.2050.plot$plot, 
       width = 10, height = 8, dpi = 450)
ggsave("figures/WT.Africa.2050.svg", plot = WT.Africa.2050.plot$plot, 
       width = 10, height = 8, dpi = 450)

###################
# SOUTHERN AFRICA
###################

WT.SAfrica.2050.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_2050$bio_10,
  ctmax_df = CT.MEANS,
  vect = southern_africa_vect,
  EXT = southernafrica_ext,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 2,
  xlim = c(10, 43),
  ylim = c(-35, -3),
  lower_bound = 18,
  upper_bound = 35,
  bin_width = 3
) ;WT.SAfrica.2050.plot$plot

ggsave("figures/WT.SAfrica.2050.png", plot = WT.SAfrica.2050.plot$plot, 
       width = 8, height = 8, dpi = 450)
ggsave("figures/WT.SAfrica.2050.svg", plot = WT.SAfrica.2050.plot$plot, 
       width = 8, height = 8, dpi = 450)

############################
# 2070 climate
############################

#########
# GLOBAL
#########

WT.global.2070.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_2070$bio_10,
  ctmax_df = CT.MEANS,
  vect = world_vect,
  EXT = world_map,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 2,
  xlim = c(-180, 180),
  ylim = c(-60, 90),
  lower_bound = 10,
  upper_bound = 32,
  bin_width = 4
) ;WT.global.2070.plot$plot

ggsave("figures/WT.global.2070.png", plot = WT.global.2070.plot$plot, 
       width = 8, height = 8, dpi = 450)
ggsave("figures/WT.global.2070.svg", plot = WT.global.2070.plot$plot, 
       width = 8, height = 8, dpi = 450)

#########
# AFRICA
#########

WT.Africa.2070.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_2070$bio_10,
  ctmax_df = CT.MEANS,
  vect = africa_vect,
  EXT = africa_ext,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 3,
  xlim = c(-25, 55),
  ylim = c(-38, 38),
  lower_bound = 10,
  upper_bound = 32,
  bin_width = 4
) ;WT.Africa.2070.plot$plot

ggsave("figures/WT.Africa.2070.png", plot = WT.Africa.2070.plot$plot, 
       width = 10, height = 8, dpi = 450)
ggsave("figures/WT.Africa.2070.svg", plot = WT.Africa.2070.plot$plot, 
       width = 8, height = 8, dpi = 450)

###################
# SOUTHERN AFRICA
###################

WT.SAfrica.2070.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_2070$bio_10,
  ctmax_df = CT.MEANS,
  vect = southern_africa_vect,
  EXT = southernafrica_ext,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 2,
  xlim = c(10, 43),
  ylim = c(-35, -3),
  lower_bound = 18,
  upper_bound = 35,
  bin_width = 3
) ;WT.SAfrica.2070.plot$plot

ggsave("figures/WT.SAfrica.2070.png", plot = WT.SAfrica.2070.plot$plot, 
       width = 8, height = 8, dpi = 450)
ggsave("figures/WT.SAfrica.2070.svg", plot = WT.SAfrica.2070.plot$plot, 
       width = 8, height = 8, dpi = 450)

############################
# 2100 climate
############################

#########
# GLOBAL
#########

WT.global.2100.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_2100$bio_10,
  ctmax_df = CT.MEANS,
  vect = world_vect,
  EXT = world_map,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 2,
  xlim = c(-180, 180),
  ylim = c(-60, 90),
  lower_bound = 10,
  upper_bound = 32,
  bin_width = 4
) ;WT.global.2100.plot$plot

ggsave("figures/WT.global.2100.png", plot = WT.global.2100.plot$plot, 
       width = 8, height = 8, dpi = 450)
ggsave("figures/WT.global.2100.svg", plot = WT.global.2100.plot$plot, 
       width = 8, height = 8, dpi = 450)

#########
# AFRICA
#########

WT.Africa.2100.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_2100$bio_10,
  ctmax_df = CT.MEANS,
  vect = africa_vect,
  EXT = africa_ext,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 3,
  xlim = c(-25, 55),
  ylim = c(-38, 38),
  lower_bound = 10,
  upper_bound = 32,
  bin_width = 4
) ;WT.Africa.2100.plot$plot

ggsave("figures/WT.Africa.2100.png", plot = WT.Africa.2100.plot$plot, 
       width = 10, height = 8, dpi = 450)
ggsave("figures/WT.Africa.2100.svg", plot = WT.Africa.2100.plot$plot, 
       width = 8, height = 8, dpi = 450)

###################
# SOUTHERN AFRICA
###################

WT.SAfrica.2100.plot = CREATE.WT.MAP(
  climate_layer = pred_clim_2100$bio_10,
  ctmax_df = CT.MEANS,
  vect = southern_africa_vect,
  EXT = southernafrica_ext,
  layer_label = "Warming Tolerance (°C)",
  facet_ncol = 2,
  xlim = c(10, 43),
  ylim = c(-35, -3),
  lower_bound = 18,
  upper_bound = 35,
  bin_width = 3
) ;WT.SAfrica.2100.plot$plot

ggsave("figures/WT.SAfrica.2100.png", plot = WT.SAfrica.2100.plot$plot, 
       width = 8, height = 8, dpi = 450)
ggsave("figures/WT.SAfrica.2100.svg", plot = WT.SAfrica.2100.plot$plot, 
       width = 8, height = 8, dpi = 450)


#################################################
# plot species across each time frame

# WT.SAfrica.current.plot$rasters$`Delottococcus aberiae`
# WT.SAfrica.2030.plot$rasters$`Delottococcus aberiae`
# WT.SAfrica.2050.plot$rasters$`Delottococcus aberiae`
# WT.SAfrica.2070.plot$rasters$`Delottococcus aberiae`

# 1. Get species list
# Manually re-assign list names from the internal raster names
names(WT.Africa.current.plot$rasters) <- sapply(WT.Africa.current.plot$rasters, names)
species_list <- names(WT.Africa.current.plot$rasters)
species_list =  species_list[species_list != "NULL"]

names(WT.Africa.current.plot$rasters) <- names(WT.Africa.current.plot$rasters) 
names(WT.Africa.2050.plot$rasters)    <- names(WT.Africa.current.plot$rasters) 
names(WT.Africa.2070.plot$rasters)    <- names(WT.Africa.current.plot$rasters)
names(WT.Africa.2100.plot$rasters)    <- names(WT.Africa.current.plot$rasters)

WT.Africa.current.plot$rasters <- WT.Africa.current.plot$rasters[names(WT.Africa.current.plot$rasters) != "NULL"]
WT.Africa.2050.plot$rasters    <- WT.Africa.2050.plot$rasters[names(WT.Africa.2050.plot$rasters) != "NULL"]
WT.Africa.2070.plot$rasters    <- WT.Africa.2070.plot$rasters[names(WT.Africa.2070.plot$rasters) != "NULL"]
WT.Africa.2100.plot$rasters    <- WT.Africa.2100.plot$rasters[names(WT.Africa.2100.plot$rasters) != "NULL"]


# 2. Function to convert a raster to df with metadata
raster_to_df <- function(ras, species, time) {
  # Only process if raster is not NULL and has values
  if (!is.null(ras) && terra::ncell(ras) > 0 && !is.null(terra::values(ras))) {
    df <- as.data.frame(ras, xy = TRUE, na.rm = TRUE)
    if (ncol(df) >= 3) {
      names(df)[3] <- "WT"  # rename the third column to WT
      df$species <- species
      df$time <- time
      return(df)
    }
  }
  return(NULL)  # Return nothing if raster is invalid
}


# 3. Iterate over species and times

# Southern Africa
# all_df <- purrr::map_dfr(species_list, function(sp) {
#   bind_rows(
#     raster_to_df(WT.SAfrica.2030.plot$rasters[[sp]], sp, "2030"),
#     raster_to_df(WT.SAfrica.2050.plot$rasters[[sp]], sp, "2050"),
#     raster_to_df(WT.SAfrica.2070.plot$rasters[[sp]], sp, "2070"),
#     raster_to_df(WT.SAfrica.2100.plot$rasters[[sp]], sp, "2100"),
#   )
# })

# Africa -> whole continent
all_df <- purrr::map_dfr(species_list, function(sp) {
  dplyr::bind_rows(
    raster_to_df(WT.Africa.current.plot$rasters[[sp]], sp, "Current"),
    raster_to_df(WT.Africa.2050.plot$rasters[[sp]], sp, "2050"),
    raster_to_df(WT.Africa.2070.plot$rasters[[sp]], sp, "2070"),
    raster_to_df(WT.Africa.2100.plot$rasters[[sp]], sp, "2100"),
  )
})

str(all_df)

# Parameters
lower_bound <- 10 # 18 and 35 for southern Africa
upper_bound <- 32 # 35
bin_width <- 4

# Create breaks dynamically
breaks <- c(-Inf, seq(lower_bound, upper_bound, by = bin_width), Inf)

# Create labels dynamically
labels <- c(
  paste0("<", lower_bound),
  paste(
    seq(lower_bound, upper_bound - bin_width, by = bin_width),
    seq(lower_bound + bin_width, upper_bound, by = bin_width) - 1,
    sep = "-"
  ),
  paste0(">", upper_bound)
)

# Your base_colors, adjusted for number of bins
base_colors <- c(
  "#d73027",  # deep red (< lower bound)
  "#fc8d59",  # orange-red
  "#fdae61",  # orange
  "#fee08b",  # yellow
  "#d9ef8b",  # light green
  "#a6d96a",  # medium green
  "#66bd63",  # darker green
  "#41ab5d",  # even darker green
  "#1a9850",  # dark green
  "#006837"   # darkest green (> upper bound)
)

# Adjust base_colors length to match number of labels (if needed)
n_bins <- length(labels)
if (n_bins <= length(base_colors)) {
  bin_colors <- base_colors[1:n_bins]
} else {
  # interpolate colors if more bins than colors
  bin_colors <- colorRampPalette(base_colors)(n_bins)
}
names(bin_colors) <- labels

# Apply binning
all_df <- all_df %>%
  dplyr::mutate(WT_bin = cut(WT, breaks = breaks, labels = labels, right = FALSE)) %>%
  dplyr::mutate(time = factor(time, levels = c("Current", "2050", "2070", "2100")))

# 1. Replace underscores with a space
all_df$species <- gsub("_", " ", all_df$species)

# Plot
combo.plot = ggplot(all_df, aes(x = x, y = y, fill = WT_bin)) +
  geom_raster() +
  geom_sf(data = africa_ext, fill = NA, color = "black", size = 0.3, inherit.aes = FALSE) +
  scale_fill_manual(
    values = bin_colors,
    name = "Warming Tolerance (°C)"
  ) +
  # southern Africa
  # coord_sf(xlim = c(10, 43),
  #          ylim = c(-36, -3),
  #          expand = FALSE) +
  # Africa
  coord_sf(xlim = c(-25, 55),
  ylim = c(-38, 38),
  expand = FALSE) +
  facet_grid(rows = vars(species), cols = vars(time)) +
  theme(
    strip.text = element_text(face = "italic", size = 7),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  theme(legend.position = "right") +
  xlab("Longitude") +
  ylab("Latitude") 

combo.plot


ggsave("figures/africa_combo.png", plot = combo.plot, 
       width = 10, height = 10, dpi = 450)
ggsave("figures/africa_combo.svg", plot = combo.plot, 
       width = 10, height = 10, dpi = 450)



# Calculate Min and 5th Percentile
wt_summary <- all_df %>%
  group_by(species, time) %>%
  summarise(
    min_wt = min(WT, na.rm = TRUE),
    p05_wt = quantile(WT, probs = 0.05, na.rm = TRUE),
    .groups = "drop" # Keeps the data frame clean for further use
  )

# View the result
print(wt_summary)
str(wt_summary)
wt_summary$species = as.factor(wt_summary$species)

ggplot(wt_summary, aes(x = time, y = p05_wt, group = 1)) +
  geom_line(color = "grey60", size = 1) +       # lines showing trend
  geom_point(color = "black", size = 2) +       # points at each period
  facet_wrap(~ species, ncol = 2) +  # panel per species
  xlab("Period") +
  ylab("5th Percentile Warming Tolerance (°C)") +
  #scale_y_continuous(limits = c(12, 18), breaks = seq(12, 18, by = 0.5)) +
  theme_classic() +
  theme(strip.text = element_text(face = "italic"))   # italic species names


WT.plot = ggplot(wt_summary, aes(x = time, group = 1)) +
  # both lines with legend
  geom_line(aes(y = p05_wt, color = "5th percentile WT"), size = 1) +
  geom_point(aes(y = p05_wt, color = "5th percentile WT"), size = 2) +
  geom_line(aes(y = min_wt, color = "Minimum WT"), size = 1, linetype = "dashed") +
  geom_point(aes(y = min_wt, color = "Minimum WT"), size = 2, shape = 17) +
  facet_wrap(~ species, ncol = 2) +
  xlab("Period") +
  ylab("Warming Tolerance (°C)") +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 2)) +
  scale_color_manual(
    name = "Metric",
    values = c("5th percentile WT" = "black", "Minimum WT" = "tomato")
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "italic"),
    legend.position = "right"
  ) +
  geom_vline(xintercept = which(levels(wt_summary$time) == "Current"), linetype = "dotted") +
  geom_vline(xintercept = which(levels(wt_summary$time) == "2100"), linetype = "dotted")

WT.plot
  
ggsave("figures/WT.plot.png", plot = WT.plot, 
       width = 8, height = 5, dpi = 450)
ggsave("figures/WT.plot.svg", plot = WT.plot, 
       width = 10, height = 10, dpi = 450)

wt_change <- wt_summary %>%
  pivot_wider(names_from = time, values_from = c(min_wt, p05_wt)) %>%
  mutate(
    min_abs_change = min_wt_Current - min_wt_2100,
    min_perc_reduction = (min_abs_change / min_wt_Current) * 100,
    
    p05_abs_change = p05_wt_Current - p05_wt_2100,
    p05_perc_reduction = (p05_abs_change / p05_wt_Current) * 100
  )

wt_change

wt_range <- wt_summary %>%
  dplyr::select(species, time, min_wt) %>%
  dplyr::filter(time %in% c("Current", "2100")) %>%
  tidyr::pivot_wider(names_from = time, values_from = min_wt) %>%
  dplyr::mutate(range_current_2100 = Current - `2100`)

wt_range
