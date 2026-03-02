source("setup.R")

# reduced predictor files per species
d.aberiae.preds = terra::rast("data/d.aberiae_current.tif") 
n.viridis.preds = terra::rast("data/n.viridis_current.tif")  
p.burnerae.preds = terra::rast("data/p.burnerae_current.tif") 
p.citri.preds = terra::rast("data/p.citri_current.tif")  
p.calceolariae.preds = terra::rast("data/p.calceolariae_current.tif")  
p.longispinus.preds = terra::rast("data/p.longispinus_current.tif") 

PREDICTOR.LIST.ALLSPP = list(d.aberiae.preds, n.viridis.preds, p.burnerae.preds,
p.citri.preds, p.calceolariae.preds, p.longispinus.preds)
names(PREDICTOR.LIST.ALLSPP) = c("d.aberiae", "n.viridis", "p.burnerae",
                          "p.citri", "p.calceolariae", "p.longispinus")

# get gps coordinates for each species (also a list)
SP.DATA.THIN = readRDS(file = "data/SP.DATA.THIN.rds")

names(PREDICTOR.LIST.ALLSPP$d.aberiae)
SP.DATA.THIN$`Delottococcus aberiae`

# Create an empty list to hold data frames for each species
allSpeciesPresPoints = list()

# Loop through each species
for (i in seq_along(PREDICTOR.LIST.ALLSPP)) {
  
  # Extract predictor raster and occurrence data
  predictor_raster = PREDICTOR.LIST.ALLSPP[[i]]
  species_points = SP.DATA.THIN[[i]]
  
  # Extract environmental data at occurrence points
  env_data = as.data.frame(raster::extract(predictor_raster, 
                          cbind(species_points$longitude, species_points$latitude)))
  
  # Combine with GPS data
  species_data = cbind(species_points, env_data)
  
  # Add species name (assumes SP.DATA.THIN[[i]] includes species name or use a separate name list)
  species_name = unique(species_points$species)[1]
  species_data = species_data %>%
    dplyr::mutate(species = species_name) %>%
    dplyr::select(species, everything()) %>%
    na.omit()
  
  # Optional: assign projection if needed (skip if you're not using spatial ops after this)
  coordinates(species_data) = ~ longitude + latitude
  crs(species_data) = "EPSG:4326"
  
  # Store the result
  allSpeciesPresPoints[[i]] = as.data.frame(species_data)
}

names(allSpeciesPresPoints) = c("d.aberiae", "n.viridis", "p.burnerae",
                                "p.citri", "p.calceolariae", "p.longispinus")

head(allSpeciesPresPoints$d.aberiae)
allSpeciesPresPoints

saveRDS(object = allSpeciesPresPoints, file = "data/PRESPOINTS.ALLSPP.rds")

PRESPOINTS = readRDS(file = "data/PRESPOINTS.ALLSPP.rds")


# Check points on a map

world_map <- rnaturalearth::ne_countries(
  scale = "medium", 
  returnclass = "sf"
) 

for(k in 1:length(PRESPOINTS)){
  
  global_distr = ggplot() +
    # Add raster layer of world map 
    geom_sf(data = world_map, alpha = 0.5) +
    # Add GPS points 
    geom_point(
      data = PRESPOINTS[[k]], 
      size = 0.5,
      colour = "darkred",
      aes(
        x = longitude, 
        y = latitude
      )
    )  +
    # Set world map CRS 
    coord_sf(
      crs = 4326,
      ylim = c(-58, 90),  # Clip below -52 degrees latitude
      expand = FALSE
    ) + 
    xlab("Longitude") + 
    ylab("Latitude") +
    ggtitle(names(PRESPOINTS)[k])
  
  ggsave(plot = global_distr, filename = paste0("figures/", names(PRESPOINTS)[k], "_global.png"),
         width = 6, height = 4, dpi = 350)
  
}

global_distr

################################
# REMOVE INCORRECT GPS POINTS
################################

# N. viridis has a point in Russia and north China
PRESPOINTS$n.viridis = PRESPOINTS$n.viridis %>%
  dplyr::filter(latitude < 38)

# P. citri - remove point in Norway
PRESPOINTS$p.citri = PRESPOINTS$p.citri %>%
  dplyr::filter(latitude < 55)

# P. citri - remove point in US, close to Canada
PRESPOINTS$p.citri = PRESPOINTS$p.citri %>%
  dplyr::filter(!(latitude > 40 & longitude < -80))

# P. longispinus - remove Alaska and norhtern Europe
PRESPOINTS$p.longispinus = PRESPOINTS$p.longispinus %>%
  dplyr::filter(latitude <= 60)


saveRDS(object = PRESPOINTS, file = "data/PRESPOINTS.ALLSPP.rds")

# check again

PRESPOINTS = readRDS(file = "data/PRESPOINTS.ALLSPP.rds")

for(k in 1:length(PRESPOINTS)){
  
  global_distr = ggplot() +
    # Add raster layer of world map 
    geom_sf(data = world_map, alpha = 0.5) +
    # Add GPS points 
    geom_point(
      data = PRESPOINTS[[k]], 
      #data = PRESPOINTS$p.citri,
      size = 0.5,
      colour = "darkred",
      aes(
        x = longitude, 
        y = latitude
      )
    )  +
    # Set world map CRS 
    coord_sf(
      crs = 4326,
      ylim = c(-58, 90),  # Clip below -52 degrees latitude
      expand = FALSE
    ) + 
    ggtitle(names(PRESPOINTS)[k]) +
    xlab("Longitude") + 
    ylab("Latitude") 
   
  
  ggsave(plot = global_distr, filename = paste0("figures/", names(PRESPOINTS)[k], "_global.png"),
         width = 6, height = 4, dpi = 350)
  
}

summary(PRESPOINTS$p.calceolariae)
summary(PRESPOINTS$p.citri)
summary(PRESPOINTS$p.longispinus)
