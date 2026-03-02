
# get gps coordinates for each species (also a list)
SP.DATA.THIN = readRDS(file = "data/PRESPOINTS.ALLSPP.rds")
PREDICTOR.LIST.ALLSPP

# read in the koppen-geiger (KG) data
kg_layer = sf::st_read("climate_layers/koppen_geiger")

# Reproject KG layer
geo_proj = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kg_layer = sf::st_transform(kg_layer, geo_proj)

names(SP.DATA.THIN)
names(PREDICTOR.LIST.ALLSPP)

BG.AREAS = list()

# Loop through each species in SP.DATA.THIN
for (k in 1:length(SP.DATA.THIN)) {
  
  cat("Processing species", k, "...\n")
  
DF = as.data.frame(SP.DATA.THIN[[k]])
  
  # Coerce GPS records into SpatialPointsDataFrame
  records_spatial = sp::SpatialPointsDataFrame(
    coords = cbind(DF$longitude, DF$latitude),
    data = DF,
    proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
  )
  
  # Select KG ecoregions in which there is at least one GPS record
  kg_sp = as(kg_layer, "Spatial")
  kg_contain = kg_sp[records_spatial, ]
  
  # Define background area by masking predictor layers to just the KG zones with at 
  # least 1 GPS record
  # Convert the KG zones containing GPS records back into an 'sf' object
  crs_wgs84 = CRS(SRS_string = "EPSG:4326")
  kg_contain = sf::st_as_sf(kg_contain)
  kg_contain = sf::st_set_crs(kg_contain, crs_wgs84)
  
  bg_area = terra::mask(PREDICTOR.LIST.ALLSPP[[k]], kg_contain)  
  BG.AREAS[[k]] = bg_area
}

BG.AREAS
names(BG.AREAS) = c("d.aberiae", "n.viridis", "p.burnerae",
                      "p.citri", "p.calceolariae", "p.longispinus")

# check
terra::plot(BG.AREAS$d.aberiae$`wc2.1_2.5m_bioc_AWI-CM-1-1-MR_ssp245_2021-2040_1`)
terra::plot(BG.AREAS$n.viridis$`wc2.1_2.5m_bioc_AWI-CM-1-1-MR_ssp245_2021-2040_1`)

# get the number of pres points per species
species_row_counts = sapply(SP.DATA.THIN, nrow)
species_row_counts[[1]]


BG.POINTS = list()

# select background points

for(k in 1:length(BG.AREAS)){
  
  set.seed(2023)
  
  bg_points = terra::spatSample(
    x = BG.AREAS[[k]],        # Raster of background area to sample points from 
    size = 10000,       # How many background points do we want?
    method = "random",  # Random points
    replace = FALSE,    # Sample without replacement
    na.rm = TRUE,       # Remove background points that have NA climate data
    as.df = TRUE,       # Return background points as data.frame object
    xy = TRUE           # Return lat/lon values for each background point
  ) %>%
    dplyr::rename(longitude = x, latitude = y)
  
  BG.POINTS[[k]] = bg_points

}

names(BG.POINTS) = c("d.aberiae", "n.viridis", "p.burnerae",
                    "p.citri", "p.calceolariae", "p.longispinus")

head(BG.POINTS)
PREDICTOR.LIST.ALLSPP

# save to file
saveRDS(object = BG.POINTS, file = "data/BACKPOINTS.ALLSPP.rds")
