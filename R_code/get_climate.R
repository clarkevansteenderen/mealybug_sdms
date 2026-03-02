#############################################################################
# CURRENT CLIMATE
#############################################################################

pred_clim_current = terra::rast( list.files(
  here::here("climate_layers/current/wc2.1_2.5m/") ,
  full.names = TRUE,
  pattern = '.tif'
))  

# set the CRS (coordinate reference system) projection for the current climate layers
terra::crs(pred_clim_current) = "epsg:4326"
terra::crs(pred_clim_current, describe = T)

names(pred_clim_current) = paste0("bio_", sub(".*_", "", names(pred_clim_current)))


#############################################################################
# FUTURE CLIMATES --> RUN THIS JUST ONCE TO DONWLOAD
#############################################################################

# Define the models you want to use
# models_to_use = c("MIROC6", "BCC-CSM2-MR", "CanESM5")

# Function to get an averaged ensemble for a specific time window
# get_future_ensemble = function(time_window, ssp_code = "245") {
#   
#   message(paste("Processing Ensemble for:", time_window))
#   
#   model_stacks = list()
#   
#   for (m in models_to_use) {
#     # Download data into a temporary directory to avoid folder clutter
#     dat = cmip6_world(
#       model = m, 
#       ssp = ssp_code, 
#       time = time_window, 
#       res = 2.5, 
#       var = "bioc", 
#       path = tempdir() 
#     )
#     model_stacks[[m]] = dat
#   }
#   
#   # Average the models (Mean of the 3 stacks)
#   ensemble = (model_stacks[[1]] + model_stacks[[2]] + model_stacks[[3]]) / 3
#   
#   # Standardize projection and names
#   terra::crs(ensemble) = "epsg:4326"
#   names(ensemble) = paste0("bio_", 1:19)
#   
#   return(ensemble)
# }

# generate the ensembles
#pred_clim_2030 = get_future_ensemble("2021-2040")
# pred_clim_2050 = get_future_ensemble("2041-2060")
# pred_clim_2070 = get_future_ensemble("2061-2080")
# pred_clim_2100 = get_future_ensemble("2081-2100")
# 
# # save them as TIFs so you don't have to download again
# terra::writeRaster(pred_clim_2050, here("climate_layers/future/ensemble_2050_ssp245.tif"))
# terra::writeRaster(pred_clim_2070, here("climate_layers/future/ensemble_2070_ssp245.tif"))
# terra::writeRaster(pred_clim_2100, here("climate_layers/future/ensemble_2100_ssp245.tif"))

read_future_env = function(file_path) {
  # Load the raster
  r = terra::rast(file_path)
  
  # Ensure names match your training data (e.g., bio_1, bio_2...)
  # This assumes you have 19 bioclim variables in order
  names(r) = paste0("bio_", 1:19)
  
  # Set the CRS just in case it wasn't preserved
  terra::crs(r) = "epsg:4326"
  
  return(r)
}

#pred_clim_2030 = read_future_env(here::here("climate_layers/future/ensemble_2030_ssp245.tif"))
pred_clim_2050 = read_future_env(here::here("climate_layers/future/ensemble_2050_ssp245.tif"))
pred_clim_2070 = read_future_env(here::here("climate_layers/future/ensemble_2070_ssp245.tif"))
pred_clim_2100 = read_future_env(here::here("climate_layers/future/ensemble_2100_ssp245.tif"))

##############################################################################
# KOPPEN-GEIGER DATA
##############################################################################

# read in the koppen-geiger (KG) data
kg_layer = sf::st_read(dsn = "climate_layers/koppen_geiger")

# Reproject KG layer
geo_proj = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kg_layer = sf::st_transform(kg_layer, geo_proj)

#############################################################################