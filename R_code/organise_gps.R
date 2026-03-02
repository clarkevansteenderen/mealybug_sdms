
# get the lat and lon values from each sheet

#######################################################################
# Reyard's data
#######################################################################

sheet_names = readxl::excel_sheets("data/Mealybugs coordinates.xlsx")

# Read each sheet and add species name as a new column
combined_data_reyard = lapply(sheet_names, function(sheet) {
  data = read_excel("data/Mealybugs coordinates.xlsx", sheet = sheet)
  
  # Make sure it has 'lat' and 'lon' columns (case-insensitive match)
  data = data %>% rename_with(tolower)
  if (all(c("latitude", "longitude") %in% names(data))) {
    data$species = sheet  # Add species name
    return(data)
  } else {
    warning(paste("Sheet", sheet, "skipped: missing lat/lon columns"))
    return(NULL)
  }
}) %>%
  bind_rows()  # Combine into one data frame

combined_data_reyard$source = "literature"

#######################################################################
# Deric's data
#######################################################################

# function to convert to decimal degree format

dms_to_dd = function(dms) {
  # Updated pattern: allow decimals in seconds
  parts = stringr::str_match(dms, "(\\d+)°(\\d+)'(\\d+(?:\\.\\d+)?)[\"′]?(\\w)")
  deg = as.numeric(parts[,2])
  min = as.numeric(parts[,3])
  sec = as.numeric(parts[,4])
  dir = parts[,5]
  
  dd = deg + min / 60 + sec / 3600
  if (dir %in% c("S", "W")) dd = -dd
  return(dd)
}

sheet_names_deric = readxl::excel_sheets("data/deric_coordinates.xlsx") 

# Read each sheet and add species name as a new column
combined_data_deric = lapply(sheet_names_deric, function(sheet) {
  data = read_excel("data/deric_coordinates.xlsx", sheet = sheet)
  
  # Make sure it has 'lat' and 'lon' columns (case-insensitive match)
  data = data %>% rename_with(tolower)
  if (all(c("latitude", "longitude") %in% names(data))) {
    data$species = sheet  # Add species name
    return(data)
  } else {
    warning(paste("Sheet", sheet, "skipped: missing lat/lon columns"))
    return(NULL)
  }
}) %>%
  bind_rows()  # Combine into one data frame

# convert to decimal format
combined_data_deric = combined_data_deric %>%
  dplyr::mutate(
    latitude = sapply(latitude, dms_to_dd),
    longitude = sapply(longitude, dms_to_dd)
  )

combined_data_deric$source = "field"

MEALY.DATA = rbind(combined_data_reyard, combined_data_deric)

############################
# READ IN GBIF DOWNLOAD
###########################

GBIF.DATA = read.csv("data/0061000-250717081556266/0061000-250717081556266.csv",
                     sep = "\t") %>%
  dplyr::select(
    latitude = decimalLatitude,
    longitude = decimalLongitude,
    species
  )

GBIF.DATA$source = "GBIF"

# bind into one DF, and remove duplicate records
SP.DATA = rbind(MEALY.DATA, GBIF.DATA) %>%
  dplyr::distinct(longitude, latitude, .keep_all= TRUE)

SP.DATA = split(SP.DATA, SP.DATA$species)

# Convert one of your environmental predictors to a raster
r = raster::raster(pred_clim_current[[1]])

# Apply thinning to each species
SP.DATA.THIN = lapply(SP.DATA, function(df) {
  xy = df %>%
    dplyr::select(longitude, latitude) %>%
    as.data.frame()
  
  set.seed(2012)
  thinned_xy = dismo::gridSample(xy = xy, r = r, n = 1)
  
  # Filter df to keep only the rows that match thinned coordinates
  df_thinned = dplyr::semi_join(df, thinned_xy, by = c("longitude", "latitude"))
  
  return(df_thinned)
})


# get a summary table of records per species in original
SP.DATA.overview = data.frame(
  species = names(SP.DATA),
  rows = sapply(SP.DATA, nrow)
)
rownames(SP.DATA.overview) = NULL

SP.DATA.overview

# get a summary table of records per species in thinned
SP.DATA.THIN.overview = data.frame(
  species = names(SP.DATA.THIN),
  rows = sapply(SP.DATA.THIN, nrow)
)
rownames(SP.DATA.THIN.overview) = NULL

SP.DATA.THIN.overview

# Save data list to file
saveRDS(SP.DATA.THIN, "data/SP.DATA.THIN.rds")

########################
# plotting GPS
########################

world_map = rnaturalearth::ne_countries(
  scale = "medium", 
  returnclass = "sf"
) 

df_all = bind_rows(SP.DATA.THIN, .id = "species")

my_species_colors = c(
  "#1b9e77",  # Teal green
  "#d95f02",  # Orange
  "#7570b3",  # Purple
  "#e7298a",  # Pink/red
  "#66a61e",  # Olive green
  "#e6ab02"   # Mustard yellow
)

global_distr = ggplot() +
  # Add raster layer of world map 
  geom_sf(data = world_map, alpha = 0.5) +
  # Add GPS points 
  geom_point(
    data = df_all, 
    size = 1,
    aes(
      x = longitude, 
      y = latitude,
      color = species
    )
  )  +
  scale_color_manual(values = my_species_colors) +
  # Set world map CRS 
  coord_sf(
    crs = 4326,
    ylim = c(-58, 90),  # Clip below -52 degrees latitude
    expand = FALSE
  ) + 
  xlab("Longitude") + 
  ylab("Latitude") +
  guides(color = guide_legend(title = "Species"))

global_distr
