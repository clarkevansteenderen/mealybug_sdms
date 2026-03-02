#################################################
# MEALYBUG SDM WORKFLOW
# Clarke van Steenderen
# July 2025
#################################################

# set preferences and load libs
source("setup.R")

# function for creating warming tolerance maps
source("create_wt_map.R")

# get clim files -> current and future, koppen-geiger
source("get_climate.R")

# get all GBIF records
source("get_gps.R")

# get all GPS records into one data frame
source("organise_gps.R")

# warming tolerance calculations and maps
source("warming_tolerance.R")

# SDM prep
source("pre_SDM.R")

# get_prespoints
# get_backpoints

