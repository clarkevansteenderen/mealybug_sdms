
# PRESENCE POINTS
PRESPOINTS = readRDS(file = "data/PRESPOINTS.ALLSPP.rds")
PRESPOINTS = map(PRESPOINTS, ~ dplyr::select(.x, longitude, latitude))

# BACKGROUND POINTS
BACKPOINTS = readRDS(file = "data/BACKPOINTS.ALLSPP.rds")
BACKPOINTS = map(BACKPOINTS, ~ dplyr::select(.x, longitude, latitude))

names(PRESPOINTS)

# write out RDS files for each species

# For presence points
purrr::iwalk(PRESPOINTS, ~ saveRDS(.x, file = paste0("data/PRES_", .y, ".rds")))

# For background points
purrr::iwalk(BACKPOINTS, ~ saveRDS(.x, file = paste0("data/BACK_", .y, ".rds")))

# quick check
readRDS("data/PRES_n.viridis.rds")

# p.burnerae not working on hpc -> try here

p.burnerae.pres = PRESPOINTS$p.burnerae
p.burnerae.back = BACKPOINTS$p.burnerae
p.burnerae.preds = terra::rast("data/p.burnerae.tif")

list_settings = list(
  fc = c("L","Q","H","LQH"), 
  rm = c(1:5)
)

options(java.parameters = "-Xmx16000m")

tuning_results = 
  ENMeval::ENMevaluate(
    occs = p.burnerae.pres,
    envs = p.burnerae.preds,
    bg = p.burnerae.back,
    tune.args = list_settings, 
    partitions = "block",
    algorithm = "maxent.jar",
    doClamp = FALSE,
    parallel = F
  )
