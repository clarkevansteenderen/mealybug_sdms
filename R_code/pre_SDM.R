#################################
# Multicollinearity testing
#################################

source("get_climate.R")

# all 6 species in a list
SP.DATA.THIN = readRDS("data/SP.DATA.THIN.rds")

###############################################################################
# head(SP.DATA.THIN$`Delottococcus aberiae`)

# Loop through each species

multicol_results = lapply(names(SP.DATA.THIN), function(sp_name) {
  cat("\nProcessing species:", sp_name, "\n")
  
  # Extract just lat/lon
  sp_data = SP.DATA.THIN[[sp_name]]
  gps_coords = sp_data[, c("longitude", "latitude")]  # adjust column names if needed
  
  # Convert to SpatVector
  sp_points = terra::vect(gps_coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  
  # Extract climate data
  clim_sp = terra::extract(pred_clim_2030, sp_points) %>% na.omit() %>%
    dplyr::select(-ID) # drop the ID column
  
  clim_sp_corr = clim_sp
  new_colnames = gsub(".*_(\\d+)$", "\\1", names(clim_sp_corr))
  names(clim_sp_corr) = new_colnames
  clim_sp_corr = clim_sp_corr[, gtools::mixedsort(names(clim_sp_corr))]
  
  cor_mat = cor(clim_sp_corr)
  
  cor_mat_plot = cor_mat
  cor_mat_plot[cor_mat_plot == 1] = 0
  corrplot::corrplot(cor_mat_plot, type = "upper", order = "original",
                     tl.col = "black", tl.srt = 0,
                     main = paste("Corrplot for", sp_name))
  
  CORR.VAL = 0.7
  cor_mat_strong = cor_mat
  cor_mat_strong[cor_mat_strong < CORR.VAL & cor_mat_strong > -CORR.VAL] = 0
  
  network = igraph::graph_from_adjacency_matrix(cor_mat_strong, weighted = TRUE,
                                                 mode = "undirected", diag = FALSE)
  plot(network, vertex.color = "lightgreen", main = paste("Network for", sp_name))
  
  list(
    species = sp_name,
    correlation_matrix = cor_mat,
    strong_corr_matrix = cor_mat_strong,
    network = network
  )
})

multicol_results[[6]]$species                 # Species name
multicol_results[[1]]$correlation_matrix      # Full correlation matrix
multicol_results[[1]]$strong_corr_matrix      # Thresholded correlation matrix
plot(multicol_results[[6]]$network, vertex.color = "lightgreen")  # igraph object

#################################
# Select bio clim variables
#################################

#########################
# Delottococcus aberiae
#########################

d.aberiae.climvars_r2 = c(
  "bio_1",
  "bio_5",
  "bio_6",
  "bio_8",
  "bio_9",
  "bio_12"
)

# reduce our original full set of predictor variables to the ones specified above
reduced_preds = terra::subset(x = pred_clim_2030, subset = d.aberiae.climvars_r2)

# extract the data for each GPS record
clim_sp_reduced = terra::extract(
  x = reduced_preds,         
  y = dplyr::select(SP.DATA.THIN$`Delottococcus aberiae`, longitude, latitude)        
) %>%
  dplyr::select(-ID) %>% # drop the ID column
  na.omit() # remove any NA values

# a common approach is to run a follow-up test, using a Variance Inflation Factor method (VIF)
# this is another method of assessing multicollinearity between variables, where a threshold value of 5
# indicates moderate collinearity
usdm::vifstep(clim_sp_reduced, th = 5)

# let's remove bio 5

d.aberiae.climvars_r2 = c(
  "bio_1",
  "bio_6",
  "bio_8",
  "bio_9",
  "bio_12"
)

# reduce our original full set of predictor variables to the ones specified above
d.aberiae.preds.current = terra::subset(x = pred_clim_current, 
                                subset = d.aberiae.climvars_r2)

d.aberiae.preds.2030 = terra::subset(x = pred_clim_2030, 
                                        subset = d.aberiae.climvars_r2)

d.aberiae.preds.2050 = terra::subset(x = pred_clim_2050, 
                                        subset = d.aberiae.climvars_r2)

d.aberiae.preds.2070 = terra::subset(x = pred_clim_2070, 
                                        subset = d.aberiae.climvars_r2)

d.aberiae.preds.2100 = terra::subset(x = pred_clim_2100, 
                                        subset = d.aberiae.climvars_r2)

#########################
# Nipaecoccus viridis
#########################

n.viridis.climvars_r2 = c(
  "bio_1",
  "bio_5",
  "bio_8",
  "bio_12",
  "bio_19"
)

# reduce our original full set of predictor variables to the ones specified above
reduced_preds = terra::subset(x = pred_clim_2030, subset = n.viridis.climvars_r2)

# extract the data for each GPS record
clim_sp_reduced = terra::extract(
  x = reduced_preds,         
  y = dplyr::select(SP.DATA.THIN$`Nipaecoccus viridis`, longitude, latitude)        
) %>%
  dplyr::select(-ID) %>% # drop the ID column
  na.omit() # remove any NA values

# a common approach is to run a follow-up test, using a Variance Inflation Factor method (VIF)
# this is another method of assessing multicollinearity between variables, where a threshold value of 5
# indicates moderate collinearity
usdm::vifstep(clim_sp_reduced, th = 5)

# don't remove any bios, as bio12 is important

# reduce our original full set of predictor variables to the ones specified above
n.viridis.preds.current = terra::subset(x = pred_clim_current, 
                                        subset = n.viridis.climvars_r2)

n.viridis.preds.2030 = terra::subset(x = pred_clim_2030, 
                                     subset = n.viridis.climvars_r2)

n.viridis.preds.2050 = terra::subset(x = pred_clim_2050, 
                                     subset = n.viridis.climvars_r2)

n.viridis.preds.2070 = terra::subset(x = pred_clim_2070, 
                                     subset = n.viridis.climvars_r2)

n.viridis.preds.2100 = terra::subset(x = pred_clim_2100, 
                                     subset = n.viridis.climvars_r2)


#########################
# Paracoccus burnerae
#########################

p.burnerae.climvars_r2 = c(
  "bio_1",
  "bio_5",
  "bio_12",
  "bio_18",
  "bio_19"
)

# reduce our original full set of predictor variables to the ones specified above
p.burnerae.preds.current = terra::subset(x = pred_clim_current, 
                                        subset = p.burnerae.climvars_r2)

p.burnerae.preds.2030 = terra::subset(x = pred_clim_2030, 
                                     subset = p.burnerae.climvars_r2)

p.burnerae.preds.2050 = terra::subset(x = pred_clim_2050, 
                                     subset = p.burnerae.climvars_r2)

p.burnerae.preds.2070 = terra::subset(x = pred_clim_2070, 
                                     subset = p.burnerae.climvars_r2)

p.burnerae.preds.2100 = terra::subset(x = pred_clim_2100, 
                                     subset = p.burnerae.climvars_r2)

#########################
# Planococcus citri
#########################

p.citri.climvars_r2 = c(
  "bio_1",
  "bio_2",
  "bio_3",
  "bio_5",
  "bio_6",
  "bio_8",
  "bio_9",
  "bio_12",
  "bio_17",
  "bio_19"
)

# reduce our original full set of predictor variables to the ones specified above
reduced_preds = terra::subset(x = pred_clim_2030, subset = p.citri.climvars_r2)

# extract the data for each GPS record
clim_sp_reduced = terra::extract(
  x = reduced_preds,         
  y = dplyr::select(SP.DATA.THIN$`Planococcus citri`, longitude, latitude)        
) %>%
  dplyr::select(-ID) %>% # drop the ID column
  na.omit() # remove any NA values

# a common approach is to run a follow-up test, using a Variance Inflation Factor method (VIF)
# this is another method of assessing multicollinearity between variables, where a threshold value of 5
# indicates moderate collinearity
usdm::vifstep(clim_sp_reduced, th = 5)

# remove bio6

p.citri.climvars_r2 = c(
  "bio_1",
  "bio_2",
  "bio_3",
  "bio_5",
  "bio_8",
  "bio_9",
  "bio_12",
  "bio_17",
  "bio_19"
)

# reduce our original full set of predictor variables to the ones specified above
p.citri.preds.current = terra::subset(x = pred_clim_current, 
                                        subset = p.citri.climvars_r2)

p.citri.preds.2030 = terra::subset(x = pred_clim_2030, 
                                     subset = p.citri.climvars_r2)

p.citri.preds.2050 = terra::subset(x = pred_clim_2050, 
                                     subset = p.citri.climvars_r2)

p.citri.preds.2070 = terra::subset(x = pred_clim_2070, 
                                     subset = p.citri.climvars_r2)

p.citri.preds.2100 = terra::subset(x = pred_clim_2100, 
                                     subset = p.citri.climvars_r2)


#########################
# Pseudococcus calceolariae
#########################

p.calceolariae.climvars_r2 = c(
  "bio_1",
  "bio_5",
  "bio_8",
  "bio_9",
  "bio_12"
  
)

# reduce our original full set of predictor variables to the ones specified above
reduced_preds = terra::subset(x = pred_clim_2030, subset = p.calceolariae.climvars_r2)

# extract the data for each GPS record
clim_sp_reduced = terra::extract(
  x = reduced_preds,         
  y = dplyr::select(SP.DATA.THIN$`Pseudococcus calceolariae`, longitude, latitude)        
) %>%
  dplyr::select(-ID) %>% # drop the ID column
  na.omit() # remove any NA values

# a common approach is to run a follow-up test, using a Variance Inflation Factor method (VIF)
# this is another method of assessing multicollinearity between variables, where a threshold value of 5
# indicates moderate collinearity
usdm::vifstep(clim_sp_reduced, th = 5)

# remove bio5

p.calceolariae.climvars_r2 = c(
  "bio_1",
  "bio_8",
  "bio_9",
  "bio_12"
  
)

# reduce our original full set of predictor variables to the ones specified above
p.calceolariae.preds.current = terra::subset(x = pred_clim_current, 
                                        subset = p.calceolariae.climvars_r2)

p.calceolariae.preds.2030 = terra::subset(x = pred_clim_2030, 
                                     subset = p.calceolariae.climvars_r2)

p.calceolariae.preds.2050 = terra::subset(x = pred_clim_2050, 
                                     subset = p.calceolariae.climvars_r2)

p.calceolariae.preds.2070 = terra::subset(x = pred_clim_2070, 
                                     subset = p.calceolariae.climvars_r2)

p.calceolariae.preds.2100 = terra::subset(x = pred_clim_2100, 
                                     subset = p.calceolariae.climvars_r2)



#########################
# Pseudococcus longispinus
#########################

p.longispinus.climvars_r2 = c(
  "bio_1",
  "bio_7",
  "bio_10",
  "bio_12",
  "bio_8",
  "bio_19"
  
)

# reduce our original full set of predictor variables to the ones specified above
reduced_preds = terra::subset(x = pred_clim_2030, subset = p.longispinus.climvars_r2)

# extract the data for each GPS record
clim_sp_reduced = terra::extract(
  x = reduced_preds,         
  y = dplyr::select(SP.DATA.THIN$`Pseudococcus longispinus`, longitude, latitude)        
) %>%
  dplyr::select(-ID) %>% # drop the ID column
  na.omit() # remove any NA values

# a common approach is to run a follow-up test, using a Variance Inflation Factor method (VIF)
# this is another method of assessing multicollinearity between variables, where a threshold value of 5
# indicates moderate collinearity
usdm::vifstep(clim_sp_reduced, th = 5)

# don't remove any bios

# reduce our original full set of predictor variables to the ones specified above
p.longispinus.preds.current = terra::subset(x = pred_clim_current, 
                                        subset = p.longispinus.climvars_r2)

p.longispinus.preds.2030 = terra::subset(x = pred_clim_2030, 
                                     subset = p.longispinus.climvars_r2)

p.longispinus.preds.2050 = terra::subset(x = pred_clim_2050, 
                                     subset = p.longispinus.climvars_r2)

p.longispinus.preds.2070 = terra::subset(x = pred_clim_2070, 
                                     subset = p.longispinus.climvars_r2)

p.longispinus.preds.2100 = terra::subset(x = pred_clim_2100, 
                                     subset = p.longispinus.climvars_r2)



mealybug.preds.current= list(d.aberiae.preds.current, n.viridis.preds.current, 
                             p.burnerae.preds.current,
                     p.citri.preds.current, p.calceolariae.preds.current, 
                     p.longispinus.preds.current)

names(mealybug.preds.current) = c("d.aberiae", "n.viridis", "p.burnerae",
                           "p.citri", "p.calceolariae", "p.longispinus")

mealybug.preds.2030= list(d.aberiae.preds.2030, n.viridis.preds.2030, 
                             p.burnerae.preds.2030,
                             p.citri.preds.2030, p.calceolariae.preds.2030, 
                             p.longispinus.preds.2030)

names(mealybug.preds.2030) = c("d.aberiae", "n.viridis", "p.burnerae",
                                  "p.citri", "p.calceolariae", "p.longispinus")

names(mealybug.preds.2030$d.aberiae)

mealybug.preds.2050= list(n.viridis.preds.2050, 
                             
                             p.citri.preds.2050, p.calceolariae.preds.2050, 
                             p.longispinus.preds.2050)

names(mealybug.preds.2050) = c("n.viridis", 
                                  "p.citri", "p.calceolariae", "p.longispinus")

mealybug.preds.2070= list( n.viridis.preds.2070, 
                             
                             p.citri.preds.2070, p.calceolariae.preds.2070, 
                             p.longispinus.preds.2070)

names(mealybug.preds.2070) = c("n.viridis", 
                                  "p.citri", "p.calceolariae", "p.longispinus")

mealybug.preds.2100= list( n.viridis.preds.2100, 
                             
                             p.citri.preds.2100, p.calceolariae.preds.2100, 
                             p.longispinus.preds.2100)

names(mealybug.preds.2100) = c( "n.viridis", 
                                  "p.citri", "p.calceolariae", "p.longispinus")


# Loop through each named element and save predictors per species as .tif
for (i in seq_along(mealybug.preds.current)) {
  name = names(mealybug.preds.current)[i]
  rast = mealybug.preds.current[[i]]
  terra::writeRaster(rast, filename = paste0("data/", name, "_current.tif"), overwrite = TRUE)
}

for (i in seq_along(mealybug.preds.2030)) {
  name = names(mealybug.preds.2030)[i]
  rast = mealybug.preds.2030[[i]]
  terra::writeRaster(rast, filename = paste0("data/", name, "_2030.tif"), overwrite = TRUE)
}

for (i in seq_along(mealybug.preds.2050)) {
  name = names(mealybug.preds.2050)[i]
  rast = mealybug.preds.2050[[i]]
  terra::writeRaster(rast, filename = paste0("data/", name, "_2050.tif"), overwrite = TRUE)
}

for (i in seq_along(mealybug.preds.2070)) {
  name = names(mealybug.preds.2070)[i]
  rast = mealybug.preds.2070[[i]]
  terra::writeRaster(rast, filename = paste0("data/", name, "_2070.tif"), overwrite = TRUE)
}

for (i in seq_along(mealybug.preds.2100)) {
  name = names(mealybug.preds.2100)[i]
  rast = mealybug.preds.2100[[i]]
  terra::writeRaster(rast, filename = paste0("data/", name, "_2100.tif"), overwrite = TRUE)
}
