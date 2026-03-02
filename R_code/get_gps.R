##################################
# Mealybug species GBIF data
##################################

# A list of species names
species_list = c("Pseudococcus longispinus",
                  "Planococcus citri",
                  "Paracoccus burnerae",
                  "Pseudococcus calceolariae",
                  "Nipaecoccus viridis",
                  "Delottococcus aberiae")

# Get keys for the list of species
results_list = rgbif::name_backbone_checklist(name_data = species_list)
species_keys = results_list$speciesKey

# download from GBIF
gbif_download = rgbif::occ_download(
  # exclude any records with geospatial issues
  rgbif::pred("hasGeospatialIssue", FALSE),
  # keep only records with available GPS coordinates
  rgbif::pred("hasCoordinate", TRUE),
  # remove absent records
  rgbif::pred("occurrenceStatus","PRESENT"), 
  # automatically looks for unique species keys (no duplication by default)
  rgbif::pred_in("speciesKey", species_keys), 
  format = "SIMPLE_CSV",
  user = "clarke.vansteenderen",
  pwd = "roxie2@!",
  email = "vsteenderen@gmail.com"
)

rgbif::occ_download_wait(gbif_download, 
                         quiet = FALSE, 
                         status_ping = 3) 

# <<gbif download metadata>>
# Status: SUCCEEDED
# DOI: 10.15468/dl.yqhhmd
# Format: SIMPLE_CSV
# Download key: 0061000-250717081556266
# Created: 2025-07-28T11:58:03.599+00:00
# Modified: 2025-07-28T12:00:13.649+00:00
# Download link: https://api.gbif.org/v1/occurrence/download/request/0061000-250717081556266.zip
# Total records: 3032

result = rgbif::occ_download_get(key = gbif_download, 
                                 overwrite = TRUE, 
                                 path = paste0("data/"))
