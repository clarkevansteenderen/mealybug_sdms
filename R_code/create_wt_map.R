
# Function to calculate and plot warming tolerance maps
CREATE.WT.MAP <- function(climate_layer,
                          ctmax_df,
                          vect,
                          EXT,
                          layer_label = "Warming Tolerance (°C)",
                          facet_ncol = 2,
                          xlim = c(-25, 55),
                          ylim = c(-38, 38),
                          lower_bound = 10,
                          upper_bound = 55,
                          bin_width = 5) {
  
  # Crop and mask climate raster to Africa
  Africa.raster <- terra::crop(climate_layer, vect) %>%
    terra::mask(vect)
  
  # Compute warming tolerance rasters per species
  WT.rasters = list()
  
  for (i in 1:nrow(ctmax_df)) {
    species_name <- ctmax_df$species[i]
    ctmax <- ctmax_df$mean_ct[i]
    
    # 1. Skip if CTmax itself is NA
    if (is.na(ctmax)) {
      cat("  ⚠️ Skipping", species_name, "- CTmax value is NA\n")
      next
    }
    
    cat("Processing", species_name, "| CTmax:", ctmax, "\n")
    
    # 2. Perform math
    wt_raster <- ctmax - Africa.raster
    
    # 3. Check if the resulting raster is entirely empty
    # (Using global to check if all values are NA)
    if (all(is.na(terra::values(wt_raster)))) {
      cat("  ❌ Warning:", species_name, "produced an all-NA raster.\n")
    }
    
    names(wt_raster) <- gsub(" ", "_", species_name)
    WT.rasters[[species_name]] <- wt_raster
  }
  
  WT.df.list <- lapply(seq_along(WT.rasters), function(i) {
    r <- WT.rasters[[i]]
    
    # 1. Convert to data frame
    df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
    
    # 2. CHECK: If the raster was all NA, df has 0 rows. 
    if (nrow(df) == 0) {
      cat("⚠️ Skipping", names(r), "- No valid data points found.\n")
      return(NULL)
    }
    
    # 3. Standardize the column name
    value_col <- names(r)
    names(df)[names(df) == value_col] <- "wt"
    
    # 4. PLACE IT HERE: Assign the species name directly from the raster name
    df$species <- names(r) 
    
    return(df)
  })
  
  # Remove the NULL entries from the list before combining
  WT.df.list <- Filter(Negate(is.null), WT.df.list)
  
  # Combine into one big data frame
  WT.df <- do.call(rbind, WT.df.list)
  
  # remove underscores if in the species name
  WT.df$species <- gsub("_", " ", WT.df$species)

  # Dynamically create breaks and labels
  breaks_dynamic <- c(-Inf, seq(lower_bound, upper_bound, by = bin_width), Inf)
  
  labels_dynamic <- c(
    paste0("<", lower_bound),
    paste(seq(lower_bound, upper_bound - bin_width, by = bin_width),
          seq(lower_bound + bin_width, upper_bound, by = bin_width) - 1,
          sep = "-"),
    paste0(">", upper_bound)
  )
  
  # Bin warming tolerance values
  WT.df <- WT.df %>%
    dplyr::mutate(
      wt_bin = cut(wt,
                   breaks = breaks_dynamic,
                   labels = labels_dynamic,
                   right = FALSE)
    )
  
  # Define distinct colors for each bin, including reds, oranges, yellows, greens, blues
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
  
  # Adjust color vector length to match bins
  n_bins <- length(labels_dynamic)
  if (n_bins <= length(base_colors)) {
    bin_colors_dynamic <- base_colors[1:n_bins]
  } else {
    # Interpolate middle colors if more bins than base colors
    middle_colors <- colorRampPalette(base_colors[2:(length(base_colors)-1)])(n_bins - 2)
    bin_colors_dynamic <- c(base_colors[1], middle_colors, base_colors[length(base_colors)])
  }
  names(bin_colors_dynamic) <- labels_dynamic
  
  # Plot
  p = ggplot(WT.df, aes(x = x, y = y, fill = wt_bin)) +
    geom_raster() +
    geom_sf(data = EXT, fill = NA, color = "black", size = 0.3, inherit.aes = FALSE) +
    scale_fill_manual(
      values = bin_colors_dynamic,
      na.value = "grey90",
      name = layer_label,
      drop = FALSE
    ) +
    # Move the global theme UP
    theme_classic() + 
    # Use coord_sf ONCE for both alignment and limits
    coord_sf(xlim = xlim, ylim = ylim, crs = 4326, expand = FALSE) +
    facet_wrap(~ species, ncol = facet_ncol) +
    labs(x = "Longitude", y = "Latitude") +
    # Specific theme overrides MUST come after theme_classic()
    theme(
      legend.position = "right",
      strip.text = element_text(face = "italic")
    )
  
  
  return(list(
    plot = p,
    rasters = WT.rasters
  ))
  
}
