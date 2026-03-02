# Load libraries
library(tidyverse)
library(viridis)
library(tidytext)  # for reorder_within

# Original data
data <- tribble(
  ~species,            ~bio,     ~value,
  "Nipaecoccus viridis",         "bio12",  46.7,
  "Nipaecoccus viridis",         "bio19",  21,
  "Nipaecoccus viridis",         "bio1",   18.4,
  "Nipaecoccus viridis",         "bio5",   13.8,
  
  "Pseudococcus calceolariae",   "bio1",   69.5,
  "Pseudococcus calceolariae",   "bio8",   19.8,
  "Pseudococcus calceolariae",   "bio9",   9.2,
  "Pseudococcus calceolariae",   "bio12",  1.5,
  
  "Planococcus citri",          "bio19",  36.8,
  "Planococcus citri",          "bio1",   21.8,
  "Planococcus citri",          "bio2",   19.3,
  "Planococcus citri",          "bio5",   16.2,
  "Planococcus citri",          "bio17",  3.4,
  "Planococcus citri",          "bio9",   2,
  "Planococcus citri",          "bio12",  0.6,
  
  "Pseudococcus longispinus",    "bio1",   45,
  "Pseudococcus longispinus",    "bio7",   37.8,
  "Pseudococcus longispinus",    "bio12",  8.5,
  "Pseudococcus longispinus",    "bio10",  4.9,
  "Pseudococcus longispinus",    "bio8",   3.7
)

# Plot: reorder bars within each species
ggplot(data, aes(x = reorder_within(bio, value, species), 
                 y = value, fill = species)) +
  geom_col(show.legend = FALSE, colour = "black") +
  scale_x_reordered() +
  scale_fill_viridis_d(option = "plasma", end = 0.9) +
  facet_wrap(~species, scales = "free_x") +  # each facet independent
  labs(
    #title = "Bioclimatic variable importance by species",
    x = "Bioclimatic variable",
    y = "Percentage contribution"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text = element_text(face = "italic")  # italicize facet labels
  )

