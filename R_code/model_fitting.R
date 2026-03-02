
# Fit MaxEnt models, and get response plots

if (!dir.exists("models")) dir.create("models")

PRESPOINTS = readRDS(file = "data/PRESPOINTS.ALLSPP.rds")
BACKPOINTS = readRDS(file = "data/BACKPOINTS.ALLSPP.rds")

sapply(PRESPOINTS, nrow)

#############################################
# D. aberiae
#############################################

d.aberiae.preds = terra::rast("data/d.aberiae_current.tif")

d.aberiae.PRES = PRESPOINTS$d.aberiae %>%
  dplyr::select(-c(species, latitude, longitude, source)) %>%
  janitor::clean_names()

d.aberiae.BACK = BACKPOINTS$d.aberiae %>%
  dplyr::select(-c(latitude, longitude)) %>%
  janitor::clean_names()

# bind the presence and absence data together into one data frame
d.aberiae.data = dplyr::bind_rows(d.aberiae.PRES, d.aberiae.BACK)

# since bio 6 and bio 8 didn't contribute (response plots flatline), remove 
d.aberiae.data = d.aberiae.data %>%
  dplyr::select(!c(bio_6,
                   bio_8))
  
rownames(d.aberiae.data) = NULL

# Create a vector containing 0 (indicating background points) and 1 (indicating presence points)
d.aberiae.vector = c(
  replicate(nrow(d.aberiae.PRES), "1"),
  replicate(nrow(d.aberiae.BACK), "0")
) 

length(d.aberiae.vector)
names(d.aberiae.data)

# FIT MAXENT MODEL
d.aberiae.model = dismo::maxent(
  x = d.aberiae.data,
  p = d.aberiae.vector,
  path = here::here("models/d.aberiae_maxent_LQH1/"),
  replicates = 10,
  args = c(
    # Insert the optimal RM value here
    'betamultiplier=1.0',
    # Turn these on/off to change FC combinations 
    # - To only use quadratic features, turn all to false except quadratic
    'linear=true',
    'quadratic=true',
    'product=false',
    'threshold=false',
    'hinge=true',
    # Don't change anything from here down 
    'threads=2',
    'doclamp=true',
    'fadebyclamping=true',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false',
    'responsecurves=true',
    'writemess=true',
    'writeplotdata=true',
    'writebackgroundpredictions=true'
  )
)

dismo::response(d.aberiae.model)

bio.response.1 = as.data.frame(dismo::response(d.aberiae.model, 
              var = "bio_1"))  %>%
  dplyr::mutate(name = "bio1") %>% dplyr::rename(value = V1)

bio.response.9 = as.data.frame(dismo::response(d.aberiae.model, 
                                                var = "bio_9"))  %>%
  dplyr::mutate(name = "bio9") %>% dplyr::rename(value = V1)

bio.response.12 = as.data.frame(dismo::response(d.aberiae.model, 
                                                var = "bio_12"))  %>%
  dplyr::mutate(name = "bio12") %>% dplyr::rename(value = V1)

# join everything together
bioresponses = rbind(bio.response.1, bio.response.9, bio.response.12) %>%
  dplyr::mutate(name = factor(name, levels = c("bio1", "bio9", "bio12")))

# make a custom ggplot
response.plots.d.aberiae = ggplot(bioresponses, aes(x = value, y = p)) +
  geom_line(linewidth = 1, colour = "forestgreen") +
  facet_wrap(~name, scales = "free_x") +  # allow only x-axis to vary
  coord_cartesian(ylim = c(0, 1)) +       # constrain y-axis across facets
  labs(
    x = "Predictor Value",
    y = "Predicted Response", 
    title = expression("a) " * italic("Delottococcus aberiae"))
  )

response.plots.d.aberiae

#############################################
# N. viridis
#############################################

n.viridis.preds = terra::rast("data/n.viridis_current.tif")
#terra::plot(n.viridis.preds)

n.viridis.PRES = PRESPOINTS$n.viridis %>%
  dplyr::select(-c(species, latitude, longitude, source)) %>%
  janitor::clean_names()

n.viridis.BACK = BACKPOINTS$n.viridis %>%
  dplyr::select(-c(latitude, longitude)) %>%
  janitor::clean_names()

# bind the presence and absence data together into one data frame
n.viridis.data = dplyr::bind_rows(n.viridis.PRES, n.viridis.BACK)

rownames(n.viridis.data) = NULL

# Create a vector containing 0 (indicating background points) and 1 (indicating presence points)
n.viridis.vector = c(
  replicate(nrow(n.viridis.PRES), "1"),
  replicate(nrow(n.viridis.BACK), "0")
) 

length(n.viridis.vector)

n.viridis.data = n.viridis.data %>%
  dplyr::select(!c(bio_8))

# FIT MAXENT MODEL
n.viridis.model = dismo::maxent(
  x = n.viridis.data,
  p = n.viridis.vector,
  path = here::here("models/n.viridis_maxent_H3/"),
  replicates = 10,
  args = c(
    # Insert the optimal RM value here
    'betamultiplier=3.0',
    # Turn these on/off to change FC combinations 
    # - To only use quadratic features, turn all to false except quadratic
    'linear=false',
    'quadratic=false',
    'product=false',
    'threshold=false',
    'hinge=true',
    # Don't change anything from here down 
    'threads=2',
    'doclamp=true',
    'fadebyclamping=true',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false',
    'responsecurves=true',
    'writemess=true',
    'writeplotdata=true',
    'writebackgroundpredictions=true'
  )
)

dismo::response(n.viridis.model)

bio.response.1 = as.data.frame(dismo::response(n.viridis.model, 
                                               var = "bio_1"))  %>%
  dplyr::mutate(name = "bio1") %>% dplyr::rename(value = V1)

bio.response.5 = as.data.frame(dismo::response(n.viridis.model, 
                                               var = "bio_5"))  %>%
  dplyr::mutate(name = "bio5") %>% dplyr::rename(value = V1)

bio.response.12 = as.data.frame(dismo::response(n.viridis.model, 
                                                var = "bio_12"))  %>%
  dplyr::mutate(name = "bio12") %>% dplyr::rename(value = V1)

bio.response.19 = as.data.frame(dismo::response(n.viridis.model, 
                                                var = "bio_19"))  %>%
  dplyr::mutate(name = "bio19") %>% dplyr::rename(value = V1)

# join everything together
bioresponses = rbind(bio.response.1, bio.response.5, bio.response.12,
                     bio.response.19) %>%
  dplyr::mutate(name = factor(name, levels = c("bio1", "bio5", "bio12", "bio19")))

# make a custom ggplot
response.plots.n.viridis = ggplot(bioresponses, aes(x = value, y = p)) +
  geom_line(linewidth = 1, colour = "forestgreen") +
  facet_wrap(~name, scales = "free_x") +  # allow only x-axis to vary
  coord_cartesian(ylim = c(0, 1)) +       # constrain y-axis across facets
  labs(
    x = "Predictor Value",
    y = "Predicted Response", 
    title = expression("a) " * italic("Nipaecoccus viridis"))
  )

response.plots.n.viridis

#############################################
# P. burnerae
#############################################

# could not tune a model on this species due to low sample size
# using default model here, LQH1

p.burnerae.preds = terra::rast("data/p.burnerae_current.tif")

p.burnerae.PRES = PRESPOINTS$p.burnerae %>%
  dplyr::select(-c(species, latitude, longitude, source)) %>%
  janitor::clean_names()

p.burnerae.BACK = BACKPOINTS$p.burnerae %>%
  dplyr::select(-c(latitude, longitude)) %>%
  janitor::clean_names()

# bind the presence and absence data together into one data frame
p.burnerae.data = dplyr::bind_rows(p.burnerae.PRES, p.burnerae.BACK)

rownames(p.burnerae.data) = NULL

# Create a vector containing 0 (indicating background points) and 1 (indicating presence points)
p.burnerae.vector = c(
  replicate(nrow(p.burnerae.PRES), "1"),
  replicate(nrow(p.burnerae.BACK), "0")
) 

length(p.burnerae.vector)
names(p.burnerae.data)

p.burnerae.data = dplyr::select(p.burnerae.data,
                                bio_19)

# FIT MAXENT MODEL
p.burnerae.model = dismo::maxent(
  x = p.burnerae.data,
  p = p.burnerae.vector,
  path = here::here("models/p.burnerae_maxent_LQH1/"),
  replicates = 10,
  args = c(
    # Insert the optimal RM value here
    'betamultiplier=1.0',
    # Turn these on/off to change FC combinations 
    # - To only use quadratic features, turn all to false except quadratic
    'linear=true',
    'quadratic=true',
    'product=false',
    'threshold=false',
    'hinge=true',
    # Don't change anything from here down 
    'threads=2',
    'doclamp=true',
    'fadebyclamping=true',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false',
    'responsecurves=true',
    'writemess=true',
    'writeplotdata=true',
    'writebackgroundpredictions=true'
  )
)

dismo::response(p.burnerae.model)

bio.response.19 = as.data.frame(dismo::response(p.burnerae.model, 
                                                var = "bio_19"))  %>%
  dplyr::mutate(name = "bio19") %>% dplyr::rename(value = V1)

# make a custom ggplot
response.plots.p.buernae = ggplot(bio.response.19, aes(x = value, y = p)) +
  geom_line(linewidth = 1, colour = "forestgreen") +
  facet_wrap(~name, scales = "free_x") +  # allow only x-axis to vary
  coord_cartesian(ylim = c(0, 1)) +       # constrain y-axis across facets
  labs(
    x = "Predictor Value",
    y = "Predicted Response", 
    title = expression("c) " * italic("Paracoccus burnerae"))
  )

response.plots.p.buernae

#############################################
# P. calceolariae
#############################################

p.calceolariae.preds = terra::rast("data/p.calceolariae_current.tif")

p.calceolariae.PRES = PRESPOINTS$p.calceolariae %>%
  dplyr::select(-c(species, latitude, longitude, source)) %>%
  janitor::clean_names()

p.calceolariae.BACK = BACKPOINTS$p.calceolariae %>%
  dplyr::select(-c(latitude, longitude)) %>%
  janitor::clean_names()

# bind the presence and absence data together into one data frame
p.calceolariae.data = dplyr::bind_rows(p.calceolariae.PRES, p.calceolariae.BACK)

rownames(p.calceolariae.data) = NULL

# Create a vector containing 0 (indicating background points) and 1 (indicating presence points)
p.calceolariae.vector = c(
  replicate(nrow(p.calceolariae.PRES), "1"),
  replicate(nrow(p.calceolariae.BACK), "0")
) 

length(p.calceolariae.vector)

# FIT MAXENT MODEL
p.calceolariae.model = dismo::maxent(
  x = p.calceolariae.data,
  p = p.calceolariae.vector,
  path = here::here("models/p.calceolariae_maxent_LQH5/"),
  replicates = 10,
  args = c(
    # Insert the optimal RM value here
    'betamultiplier=5.0',
    # Turn these on/off to change FC combinations 
    # - To only use quadratic features, turn all to false except quadratic
    'linear=true',
    'quadratic=true',
    'product=false',
    'threshold=false',
    'hinge=true',
    # Don't change anything from here down 
    'threads=2',
    'doclamp=true',
    'fadebyclamping=true',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false',
    'responsecurves=true',
    'writemess=true',
    'writeplotdata=true',
    'writebackgroundpredictions=true'
  )
)

dismo::response(p.calceolariae.model)

bio.response.1 = as.data.frame(dismo::response(p.calceolariae.model, 
                                               var = "bio_1"))  %>%
  dplyr::mutate(name = "bio1") %>% dplyr::rename(value = V1)

bio.response.8 = as.data.frame(dismo::response(p.calceolariae.model, 
                                               var = "bio_8"))  %>%
  dplyr::mutate(name = "bio8") %>% dplyr::rename(value = V1)

bio.response.9 = as.data.frame(dismo::response(p.calceolariae.model, 
                                                var = "bio_9"))  %>%
  dplyr::mutate(name = "bio9") %>% dplyr::rename(value = V1)

bio.response.12 = as.data.frame(dismo::response(p.calceolariae.model, 
                                                var = "bio_12"))  %>%
  dplyr::mutate(name = "bio12") %>% dplyr::rename(value = V1)


# join everything together
bioresponses = rbind(bio.response.1, bio.response.8, bio.response.9, bio.response.12) %>%
  dplyr::mutate(name = factor(name, levels = c("bio1", "bio8", "bio9", "bio12")))

# make a custom ggplot
response.plots.p.calceolariae = ggplot(bioresponses, aes(x = value, y = p)) +
  geom_line(linewidth = 1, colour = "forestgreen") +
  facet_wrap(~name, scales = "free_x") +  # allow only x-axis to vary
  coord_cartesian(ylim = c(0, 1)) +       # constrain y-axis across facets
  labs(
    x = "Predictor Value",
    y = "Predicted Response", 
    title = expression("c) " * italic("Pseudococcus calceolariae"))
  )

response.plots.p.calceolariae

#############################################
# P. citri
#############################################

p.citir.preds = terra::rast("data/p.citri_current.tif")

p.citri.PRES = PRESPOINTS$p.citri %>%
  dplyr::select(-c(species, latitude, longitude, source)) %>%
  janitor::clean_names()

p.citri.BACK = BACKPOINTS$p.citri %>%
  dplyr::select(-c(latitude, longitude)) %>%
  janitor::clean_names()

# bind the presence and absence data together into one data frame
p.citri.data = dplyr::bind_rows(p.citri.PRES, p.citri.BACK)

rownames(p.citri.data) = NULL

# Create a vector containing 0 (indicating background points) and 1 (indicating presence points)
p.citri.vector = c(
  replicate(nrow(p.citri.PRES), "1"),
  replicate(nrow(p.citri.BACK), "0")
) 

length(p.citri.vector)
names(p.citri.data)

p.citri.data = dplyr::select(p.citri.data, 
                             !c(bio_3, bio_8))

# FIT MAXENT MODEL
p.citri.model = dismo::maxent(
  x = p.citri.data,
  p = p.citri.vector,
  path = here::here("models/p.citri_maxent_LQH2/"),
  replicates = 10,
  args = c(
    # Insert the optimal RM value here
    'betamultiplier=2.0',
    # Turn these on/off to change FC combinations 
    # - To only use quadratic features, turn all to false except quadratic
    'linear=true',
    'quadratic=true',
    'product=false',
    'threshold=false',
    'hinge=true',
    # Don't change anything from here down 
    'threads=2',
    'doclamp=true',
    'fadebyclamping=true',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false',
    'responsecurves=true',
    'writemess=true',
    'writeplotdata=true',
    'writebackgroundpredictions=true'
  )
)

dismo::response(p.citri.model)

bio.response.1 = as.data.frame(dismo::response(p.citri.model, 
                                               var = "bio_1"))  %>%
  dplyr::mutate(name = "bio1") %>% dplyr::rename(value = V1)

bio.response.2 = as.data.frame(dismo::response(p.citri.model, 
                                               var = "bio_2"))  %>%
  dplyr::mutate(name = "bio2") %>% dplyr::rename(value = V1)

bio.response.3 = as.data.frame(dismo::response(p.citri.model, 
                                               var = "bio_3"))  %>%
  dplyr::mutate(name = "bio3") %>% dplyr::rename(value = V1)

bio.response.5 = as.data.frame(dismo::response(p.citri.model, 
                                               var = "bio_5"))  %>%
  dplyr::mutate(name = "bio5") %>% dplyr::rename(value = V1)

bio.response.9 = as.data.frame(dismo::response(p.citri.model, 
                                               var = "bio_9"))  %>%
  dplyr::mutate(name = "bio9") %>% dplyr::rename(value = V1)

bio.response.12 = as.data.frame(dismo::response(p.citri.model, 
                                                var = "bio_12"))  %>%
  dplyr::mutate(name = "bio12") %>% dplyr::rename(value = V1)

bio.response.17 = as.data.frame(dismo::response(p.citri.model, 
                                                var = "bio_17"))  %>%
  dplyr::mutate(name = "bio17") %>% dplyr::rename(value = V1)

bio.response.19 = as.data.frame(dismo::response(p.citri.model, 
                                                var = "bio_19"))  %>%
  dplyr::mutate(name = "bio19") %>% dplyr::rename(value = V1)


# join everything together
bioresponses = rbind(bio.response.1, bio.response.2,
                     bio.response.5, bio.response.9, bio.response.12, bio.response.17,
                     bio.response.19) %>%
  dplyr::mutate(name = factor(name, levels = c("bio1", "bio2", "bio3",
                                              "bio5", "bio9",  "bio12", "bio17", "bio19")))

# make a custom ggplot
response.plots.p.citri = ggplot(bioresponses, aes(x = value, y = p)) +
  geom_line(linewidth = 1, colour = "forestgreen") +
  facet_wrap(~name, scales = "free_x") +  # allow only x-axis to vary
  coord_cartesian(ylim = c(0, 1)) +       # constrain y-axis across facets
  labs(
    x = "Predictor Value",
    y = "Predicted Response", 
    title = expression("b) " * italic("Planococcus citri"))
  )

response.plots.p.citri

#############################################
# P. longispinus
#############################################

p.longispinus.preds = terra::rast("data/p.longispinus_current.tif")

p.longispinus.PRES = PRESPOINTS$p.longispinus %>%
  dplyr::select(-c(species, latitude, longitude, source)) %>%
  janitor::clean_names()

p.longispinus.BACK = BACKPOINTS$p.longispinus %>%
  dplyr::select(-c(latitude, longitude)) %>%
  janitor::clean_names()

# bind the presence and absence data together into one data frame
p.longispinus.data = dplyr::bind_rows(p.longispinus.PRES, p.longispinus.BACK)

rownames(p.longispinus.data) = NULL

# Create a vector containing 0 (indicating background points) and 1 (indicating presence points)
p.longispinus.vector = c(
  replicate(nrow(p.longispinus.PRES), "1"),
  replicate(nrow(p.longispinus.BACK), "0")
) 

length(p.longispinus.vector)
names(p.longispinus.data)

p.longispinus.data = 
  dplyr::select(p.longispinus.data, !(bio_19))

# FIT MAXENT MODEL
p.longispinus.model = dismo::maxent(
  x = p.longispinus.data,
  p = p.longispinus.vector,
  path = here::here("models/p.longispinus_maxent_LQH2/"),
  replicates = 10,
  args = c(
    # Insert the optimal RM value here
    'betamultiplier=2.0',
    # Turn these on/off to change FC combinations 
    # - To only use quadratic features, turn all to false except quadratic
    'linear=true',
    'quadratic=true',
    'product=false',
    'threshold=false',
    'hinge=true',
    # Don't change anything from here down 
    'threads=2',
    'doclamp=true',
    'fadebyclamping=true',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false',
    'responsecurves=true',
    'writemess=true',
    'writeplotdata=true',
    'writebackgroundpredictions=true'
  )
)

dismo::response(p.longispinus.model)

bio.response.1 = as.data.frame(dismo::response(p.longispinus.model, 
                                               var = "bio_1"))  %>%
  dplyr::mutate(name = "bio1") %>% dplyr::rename(value = V1)

bio.response.7 = as.data.frame(dismo::response(p.longispinus.model, 
                                               var = "bio_7"))  %>%
  dplyr::mutate(name = "bio7") %>% dplyr::rename(value = V1)

bio.response.8 = as.data.frame(dismo::response(p.longispinus.model, 
                                               var = "bio_8"))  %>%
  dplyr::mutate(name = "bio8") %>% dplyr::rename(value = V1)

bio.response.10 = as.data.frame(dismo::response(p.longispinus.model, 
                                               var = "bio_10"))  %>%
  dplyr::mutate(name = "bio10") %>% dplyr::rename(value = V1)

bio.response.12 = as.data.frame(dismo::response(p.longispinus.model, 
                                               var = "bio_12"))  %>%
  dplyr::mutate(name = "bio12") %>% dplyr::rename(value = V1)


# join everything together
bioresponses = rbind(bio.response.1, bio.response.7, bio.response.8, 
                     bio.response.10, bio.response.12) %>%
  dplyr::mutate(name = factor(name, levels = c("bio1", "bio7", "bio8",
                                               "bio10", "bio12")))

# make a custom ggplot
response.plots.p.longispinus = ggplot(bioresponses, aes(x = value, y = p)) +
  geom_line(linewidth = 1, colour = "forestgreen") +
  facet_wrap(~name, scales = "free_x") +  # allow only x-axis to vary
  coord_cartesian(ylim = c(0, 1)) +       # constrain y-axis across facets
  labs(
    x = "Predictor Value",
    y = "Predicted Response", 
    title = expression("d) " * italic("Pseudococcus longispinus"))
  )

response.plots.p.longispinus

response.plots.all = gridExtra::grid.arrange(response.plots.n.viridis,
                        response.plots.p.citri,
                        response.plots.p.calceolariae, response.plots.p.longispinus)

ggsave("figures/response_plots/response_plots_all.png", plot = response.plots.all, 
       width = 10, height = 8, dpi = 450)
