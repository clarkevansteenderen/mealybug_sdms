
# read in tuning p.citri.tuningults

p.citri.tuning = read.csv("tuning/p.citri/model_tuning_enmeval_full.csv")

# Select the model settings (RM and FC) that optimised AICc (delta AICc == 0)
best_delta_aicc = p.citri.tuning %>% 
  dplyr::filter(delta.AICc == 0) ;best_delta_aicc
best_delta_aicc_output = best_delta_aicc %>% t(.)
write.table(best_delta_aicc_output, "tuning/p.citri/bestAICc_LQH2_model.txt", quote = FALSE)

# select the model that optimised AUC (highest AUC value)
best_auc = p.citri.tuning %>% 
  dplyr::filter(auc.val.avg == max(auc.val.avg)) ;best_auc
best_auc_test_output = best_auc %>% t(.)
write.table(best_auc_test_output, "tuning/p.citri/bestAUC_LQH5_model.txt", quote = FALSE)

# select the model that optimised CBI (highest CBI value)
best_cbi = p.citri.tuning %>% 
  dplyr::filter(cbi.val.avg == max(cbi.val.avg)) ;best_cbi
best_cbi_output = best_cbi %>% t(.)
write.table(best_cbi_output, "tuning/p.citri/bestCBI_Q2_model.txt", quote = FALSE)

# select the model that optimised the 10% omission rate (lowest or.10p value)
best_or.10p.avg = p.citri.tuning %>% 
  dplyr::filter(or.10p.avg == min(or.10p.avg)) ;best_or.10p.avg
best_or.10p.avg_output = best_or.10p.avg %>% t(.)
write.table(best_or.10p.avg_output, "tuning/p.citri/bestOR10_L2_model.txt", quote = FALSE)

# default model output
default_mod_p.citri.tuningults = p.citri.tuning %>% 
  dplyr::filter(tune.args == "fc.LQH_rm.1") 
default_mod_p.citri.tuningults = default_mod_p.citri.tuningults %>% t(.)
write.table(default_mod_p.citri.tuningults, "tuning/p.citri/best_model_default.txt", quote = FALSE)



####################
# PCA
####################

PRESPOINTS = readRDS(file = "data/PRESPOINTS.ALLSPP.rds")

d.aberiae.pca = PRESPOINTS$d.aberiae %>%
  dplyr::select(!c(species, latitude, longitude, source))
head(d.aberiae.pca)

new_names = c(paste0("bio_", c("1", "6", "8", "9", "12")))

# Set the new column names
colnames(d.aberiae.pca) = new_names
head(d.aberiae.pca)

pca_res = prcomp(d.aberiae.pca)
summary(pca_res)

loadings = pca_res$x
loadings

PCA = factoextra::fviz_pca_var(pca_res,
                               col.var = "contrib", # Color by contributions to the PC
                               gradient.cols = c("blue", "gold", "red"),
                               repel = TRUE) + theme_classic(); PCA

# Contributions of variables to PC1
a = factoextra::fviz_contrib(pca_res, choice = "var", axes = 1)
# Contributions of variables to PC2
b = factoextra::fviz_contrib(pca_res, choice = "var", axes = 2)
pca_contribs = gridExtra::grid.arrange(a, b, ncol=2, top='Contribution of the variables to the first two PCs')
