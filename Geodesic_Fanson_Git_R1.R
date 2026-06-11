############################################################
###### FANSON ET AL. VALIDATION #############################
###### Geodesic vs Euclidean costs ##########################
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(terra)
library(patchwork)

# --- Load data and nutrigonometry functions ---
source("/...General_functions_Nutrigonometry_03_December_2025_Geodesics.R")
source("...Geodesic_consolidated_Functions_6_December_2025.R")
set.seed(321909)
############################################################
###### FANSON DATA CLEANING - INDIVIDUAL LEVEL #############
############################################################
fansondt_raw <- read.csv(
  "...Fansonetal2008AgeingCell-NoChoiceData.csv",
  strip.white = TRUE,
  stringsAsFactors = FALSE,
  header = TRUE
)

# Adjusting Ratio levels
fansondt_raw$Ratio <- ifelse(fansondt_raw$Ratio == "0:01", "0:1",
                             ifelse(fansondt_raw$Ratio == "1:00", "1:0", 
                                    ifelse(fansondt_raw$Ratio == "1:01", "1:1",
                                           ifelse(fansondt_raw$Ratio == "1:02", "1:2",
                                                  ifelse(fansondt_raw$Ratio == "1:04", "1:4",
                                                         ifelse(fansondt_raw$Ratio == "1:08", "1:8", "1:16"))))))

# Adjusting the Food column
fansondt_raw <- tidyr::separate(
  fansondt_raw,
  Food,
  into = c("Food", "unit"),
  sep = " ",
  remove = TRUE
)

fansondt_raw$Food <- as.integer(fansondt_raw$Food)

# Calculating carbohydrate intake
fansondt_raw$Ceaten <- ifelse(fansondt_raw$Ratio == "0:1", fansondt_raw$totaleaten * fansondt_raw$Food,
                              ifelse(fansondt_raw$Ratio == "1:0", 0,
                                     ifelse(fansondt_raw$Ratio == "1:1", (fansondt_raw$totaleaten / 2) * fansondt_raw$Food,
                                            ifelse(fansondt_raw$Ratio == "1:2", (fansondt_raw$totaleaten / 3) * 2 * fansondt_raw$Food,
                                                   ifelse(fansondt_raw$Ratio == "1:4", (fansondt_raw$totaleaten / 5) * 4 * fansondt_raw$Food,
                                                          ifelse(fansondt_raw$Ratio == "1:8", (fansondt_raw$totaleaten / 9) * 8 * fansondt_raw$Food,
                                                                 (fansondt_raw$totaleaten / 17) * 16 * fansondt_raw$Food))))))

# Calculating protein intake
fansondt_raw$Peaten <- (fansondt_raw$totaleaten * fansondt_raw$Food) - fansondt_raw$Ceaten

# Getting mg of macronutrients
fansondt_raw$Ceaten <- fansondt_raw$Ceaten / 1000
fansondt_raw$Peaten <- fansondt_raw$Peaten / 1000

# Cleaning the data
fansondt_raw <- subset(fansondt_raw, Ceaten >= 0 & Peaten >= 0)


# Collapse to individual/fly level
fansondt <- fansondt_raw %>%
  group_by(trial, DietID, replicate, flyID, Ratio, Food) %>%
  summarise(
    lifespan = max(lifespan, na.rm = TRUE),
    eggs = sum(eggs, na.rm = TRUE),
    Peaten = sum(Peaten, na.rm = TRUE),
    Ceaten = sum(Ceaten, na.rm = TRUE),
    totaleaten = sum(totaleaten, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    dailyeggs = eggs / lifespan,
    lifetimeegg = eggs,
    Ratio = as.factor(Ratio),
    P_scaled = scale(Peaten, scale = TRUE)[, 1],
    C_scaled = scale(Ceaten, scale = TRUE)[, 1],
    #lifespan_std = as.numeric(scale(lifespan)),
    #lifetimeegg_std = as.numeric(scale(lifetimeegg)),
    #dailyeggs_std = as.numeric(scale(dailyeggs))
    lifespan_std = as.numeric(lifespan/mean(lifespan, na.rm = TRUE)),
    lifetimeegg_std = as.numeric(lifetimeegg/mean(lifetimeegg, na.rm = TRUE)),
    dailyeggs_std = as.numeric(dailyeggs/mean(dailyeggs, na.rm = TRUE))
  ) %>%
  filter(
    is.finite(P_scaled),
    is.finite(C_scaled),
    is.finite(lifespan),
    is.finite(lifetimeegg),
    is.finite(dailyeggs)
  )
############################################################
###### 2. Fit quadratic performance landscapes #############
############################################################

fit_fanson_lifespan <- fit_quadratic_surface(
  fansondt,
  "P_scaled",
  "C_scaled",
  "lifespan_std"
)

fit_fanson_lifetimeegg <- fit_quadratic_surface(
  fansondt,
  "P_scaled",
  "C_scaled",
  "lifetimeegg_std"
)

fit_fanson_dailyeggs <- fit_quadratic_surface(
  fansondt,
  "P_scaled",
  "C_scaled",
  "dailyeggs_std"
)

############################################################
###### 3. Identify peak regions ############################
############################################################

prediction_fanson_lifespan <- predict_region_GF_lm(
  fansondt,
  Pcol = "P_scaled",
  Ccol = "C_scaled",
  Zcol = "lifespan_std",
  peak = "yes"
)

prediction_fanson_lifetimeegg <- predict_region_GF_lm(
  fansondt,
  Pcol = "P_scaled",
  Ccol = "C_scaled",
  Zcol = "lifetimeegg_std",
  peak = "yes"
)

prediction_fanson_dailyeggs <- predict_region_GF_lm(
  fansondt,
  Pcol = "P_scaled",
  Ccol = "C_scaled",
  Zcol = "dailyeggs_std",
  peak = "yes"
)

############################################################
###### 4. Peak-region hulls and centroids ##################
############################################################

# Lifespan
hull_fanson_lifespan <- design_hull(prediction_fanson_lifespan)

Fanson_lifespan_centroid <- as.data.frame(
  terra::crds(
    terra::centroids(
      convHull(
        vect(cbind(
          x = hull_fanson_lifespan$PROT,
          y = hull_fanson_lifespan$CARB
        ))
      )
    )
  )
)
colnames(Fanson_lifespan_centroid) <- c("PROT", "CARB")

# Lifetime eggs
hull_fanson_lifetimeegg <- design_hull(prediction_fanson_lifetimeegg)

Fanson_lifetimeegg_centroid <- as.data.frame(
  terra::crds(
    terra::centroids(
      convHull(
        vect(cbind(
          x = hull_fanson_lifetimeegg$PROT,
          y = hull_fanson_lifetimeegg$CARB
        ))
      )
    )
  )
)
colnames(Fanson_lifetimeegg_centroid) <- c("PROT", "CARB")

# Repr. rate
hull_fanson_dailyeggs <- design_hull(prediction_fanson_dailyeggs)

Fanson_dailyeggs_centroid <- as.data.frame(
  terra::crds(
    terra::centroids(
      convHull(
        vect(cbind(
          x = hull_fanson_dailyeggs$PROT,
          y = hull_fanson_dailyeggs$CARB
        ))
      )
    )
  )
)
colnames(Fanson_dailyeggs_centroid) <- c("PROT", "CARB")



hull_fanson_lifespan
hull_fanson_dailyeggs
hull_fanson_lifetimeegg


### testing for nutritional trade-off lifetimeegg vs repr. rate
mean((hull_fanson_lifetimeegg$CARB)/(hull_fanson_lifetimeegg$PROT))
mean((hull_fanson_dailyeggs$CARB)/(hull_fanson_dailyeggs$PROT))

wilcox.test((hull_fanson_lifetimeegg$CARB)/(hull_fanson_lifetimeegg$PROT),
            (hull_fanson_dailyeggs$CARB)/(hull_fanson_dailyeggs$PROT))
############################################################
###### 5. Predict surfaces #################################
############################################################

Fanson_lifespan_surface <- predict_surface(
  fansondt,
  x1 = "P_scaled",
  x2 = "C_scaled",
  y = "lifespan_std"
)

Fanson_lifetimeegg_surface <- predict_surface(
  fansondt,
  x1 = "P_scaled",
  x2 = "C_scaled",
  y = "lifetimeegg_std",
  lambda = 0.15
)

Fanson_dailyeggs_surface <- predict_surface(
  fansondt,
  x1 = "P_scaled",
  x2 = "C_scaled",
  y = "dailyeggs_std",
  lambda = 0.15
)

############################################################
###### 6. Basic landscape plots ############################
############################################################

Fanson_lifespan_plot <- ggplot(
  data = Fanson_lifespan_surface,
  aes(x = PROT, y = CARB, z = tps, fill = tps)
) +
  geom_tile(alpha = 0.9) +
  geom_contour(
    data = Fanson_lifespan_surface,
    aes(x = PROT, y = CARB, z = tps),
    color = "black",
    size = 0.08
  ) +
  geom_polygon(
    hull_fanson_lifespan,
    mapping = aes(x = PROT, y = CARB),
    fill = "black",
    col = "black",
    inherit.aes = FALSE,
    alpha = 0.1,
    size = 0.5
  ) +
  geom_point(
    data = Fanson_lifespan_centroid,
    mapping = aes(x = PROT, y = CARB),
    pch = 21,
    fill = "steelblue2",
    col = "black",
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_fill_distiller(
    palette = "Spectral",
    name = "Lifespan"
  ) +
  xlab("Protein intake (mg)") +
  ylab("Carbohydrate intake (mg)") +
  ggtitle("Lifespan") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 1.2,
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 10)
  )

Fanson_lifetimeegg_plot <- ggplot(
  data = Fanson_lifetimeegg_surface,
  aes(x = PROT, y = CARB, z = tps, fill = tps)
) +
  geom_tile(alpha = 0.9) +
  geom_contour(
    data = Fanson_lifetimeegg_surface,
    aes(x = PROT, y = CARB, z = tps),
    color = "black",
    size = 0.08
  ) +
  geom_polygon(
    hull_fanson_lifetimeegg,
    mapping = aes(x = PROT, y = CARB),
    fill = "black",
    col = "black",
    inherit.aes = FALSE,
    alpha = 0.1,
    size = 0.5
  ) +
  geom_point(
    data = Fanson_lifetimeegg_centroid,
    mapping = aes(x = PROT, y = CARB),
    pch = 21,
    fill = "orange",
    col = "black",
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_fill_distiller(
    palette = "Spectral",
    name = "Lifetime eggs"
  ) +
  xlab("Protein intake (mg)") +
  ylab("Carbohydrate intake (mg)") +
  ggtitle("Lifetime eggs") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 1.2,
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 10)
  )

Fanson_dailyeggs_plot <- ggplot(
  data = Fanson_dailyeggs_surface,
  aes(x = PROT, y = CARB, z = tps, fill = tps)
) +
  geom_tile(alpha = 0.9) +
  geom_contour(
    data = Fanson_dailyeggs_surface,
    aes(x = PROT, y = CARB, z = tps),
    color = "black",
    size = 0.08
  ) +
  geom_polygon(
    hull_fanson_dailyeggs,
    mapping = aes(x = PROT, y = CARB),
    fill = "black",
    col = "black",
    inherit.aes = FALSE,
    alpha = 0.1,
    size = 0.5
  ) +
  geom_point(
    data = Fanson_dailyeggs_centroid,
    mapping = aes(x = PROT, y = CARB),
    pch = 21,
    fill = "firebrick2",
    col = "black",
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_fill_distiller(
    palette = "Spectral",
    name = "Repr. rate"
  ) +
  xlab("Protein intake (mg)") +
  ylab("Carbohydrate intake (mg)") +
  ggtitle("Repr. rate") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 1.2,
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 10)
  )

Fanson_lifespan_plot + Fanson_lifetimeegg_plot + Fanson_dailyeggs_plot

############################################################
###### 7. Generate geodesic paths ##########################
############################################################

set.seed(1235)

N <- 25

lifespan_peak_small_fanson <- prediction_fanson_lifespan$predicted_df %>% sample_n(N)
lifetime_peak_small_fanson <- prediction_fanson_lifetimeegg$predicted_df %>% sample_n(N)
dailyeggs_peak_small_fanson <- prediction_fanson_dailyeggs$predicted_df %>% sample_n(N)

############################################################
###### Lifespan vs Repr. rate on lifespan surface ##########
############################################################

a <- fit_fanson_lifespan$coefs[["a"]]
b <- fit_fanson_lifespan$coefs[["b"]]
c <- fit_fanson_lifespan$coefs[["c"]]
d <- fit_fanson_lifespan$coefs[["d"]]
e <- fit_fanson_lifespan$coefs[["e"]]

all_paths_fanson_lifespan <- list()

for (i in 1:N) {
  
  x0 <- lifespan_peak_small_fanson$PROT[i]
  y0 <- lifespan_peak_small_fanson$CARB[i]
  x1 <- dailyeggs_peak_small_fanson$PROT[i]
  y1 <- dailyeggs_peak_small_fanson$CARB[i]
  
  cat(sprintf("Fanson lifespan surface: sample %d of %d...\n", i, N))
  
  path <- tryCatch({
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    df <- data.frame(integrate_geodesic(vT, x0, y0))
    df$sample <- i
    all_paths_fanson_lifespan[[i]] <- df
  }, error = function(e) NULL)
}

paths_fanson_lifespan <- bind_rows(all_paths_fanson_lifespan) %>%
  mutate(
    trait = "Lifespan-Repr. rate",
    surface = "Lifespan",
    type = "geodesic"
  )

############################################################
###### Lifespan vs Repr. rate on Repr. rate surface ########
############################################################

a <- fit_fanson_dailyeggs$coefs[["a"]]
b <- fit_fanson_dailyeggs$coefs[["b"]]
c <- fit_fanson_dailyeggs$coefs[["c"]]
d <- fit_fanson_dailyeggs$coefs[["d"]]
e <- fit_fanson_dailyeggs$coefs[["e"]]

all_paths_fanson_dailyeggs <- list()

for (i in 1:N) {
  
  x0 <- lifespan_peak_small_fanson$PROT[i]
  y0 <- lifespan_peak_small_fanson$CARB[i]
  x1 <- dailyeggs_peak_small_fanson$PROT[i]
  y1 <- dailyeggs_peak_small_fanson$CARB[i]
  
  cat(sprintf("Fanson Repr. rate surface: sample %d of %d...\n", i, N))
  
  path <- tryCatch({
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    df <- data.frame(integrate_geodesic(vT, x0, y0))
    df$sample <- i
    all_paths_fanson_dailyeggs[[i]] <- df
  }, error = function(e) NULL)
}

paths_fanson_dailyeggs <- bind_rows(all_paths_fanson_dailyeggs) %>%
  mutate(
    trait = "Lifespan-Repr. rate",
    surface = "Repr. rate",
    type = "geodesic"
  )

############################################################
###### Lifetime eggs vs lifespan on lifetime surface #######
############################################################

a <- fit_fanson_lifetimeegg$coefs[["a"]]
b <- fit_fanson_lifetimeegg$coefs[["b"]]
c <- fit_fanson_lifetimeegg$coefs[["c"]]
d <- fit_fanson_lifetimeegg$coefs[["d"]]
e <- fit_fanson_lifetimeegg$coefs[["e"]]

all_paths_fanson_lifetime_lifespan <- list()

for (i in 1:N) {
  
  x0 <- lifetime_peak_small_fanson$PROT[i]
  y0 <- lifetime_peak_small_fanson$CARB[i]
  x1 <- lifespan_peak_small_fanson$PROT[i]
  y1 <- lifespan_peak_small_fanson$CARB[i]
  
  cat(sprintf("Fanson lifetime surface: lifespan comparison sample %d of %d...\n", i, N))
  
  path <- tryCatch({
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    df <- data.frame(integrate_geodesic(vT, x0, y0))
    df$sample <- i
    all_paths_fanson_lifetime_lifespan[[i]] <- df
  }, error = function(e) NULL)
}

paths_fanson_lifetime_lifespan <- bind_rows(all_paths_fanson_lifetime_lifespan) %>%
  mutate(
    trait = "Lifespan-lifetime eggs",
    surface = "Lifetime eggs",
    type = "geodesic"
  )

############################################################
###### Lifetime eggs vs Repr. rate on lifetime surface #####
############################################################

all_paths_fanson_lifetime_dailyeggs <- list()

for (i in 1:N) {
  
  x0 <- lifetime_peak_small_fanson$PROT[i]
  y0 <- lifetime_peak_small_fanson$CARB[i]
  x1 <- dailyeggs_peak_small_fanson$PROT[i]
  y1 <- dailyeggs_peak_small_fanson$CARB[i]
  
  cat(sprintf("Fanson lifetime surface: Repr. rate comparison sample %d of %d...\n", i, N))
  
  path <- tryCatch({
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    df <- data.frame(integrate_geodesic(vT, x0, y0))
    df$sample <- i
    all_paths_fanson_lifetime_dailyeggs[[i]] <- df
  }, error = function(e) NULL)
}

paths_fanson_lifetime_dailyeggs <- bind_rows(all_paths_fanson_lifetime_dailyeggs) %>%
  mutate(
    trait = "Lifetime eggs-Repr. rate",
    surface = "Lifetime eggs",
    type = "geodesic"
  )

############################################################
###### Lifetime eggs vs lifespan on lifespan surface #######
############################################################

a <- fit_fanson_lifespan$coefs[["a"]]
b <- fit_fanson_lifespan$coefs[["b"]]
c <- fit_fanson_lifespan$coefs[["c"]]
d <- fit_fanson_lifespan$coefs[["d"]]
e <- fit_fanson_lifespan$coefs[["e"]]

all_paths_fanson_lifetime_lifespan_surflifespan <- list()

for (i in 1:N) {
  
  x0 <- lifetime_peak_small_fanson$PROT[i]
  y0 <- lifetime_peak_small_fanson$CARB[i]
  x1 <- lifespan_peak_small_fanson$PROT[i]
  y1 <- lifespan_peak_small_fanson$CARB[i]
  
  cat(sprintf("Fanson lifespan surface: lifetime comparison sample %d of %d...\n", i, N))
  
  path <- tryCatch({
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    df <- data.frame(integrate_geodesic(vT, x0, y0))
    df$sample <- i
    all_paths_fanson_lifetime_lifespan_surflifespan[[i]] <- df
  }, error = function(e) NULL)
}

paths_fanson_lifetime_lifespan_surflifespan <- bind_rows(all_paths_fanson_lifetime_lifespan_surflifespan) %>%
  mutate(
    trait = "Lifespan-lifetime eggs",
    surface = "Lifespan",
    type = "geodesic"
  )

############################################################
###### Lifetime eggs vs Repr. rate on Repr. rate surface ###
############################################################

a <- fit_fanson_dailyeggs$coefs[["a"]]
b <- fit_fanson_dailyeggs$coefs[["b"]]
c <- fit_fanson_dailyeggs$coefs[["c"]]
d <- fit_fanson_dailyeggs$coefs[["d"]]
e <- fit_fanson_dailyeggs$coefs[["e"]]

all_paths_fanson_lifetime_dailyeggs_surfdailyeggs <- list()

for (i in 1:N) {
  
  x0 <- lifetime_peak_small_fanson$PROT[i]
  y0 <- lifetime_peak_small_fanson$CARB[i]
  x1 <- dailyeggs_peak_small_fanson$PROT[i]
  y1 <- dailyeggs_peak_small_fanson$CARB[i]
  
  cat(sprintf("Fanson Repr. rate surface: lifetime comparison sample %d of %d...\n", i, N))
  
  path <- tryCatch({
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    df <- data.frame(integrate_geodesic(vT, x0, y0))
    df$sample <- i
    all_paths_fanson_lifetime_dailyeggs_surfdailyeggs[[i]] <- df
  }, error = function(e) NULL)
}

paths_fanson_lifetime_dailyeggs_surfdailyeggs <- bind_rows(all_paths_fanson_lifetime_dailyeggs_surfdailyeggs) %>%
  mutate(
    trait = "Lifetime eggs-Repr. rate",
    surface = "Repr. rate",
    type = "geodesic"
  )

############################################################
###### Combine and save paths ##############################
############################################################

paths_fanson_all <- bind_rows(
  paths_fanson_lifespan,
  paths_fanson_dailyeggs,
  paths_fanson_lifetime_lifespan,
  paths_fanson_lifetime_dailyeggs,
  paths_fanson_lifetime_lifespan_surflifespan,
  paths_fanson_lifetime_dailyeggs_surfdailyeggs
)

write.csv(
  paths_fanson_all,
  "/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Nutritional Geometry/Nutrigonometry/23. Geodesics/Final Scripts/Final codes/Revisions 1/Fanson analysis/Fanson_et_al_Geodesic_paths_Lifespan_DailyEggs_Lifetime_N_50.csv",
  row.names = FALSE
)

############################################################
###### 8. Area / distance estimates ########################
############################################################

area_df_fanson_trait_surface <- c()

for (j in 1:N) {
  
  cat(sprintf("Fanson area sample %d of %d...\n", j, N))
  
  ##########################################################
  ## 1. Lifespan-Repr. rate on Lifespan surface ############
  ##########################################################
  
  path_df <- subset(
    paths_fanson_all,
    trait == "Lifespan-Repr. rate" &
      surface == "Lifespan" &
      sample == j
  )
  
  if (nrow(path_df) > 0) {
    
    inter_geo <- area_geodesic_with_peak(
      surface_df = Fanson_lifespan_surface,
      x_col_surf = "PROT",
      y_col_surf = "CARB",
      z_col_surf = "tps",
      path_df = path_df,
      standardize = TRUE,
      use_metric = TRUE,
      signed_area = FALSE
    )
    
    inter_euc <- area_geodesic_with_peak(
      surface_df = Fanson_lifespan_surface,
      x_col_surf = "PROT",
      y_col_surf = "CARB",
      z_col_surf = "tps",
      path_df = path_df,
      standardize = TRUE,
      use_metric = FALSE,
      signed_area = FALSE
    )
    
    area_df_fanson_trait_surface <- bind_rows(
      area_df_fanson_trait_surface,
      data.frame(
        area_below_geo = inter_geo$area_below_geo_baseline,
        area_between_geo = inter_geo$area_between_geo_peak,
        path_length_geo = inter_geo$path_length,
        area_below_euc = inter_euc$area_below_geo_baseline,
        area_between_euc = inter_euc$area_between_geo_peak,
        path_length_euc = inter_euc$path_length,
        surface = "Lifespan",
        trait = "Lifespan-Repr. rate",
        sample = j
      )
    )
  }
  
  ##########################################################
  ## 2. Lifespan-Repr. rate on Repr. rate surface ##########
  ##########################################################
  
  path_df <- subset(
    paths_fanson_all,
    trait == "Lifespan-Repr. rate" &
      surface == "Repr. rate" &
      sample == j
  )
  
  if (nrow(path_df) > 0) {
    
    inter_geo <- area_geodesic_with_peak(
      surface_df = Fanson_dailyeggs_surface,
      x_col_surf = "PROT",
      y_col_surf = "CARB",
      z_col_surf = "tps",
      path_df = path_df,
      standardize = TRUE,
      use_metric = TRUE,
      signed_area = FALSE
    )
    
    inter_euc <- area_geodesic_with_peak(
      surface_df = Fanson_dailyeggs_surface,
      x_col_surf = "PROT",
      y_col_surf = "CARB",
      z_col_surf = "tps",
      path_df = path_df,
      standardize = TRUE,
      use_metric = FALSE,
      signed_area = FALSE
    )
    
    area_df_fanson_trait_surface <- bind_rows(
      area_df_fanson_trait_surface,
      data.frame(
        area_below_geo = inter_geo$area_below_geo_baseline,
        area_between_geo = inter_geo$area_between_geo_peak,
        path_length_geo = inter_geo$path_length,
        area_below_euc = inter_euc$area_below_geo_baseline,
        area_between_euc = inter_euc$area_between_geo_peak,
        path_length_euc = inter_euc$path_length,
        surface = "Repr. rate",
        trait = "Lifespan-Repr. rate",
        sample = j
      )
    )
  }
  
  ##########################################################
  ## 3. Lifespan-lifetime eggs on Lifespan surface #########
  ##########################################################
  
  path_df <- subset(
    paths_fanson_all,
    trait == "Lifespan-lifetime eggs" &
      surface == "Lifespan" &
      sample == j
  )
  
  if (nrow(path_df) > 0) {
    
    inter_geo <- area_geodesic_with_peak(
      surface_df = Fanson_lifespan_surface,
      x_col_surf = "PROT",
      y_col_surf = "CARB",
      z_col_surf = "tps",
      path_df = path_df,
      standardize = TRUE,
      use_metric = TRUE,
      signed_area = FALSE
    )
    
    inter_euc <- area_geodesic_with_peak(
      surface_df = Fanson_lifespan_surface,
      x_col_surf = "PROT",
      y_col_surf = "CARB",
      z_col_surf = "tps",
      path_df = path_df,
      standardize = TRUE,
      use_metric = FALSE,
      signed_area = FALSE
    )
    
    area_df_fanson_trait_surface <- bind_rows(
      area_df_fanson_trait_surface,
      data.frame(
        area_below_geo = inter_geo$area_below_geo_baseline,
        area_between_geo = inter_geo$area_between_geo_peak,
        path_length_geo = inter_geo$path_length,
        area_below_euc = inter_euc$area_below_geo_baseline,
        area_between_euc = inter_euc$area_between_geo_peak,
        path_length_euc = inter_euc$path_length,
        surface = "Lifespan",
        trait = "Lifespan-lifetime eggs",
        sample = j
      )
    )
  }
  
  ##########################################################
  ## 4. Lifespan-lifetime eggs on Lifetime eggs surface ####
  ##########################################################
  
  path_df <- subset(
    paths_fanson_all,
    trait == "Lifespan-lifetime eggs" &
      surface == "Lifetime eggs" &
      sample == j
  )
  
  if (nrow(path_df) > 0) {
    
    inter_geo <- area_geodesic_with_peak(
      surface_df = Fanson_lifetimeegg_surface,
      x_col_surf = "PROT",
      y_col_surf = "CARB",
      z_col_surf = "tps",
      path_df = path_df,
      standardize = TRUE,
      use_metric = TRUE,
      signed_area = FALSE
    )
    
    inter_euc <- area_geodesic_with_peak(
      surface_df = Fanson_lifetimeegg_surface,
      x_col_surf = "PROT",
      y_col_surf = "CARB",
      z_col_surf = "tps",
      path_df = path_df,
      standardize = TRUE,
      use_metric = FALSE,
      signed_area = FALSE
    )
    
    area_df_fanson_trait_surface <- bind_rows(
      area_df_fanson_trait_surface,
      data.frame(
        area_below_geo = inter_geo$area_below_geo_baseline,
        area_between_geo = inter_geo$area_between_geo_peak,
        path_length_geo = inter_geo$path_length,
        area_below_euc = inter_euc$area_below_geo_baseline,
        area_between_euc = inter_euc$area_between_geo_peak,
        path_length_euc = inter_euc$path_length,
        surface = "Lifetime eggs",
        trait = "Lifespan-lifetime eggs",
        sample = j
      )
    )
  }
  
  ##########################################################
  ## 5. Lifetime eggs-Repr. rate on Lifetime surface #######
  ##########################################################
  
  path_df <- subset(
    paths_fanson_all,
    trait == "Lifetime eggs-Repr. rate" &
      surface == "Lifetime eggs" &
      sample == j
  )
  
  if (nrow(path_df) > 0) {
    
    inter_geo <- area_geodesic_with_peak(
      surface_df = Fanson_lifetimeegg_surface,
      x_col_surf = "PROT",
      y_col_surf = "CARB",
      z_col_surf = "tps",
      path_df = path_df,
      standardize = TRUE,
      use_metric = TRUE,
      signed_area = FALSE
    )
    
    inter_euc <- area_geodesic_with_peak(
      surface_df = Fanson_lifetimeegg_surface,
      x_col_surf = "PROT",
      y_col_surf = "CARB",
      z_col_surf = "tps",
      path_df = path_df,
      standardize = TRUE,
      use_metric = FALSE,
      signed_area = FALSE
    )
    
    area_df_fanson_trait_surface <- bind_rows(
      area_df_fanson_trait_surface,
      data.frame(
        area_below_geo = inter_geo$area_below_geo_baseline,
        area_between_geo = inter_geo$area_between_geo_peak,
        path_length_geo = inter_geo$path_length,
        area_below_euc = inter_euc$area_below_geo_baseline,
        area_between_euc = inter_euc$area_between_geo_peak,
        path_length_euc = inter_euc$path_length,
        surface = "Lifetime eggs",
        trait = "Lifetime eggs-Repr. rate",
        sample = j
      )
    )
  }
  
  ##########################################################
  ## 6. Lifetime eggs-Repr. rate on Repr. rate surface #####
  ##########################################################
  
  path_df <- subset(
    paths_fanson_all,
    trait == "Lifetime eggs-Repr. rate" &
      surface == "Repr. rate" &
      sample == j
  )
  
  if (nrow(path_df) > 0) {
    
    inter_geo <- area_geodesic_with_peak(
      surface_df = Fanson_dailyeggs_surface,
      x_col_surf = "PROT",
      y_col_surf = "CARB",
      z_col_surf = "tps",
      path_df = path_df,
      standardize = TRUE,
      use_metric = TRUE,
      signed_area = FALSE
    )
    
    inter_euc <- area_geodesic_with_peak(
      surface_df = Fanson_dailyeggs_surface,
      x_col_surf = "PROT",
      y_col_surf = "CARB",
      z_col_surf = "tps",
      path_df = path_df,
      standardize = TRUE,
      use_metric = FALSE,
      signed_area = FALSE
    )
    
    area_df_fanson_trait_surface <- bind_rows(
      area_df_fanson_trait_surface,
      data.frame(
        area_below_geo = inter_geo$area_below_geo_baseline,
        area_between_geo = inter_geo$area_between_geo_peak,
        path_length_geo = inter_geo$path_length,
        area_below_euc = inter_euc$area_below_geo_baseline,
        area_between_euc = inter_euc$area_between_geo_peak,
        path_length_euc = inter_euc$path_length,
        surface = "Repr. rate",
        trait = "Lifetime eggs-Repr. rate",
        sample = j
      )
    )
  }
}

write.csv(
  area_df_fanson_trait_surface,
  "/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Nutritional Geometry/Nutrigonometry/23. Geodesics/Final Scripts/Final codes/Revisions 1/Fanson analysis/Fanson_et_al_Geodesic_area_between_all_traits.csv",
  row.names = FALSE
)

############################################################
###### 9. Standardise costs by mean trait values ###########
############################################################

#mean_fanson_lifespan <- mean(fansondt$lifespan, na.rm = TRUE)
#mean_fanson_lifetimeegg <- mean(fansondt$lifetimeegg, na.rm = TRUE)
#mean_fanson_dailyeggs <- mean(fansondt$dailyeggs, na.rm = TRUE)

## because data has been standardized prior to calculations already
mean_fanson_lifespan <- 1
mean_fanson_lifetimeegg <- 1
mean_fanson_dailyeggs <- 1

area_df_fanson_std <- area_df_fanson_trait_surface %>%
  mutate(
    across(
      c(area_below_geo:path_length_euc),
      ~ case_when(
        surface == "Lifespan" ~ .x / mean_fanson_lifespan,
        surface == "Lifetime eggs" ~ .x / mean_fanson_lifetimeegg,
        surface == "Repr. rate" ~ .x / mean_fanson_dailyeggs,
        TRUE ~ .x
      )
    )
  )

############################################################
###### 10. Wilcoxon tests ##################################
############################################################

wilcox_fanson_lifespan <- wilcox.test(
  area_df_fanson_std$area_between_geo[area_df_fanson_std$surface == "Lifespan"],
  area_df_fanson_std$area_between_euc[area_df_fanson_std$surface == "Lifespan"],
  paired = TRUE
)

wilcox_fanson_lifetimeegg <- wilcox.test(
  area_df_fanson_std$area_between_geo[area_df_fanson_std$surface == "Lifetime eggs"],
  area_df_fanson_std$area_between_euc[area_df_fanson_std$surface == "Lifetime eggs"],
  paired = TRUE
)

wilcox_fanson_dailyeggs <- wilcox.test(
  area_df_fanson_std$area_between_geo[area_df_fanson_std$surface == "Repr. rate"],
  area_df_fanson_std$area_between_euc[area_df_fanson_std$surface == "Repr. rate"],
  paired = TRUE
)

wilcox_fanson_lifespan
wilcox_fanson_lifetimeegg
wilcox_fanson_dailyeggs


wilcox.test(
  area_df_fanson_std$area_between_geo[area_df_fanson_std$surface == "Repr. rate"],
  area_df_fanson_std$area_between_geo[area_df_fanson_std$surface == "Lifetime eggs"],
  paired = FALSE
)

wilcox.test(
  area_df_fanson_std$area_between_euc[area_df_fanson_std$surface == "Repr. rate"],
  area_df_fanson_std$area_between_euc[area_df_fanson_std$surface == "Lifetime eggs"],
  paired = FALSE
)
############################################################
###### 11. Summary table ###################################
############################################################

area_fanson_long <- area_df_fanson_std %>%
  pivot_longer(
    cols = c(area_between_geo, area_between_euc),
    names_to = "Type",
    values_to = "Area"
  ) %>%
  mutate(
    Type = ifelse(Type == "area_between_geo", "Geodesic", "Euclidean"),
    surface = case_when(
      surface == "Lifespan" ~ "Surface:\nLifespan",
      surface == "Lifetime eggs" ~ "Surface:\nLifetime eggs",
      surface == "Repr. rate" ~ "Surface:\nRepr. rate"
    ),
    trait = case_when(
      trait == "Lifespan-Repr. rate" ~ "Lifespan\nvs\nRepr. rate",
      trait == "Lifespan-lifetime eggs" ~ "Lifespan\nvs\nLifetime eggs",
      trait == "Lifetime eggs-Repr. rate" ~ "Lifetime eggs\nvs\nRepr. rate"
    )
  )

area_fanson_summary <- area_fanson_long %>%
  group_by(Type, trait, surface) %>%
  summarise(
    avg = mean(Area, na.rm = TRUE),
    std = sd(Area, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

area_fanson_summary

write.csv(
  area_fanson_summary,
  "/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Nutritional Geometry/Nutrigonometry/23. Geodesics/Final Scripts/Final codes/Revisions 1/Fanson analysis/Fanson_et_al_Geodesic_Euclidean_cost_summary.csv",
  row.names = FALSE
)

area_fanson_summary %>%
  arrange(Type, surface) 
############################################################
###### 12. Cost plot #######################################
############################################################

costs_plot_fanson <- ggplot(
  area_fanson_summary,
  aes(x = Type, y = avg, fill = Type)
) +
  facet_grid(surface ~ trait) +
  geom_col(col = "black", width = 0.9) +
  geom_errorbar(
    aes(ymin = avg - std, ymax = avg + std),
    width = 0.1
  ) +
  ylab(expression(italic(C[xy]))) +
  xlab("") +
  theme_linedraw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
    axis.title = element_text(size = 13),
    panel.grid = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "lavenderblush2"),
    strip.text = element_text(size = 10, color = "black")
  ) +
  scale_fill_manual(
    "Type",
    values = c("lightskyblue2", "indianred2")
  )

costs_plot_fanson


############################################################
###### 13. Representative centroid geodesics ###############
############################################################

# Euclidean paths
euclidean_path_fanson_lifespan_dailyeggs <- euclidean_path_2d(
  p1 = Fanson_lifespan_centroid,
  p2 = Fanson_dailyeggs_centroid
)

euclidean_path_fanson_lifetime_dailyeggs <- euclidean_path_2d(
  p1 = Fanson_lifetimeegg_centroid,
  p2 = Fanson_dailyeggs_centroid
)

euclidean_path_fanson_lifetime_lifespan <- euclidean_path_2d(
  p1 = Fanson_lifetimeegg_centroid,
  p2 = Fanson_lifespan_centroid
)

# Lifespan vs Repr. rate on lifespan surface
a <- fit_fanson_lifespan$coefs[["a"]]
b <- fit_fanson_lifespan$coefs[["b"]]
c <- fit_fanson_lifespan$coefs[["c"]]
d <- fit_fanson_lifespan$coefs[["d"]]
e <- fit_fanson_lifespan$coefs[["e"]]

best_fanson_lifespan_dailyeggs <- robust_optimize_geodesic(
  Fanson_lifespan_centroid[[1]],
  Fanson_lifespan_centroid[[2]],
  Fanson_dailyeggs_centroid[[1]],
  Fanson_dailyeggs_centroid[[2]]
)

sol_fanson_lifespan_dailyeggs <- integrate_geodesic(
  best_fanson_lifespan_dailyeggs,
  Fanson_lifespan_centroid[[1]],
  Fanson_lifespan_centroid[[2]]
)

# Lifetime eggs vs lifespan on lifetime surface
a <- fit_fanson_lifetimeegg$coefs[["a"]]
b <- fit_fanson_lifetimeegg$coefs[["b"]]
c <- fit_fanson_lifetimeegg$coefs[["c"]]
d <- fit_fanson_lifetimeegg$coefs[["d"]]
e <- fit_fanson_lifetimeegg$coefs[["e"]]

best_fanson_lifetime_lifespan <- robust_optimize_geodesic(
  Fanson_lifetimeegg_centroid[[1]],
  Fanson_lifetimeegg_centroid[[2]],
  Fanson_lifespan_centroid[[1]],
  Fanson_lifespan_centroid[[2]]
)

sol_fanson_lifetime_lifespan <- integrate_geodesic(
  best_fanson_lifetime_lifespan,
  Fanson_lifetimeegg_centroid[[1]],
  Fanson_lifetimeegg_centroid[[2]]
)

# Lifetime eggs vs Repr. rate on lifetime surface
best_fanson_lifetime_dailyeggs <- robust_optimize_geodesic(
  Fanson_lifetimeegg_centroid[[1]],
  Fanson_lifetimeegg_centroid[[2]],
  Fanson_dailyeggs_centroid[[1]],
  Fanson_dailyeggs_centroid[[2]]
)

sol_fanson_lifetime_dailyeggs <- integrate_geodesic(
  best_fanson_lifetime_dailyeggs,
  Fanson_lifetimeegg_centroid[[1]],
  Fanson_lifetimeegg_centroid[[2]]
)

############################################################
###### 14. Path plots ######################################
############################################################

Fanson_plot_lifespan_geo <- Fanson_lifespan_plot +
  geom_point(
    data = Fanson_dailyeggs_centroid,
    mapping = aes(x = PROT, y = CARB),
    pch = 21,
    fill = "firebrick2",
    col = "black",
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_path(
    data = sol_fanson_lifespan_dailyeggs,
    mapping = aes(x = x, y = y),
    color = "indianred2",
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_path(
    data = euclidean_path_fanson_lifespan_dailyeggs,
    mapping = aes(x = x, y = y),
    color = "lightskyblue2",
    linetype = "dashed",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  ggtitle("Lifespan vs Repr. rate", subtitle = "Landscape: Lifespan") + 
  theme(plot.title = element_text(hjust = 0.5),   # Centers Title
        plot.subtitle = element_text(hjust = 0.5))+ 
  geo_panel_theme 

Fanson_plot_lifetime_lifespan_geo <- Fanson_lifetimeegg_plot +
  geom_point(
    data = Fanson_lifespan_centroid,
    mapping = aes(x = PROT, y = CARB),
    pch = 21,
    fill = "steelblue2",
    col = "black",
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_path(
    data = sol_fanson_lifetime_lifespan,
    mapping = aes(x = x, y = y),
    color = "indianred2",
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_path(
    data = euclidean_path_fanson_lifetime_lifespan,
    mapping = aes(x = x, y = y),
    color = "lightskyblue2",
    linetype = "dashed",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  ggtitle("Lifetime eggs vs Lifespan", subtitle = "Landscape: Lifetime eggs") + 
  theme(plot.title = element_text(hjust = 0.5),   # Centers Title
        plot.subtitle = element_text(hjust = 0.5))  + 
  geo_panel_theme 

Fanson_plot_lifetime_dailyeggs_geo <- Fanson_lifetimeegg_plot +
  geom_point(
    data = Fanson_dailyeggs_centroid,
    mapping = aes(x = PROT, y = CARB),
    pch = 21,
    fill = "firebrick2",
    col = "black",
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_path(
    data = sol_fanson_lifetime_dailyeggs,
    mapping = aes(x = x, y = y),
    color = "indianred2",
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_path(
    data = euclidean_path_fanson_lifetime_dailyeggs,
    mapping = aes(x = x, y = y),
    color = "lightskyblue2",
    linetype = "dashed",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  ggtitle("Lifetime eggs vs Repr. rate", subtitle = "Landscape: Lifetime eggs") + 
  theme(plot.title = element_text(hjust = 0.5),   # Centers Title
        plot.subtitle = element_text(hjust = 0.5)) + 
  geo_panel_theme  

Fanson_path_figure <- (
  Fanson_plot_lifespan_geo +
    Fanson_plot_lifetime_lifespan_geo +
    Fanson_plot_lifetime_dailyeggs_geo
)

Fanson_path_figure









### For reviews
Fanson_path_figure_Final <- (
    Fanson_plot_lifetime_lifespan_geo +
    Fanson_plot_lifetime_dailyeggs_geo
)
Fanson_path_figure_Final


area_fanson_summary$trait <- ifelse(area_fanson_summary$trait == "Lifetime eggs\nvs\nRepr. rate",
                                    "Lifetime eggs\nvs\nRepr. rate",
                                    area_fanson_summary$trait)

area_fanson_summary$trait <- ifelse(area_fanson_summary$trait == "Lifespan\nvs\nRepr. rate",
                                    "Lifespan\nvs\nRepr. rate",
                                    area_fanson_summary$trait)


costs_plot_fanson_Final <- ggplot(
  area_fanson_summary[area_fanson_summary$surface == "Surface:\nLifetime eggs",],
  aes(x = Type, y = avg, fill = Type)
) +
  facet_grid( ~ trait) +
  geom_col(col = "black", width = 0.9) +
  geom_errorbar(
    aes(ymin = avg - std, ymax = avg + std),
    width = 0.1
  ) +
  #ylab(expression(italic(C[xy]))) +
  ylab("Cost") +
  xlab("Type") +
  theme_linedraw() + 
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5),
        axis.title  = element_text(size = 16),
        panel.grid  = element_blank(),
        legend.position = "none",
        aspect.ratio = 2,
        strip.background = element_rect(fill = "lavenderblush2"),
        strip.text      = element_text(size = 13, color = "black")) + 
  scale_fill_manual('Type', values = c("lightskyblue2", "indianred2"))

costs_plot_fanson_Final


Fanson_path_figure_Final / costs_plot_fanson_Final




area_fanson_summary %>%
  group_by(surface, Type) %>%
  summarise(mean(avg)) %>%
  arrange(Type)


dir()

############################################################
###### Save workspace ######################################
############################################################

save.image(
  file = "Fanson_geodesic_workspace.RData"
)

