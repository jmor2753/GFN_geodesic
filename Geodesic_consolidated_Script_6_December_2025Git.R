# --- Load packages ---
library(ggplot2)
library(dplyr)
library(optimx)
library(deSolve)
library(terra)
library(sf)

# --- Load data and nutrigonometry functions ---
source("/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Nutritional Geometry/Nutrigonometry/23. Geodesics/Final Scripts/Final codes/Geodesic_consolidated_Functions_6_December_2025.R")
set.seed(321908)
# --- PARAMETERS ---
a <- -0.2; b <- -0.2; c <- 0.2; d <- 0.2; e <- 0.3
x0 <- 0.5; y0 <- 9.5
target <- c(9.5, 0.5)


# --- SCALAR FIELD --
# --- SURFACE P(x, y) BACKGROUND ---
x_vals <- seq(0, 10, by = 0.05)
y_vals <- seq(0, 10, by = 0.05)
grid <- expand.grid(x = x_vals, y = y_vals)
grid$P <- with(grid, P(x, y))

## PLOT 
ggplot() +
  geom_tile(data = grid, aes(x = x, y = y, fill = P), color = NA) +
  scale_fill_viridis_c()

### functions
# --- USAGE EXAMPLE ---
x0 <- 2.5; y0 <- 7.5
target <- c(7.5, 2.5)


best_vT <- robust_optimize_geodesic(x0, y0, target[1], target[2]) ## may take a while
sol_df <- integrate_geodesic(best_vT, x0, y0)





# --- SURFACE P(x, y) BACKGROUND ---
x_vals <- seq(0, 10, by = 0.05)
y_vals <- seq(0, 10, by = 0.05)
grid <- expand.grid(x = x_vals, y = y_vals)
grid$P <- with(grid, P(x, y))




# --- FINAL PLOT ---
ggplot() +
  geom_tile(data = grid, aes(x = x, y = y, fill = P), color = NA) +
  geom_contour(data = grid, aes(x = x, y = y, z = P), col = "grey50", linewidth = 0.2) + 
  scale_fill_viridis_c() +
  geom_path(data = sol_df, aes(x = x, y = y), color = "red", linewidth = 0.5) +
  geom_point(aes(x = x0, y = y0), color = "white", size = 2) +
  geom_point(aes(x = target[1], y = target[2]), color = "orange", size = 2) +
  coord_equal() +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        legend.position = "none")






# --- EUCLIDEAN PATH ---
euclidean_path <- data.frame(
  x = seq(x0, target[1], length.out = 100),
  y = seq(y0, target[2], length.out = 100)
)



ggplot() +
  geom_tile(data = grid, aes(x = x, y = y, fill = P), color = NA) +
  geom_contour(data = grid, aes(x = x, y = y, z = P), col = "grey50", linewidth = 0.2) + 
  scale_fill_viridis_c() +
  geom_path(data = sol_df, aes(x = x, y = y), color = "red", linewidth = 0.5) +
  geom_point(aes(x = x0, y = y0), color = "white", size = 2) +
  geom_point(aes(x = target[1], y = target[2]), color = "orange", size = 2) +
  geom_path(data = euclidean_path, aes(x = x, y = y), color = "black", linetype = "dashed", linewidth = 0.5) +
  coord_equal() +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        legend.position = "none")




# --- UPDATE DISTANCE STATS ---


geo_dist <- geodesic_distance(sol_df)
euc_dist <- sqrt(sum((target - c(x0, y0))^2))
euclidean_path <- data.frame(
  x = seq(x0, target[1], length.out = 100),
  y = seq(y0, target[2], length.out = 100)
)




# Euclidean distance
cat(sprintf("Geodesic Distance: %.4f\n", geo_dist))
cat(sprintf("Euclidean Distance: %.4f\n", euc_dist))
cat(sprintf("Geodesic / Euclidean Ratio: %.4f\n", geo_dist / euc_dist))





# --- DISTANCE ANNOTATION ---
label_text <- sprintf("Geodesic: %.2f\nEuclidean: %.2f\nRatio: %.2f", 
                      geo_dist, euc_dist, geo_dist / euc_dist)

# Updated plot with annotation


ggplot() +
  geom_tile(data = grid, aes(x = x, y = y, fill = P), color = NA) +
  geom_contour(data = grid, aes(x = x, y = y, z = P), col = "grey50", linewidth = 0.2) + 
  scale_fill_viridis_c() +
  geom_path(data = sol_df, aes(x = x, y = y), color = "red", linewidth = 0.5) +
  geom_point(aes(x = x0, y = y0), color = "white", size = 2) +
  geom_point(aes(x = target[1], y = target[2]), color = "orange", size = 2) +
  geom_path(data = euclidean_path, aes(x = x, y = y), color = "black", linetype = "dashed", linewidth = 0.5) +
  coord_equal() +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") + 
  annotate("text", x = x0 + 1, y = y0 - 1, label = label_text, hjust = 0, size = 4)



##### ~~~ simulating multiple landscapes and calculating geodesic, then plotting

generate_landscape_batch <- function(n = 15, seed = 42) {
  set.seed(seed)
  
  landscapes <- list()
  
  for (i in 1:n) {
    # 1. Random parameters for surface
    a <- runif(1, -0.3, 0.3)
    b <- runif(1, -0.3, 0.3)
    c <- runif(1, -0.3, 0.3)
    d <- runif(1, -0.3, 0.3)
    e <- runif(1, -0.3, 0.3)
    
    # 2. Start/target points (fixed or randomized)
    x0 <- runif(1, 1, 3)
    y0 <- runif(1, 7, 9)
    x1 <- runif(1, 7, 9)
    y1 <- runif(1, 1, 3)
    
    # 3. Define scalar field and gradient functions
    P <- function(x, y) a*x^2 + b*y^2 + c*x + d*y + e*x*y
    P_grad <- function(x, y) c(2*a*x + c + e*y, 2*b*y + d + e*x)
    metric_tensor <- function(x, y) {
      grad <- P_grad(x, y)
      Px <- grad[1]; Py <- grad[2]
      matrix(c(1 + Px^2, Px * Py, Px * Py, 1 + Py^2), nrow = 2)
    }
    normalize_direction_signed <- function(x, y, dir) {
      g <- metric_tensor(x, y)
      len2 <- as.numeric(t(dir) %*% g %*% dir)
      scale <- 1 / sqrt(len2)
      signed_dir <- as.numeric(scale * dir)
      return(signed_dir)
    }
    
    generate_initial_guess <- function(x0, y0, x1, y1) {
      dir <- c(x1 - x0, y1 - y0)
      norm <- sqrt(sum(dir^2))
      if (norm == 0) norm <- 1
      v_unit <- dir / norm
      T0 <- norm * 2
      return(c(v_unit * T0, T0))
    }
    
    geodesic_eqs <- function(t, state, parameters) {
      x <- state[1]; dx <- state[2]
      y <- state[3]; dy <- state[4]
      
      Px <- 2*a*x + c + e*y
      Py <- 2*b*y + d + e*x
      
      g <- metric_tensor(x, y)
      speed2 <- dx*(g[1,1]*dx + g[1,2]*dy) + dy*(g[2,1]*dx + g[2,2]*dy)
      speed <- sqrt(speed2)
      if (abs(speed - 1) > 1e-6 || t %% 0.01 < 1e-5) {
        dx <- dx / speed
        dy <- dy / speed
      }
      
      D <- 4*a^2*x^2 + 4*a*c*x + 4*a*e*x*y +
        4*b^2*y^2 + 4*b*d*y + 4*b*e*x*y +
        c^2 + 2*c*e*y + d^2 + 2*d*e*x +
        e^2*x^2 + e^2*y^2 + 1
      D <- ifelse(abs(D) < 1e-6, 1e-6, D)
      D <- min(D, 1e6)
      
      num_x <- a * dx^2 * (4*a*x + 2*c + 2*e*y) +
        b * dy^2 * (4*a*x + 2*c + 2*e*y) +
        2*e * dx * dy * (2*a*x + c + e*y)
      
      num_y <- a * dx^2 * (4*b*y + 2*d + 2*e*x) +
        b * dy^2 * (4*b*y + 2*d + 2*e*x) +
        2*e * dx * dy * (2*b*y + d + e*x)
      
      ddx <- -num_x / D
      ddy <- -num_y / D
      list(c(dx, ddx, dy, ddy))
    }
    
    integrate_geodesic <- function(vT, x0, y0) {
      v0 <- vT[1:2]
      T <- abs(vT[3])
      dir_metric <- normalize_direction_signed(x0, y0, v0)
      state0 <- c(x = x0, dx = dir_metric[1], y = y0, dy = dir_metric[2])
      times <- seq(0, T, by = 0.001)
      ode(y = state0, times = times, func = geodesic_eqs, parms = NULL, method = "lsoda") %>%
        as.data.frame()
    }
    
    geodesic_objective <- function(vT, x0, y0, x1, y1) {
      df <- tryCatch(integrate_geodesic(vT, x0, y0), error = function(e) return(NULL))
      if (is.null(df)) return(Inf)
      endpoint <- df[nrow(df), c("x", "y")]
      sum((endpoint - c(x1, y1))^2)
    }
    
    robust_optimize_geodesic <- function(x0, y0, x1, y1) {
      guess <- generate_initial_guess(x0, y0, x1, y1)
      res <- tryCatch({
        optimx::optimx(
          par = guess,
          fn = function(vT) geodesic_objective(vT, x0, y0, x1, y1),
          method = "L-BFGS-B",
          lower = c(-20, -20, 1),
          upper = c(20, 20, 50),
          control = list(maxit = 100)
        )
      }, error = function(e) return(NULL))
      if (is.null(res) || any(is.na(res[1, c("p1", "p2", "p3")]))) return(NULL)
      as.numeric(res[1, c("p1", "p2", "p3")])
    }
    
    # Surface grid
    x_vals <- seq(0, 10, by = 0.05)
    y_vals <- seq(0, 10, by = 0.05)
    grid <- expand.grid(x = x_vals, y = y_vals)
    grid$P <- with(grid, P(x, y))
    
    # Geodesic path
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    if (is.null(vT)) next
    geodesic_df <- integrate_geodesic(vT, x0, y0)
    
    # Euclidean path
    euclidean_df <- data.frame(
      x = seq(x0, x1, length.out = 100),
      y = seq(y0, y1, length.out = 100)
    )
    
    # Metric distance function
    geodesic_distance <- function(df) {
      total <- 0
      for (i in 2:nrow(df)) {
        dx_vec <- c(df$x[i] - df$x[i-1], df$y[i] - df$y[i-1])
        mid_x <- (df$x[i] + df$x[i-1]) / 2
        mid_y <- (df$y[i] + df$y[i-1]) / 2
        g <- metric_tensor(mid_x, mid_y)
        ds2 <- t(dx_vec) %*% g %*% dx_vec
        total <- total + sqrt(ds2)
      }
      return(as.numeric(total))
    }
    
    geo_dist <- geodesic_distance(geodesic_df)
    euc_dist <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
    
    # Save everything
    landscapes[[i]] <- list(
      params = c(a, b, c, d, e),
      start = c(x0, y0),
      target = c(x1, y1),
      surface = grid,
      geodesic = geodesic_df,
      euclidean = euclidean_df,
      geo_dist = geo_dist,
      euc_dist = euc_dist
    )
  }
  
  return(landscapes)
}



# Generate data
## from upper right to lower left
##landscapes <- generate_landscape_batch(n = 15)
##saveRDS(landscapes, file = "landscapes_batch_1.rds")


## from lower left to upper right
##landscapes_2 <- generate_landscape_batch(n = 5, seed = 43)
##saveRDS(landscapes_2, file = "landscapes_batch_2.rds")


## random
##landscapes_3 <- generate_landscape_batch(n = 5, seed = 44)
##saveRDS(landscapes_3, file = "landscapes_batch_3.rds")

## total random
##landscapes_4 <- generate_landscape_batch(n = 7, seed = 45)
##saveRDS(landscapes_4, file = "landscapes_batch_4.rds")
#landscapes_4 <- readRDS("/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Nutritional Geometry/Nutrigonometry/23. Geodesics/rds landscapes/landscapes_batch_4.rds")

# Plotting
for (i in seq_along(landscapes_4)) {
  l <- landscapes_4[[i]]
  
  param_str <- paste0("a=", round(l$params[1], 2), ", b=", round(l$params[2], 2),
                      ", c=", round(l$params[3], 2), ", d=", round(l$params[4], 2),
                      ", e=", round(l$params[5], 2))
  
  dist_str <- sprintf("Start:(%.1f,%.1f)  →  Target:(%.1f,%.1f)\nGeo=%.2f, Euc=%.2f, Ratio=%.2f",
                      l$start[1], l$start[2], l$target[1], l$target[2],
                      l$geo_dist, l$euc_dist, l$geo_dist / l$euc_dist)
  
  p <- ggplot() +
    geom_tile(data = l$surface, aes(x = x, y = y, fill = P), color = NA) +
    geom_contour(data = l$surface, aes(x = x, y = y, z = P), color = "grey40", linewidth = 0.2) +
    geom_path(data = l$geodesic, aes(x = x, y = y), color = "red", linewidth = 0.7) +
    geom_path(data = l$euclidean, aes(x = x, y = y), color = "black", linetype = "dashed") +
    geom_point(aes(x = l$start[1], y = l$start[2]), color = "white", size = 2) +
    geom_point(aes(x = l$target[1], y = l$target[2]), color = "orange", size = 2) +
    scale_fill_viridis_c() +
    coord_equal() +
    theme_minimal() +
    theme(panel.grid = element_blank(), legend.position = "none") +
    ggtitle(paste("Landscape", i), subtitle = paste(param_str, "\n", dist_str))
  
  print(p)
}











###### **  LEE ET AL CALIDATION ** ############
leedt <- read.csv("/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Nutritional Geometry/Nutrigonometry/1-Right-angle-triangle/Dataset_Lee_Nutrigonometry.csv", 
                  header = TRUE, 
                  strip.white = TRUE, 
                  stringsAsFactors = FALSE,
                  na.string = 'NA') %>%
  na.omit()

### Data pre-processing
leedt$Food <- as.factor(leedt$Food)
leedt$CONC <- as.factor(ifelse(leedt$Food == "45", 'x0.25', 
                               ifelse(leedt$Food == "90", "x0.5",
                                      ifelse(leedt$Food == "180", 'x1', 'x2'))))


# Adjusting Ratio 
leedt$RATIO <- as.factor(stringr::str_replace_all(leedt$Ratio, '\\(|\\)', ''))
levels(leedt$RATIO)

# Adjusting names of carb and prot variables
leedt$carb_eaten <- as.numeric(leedt$carb_eaten)
leedt$protein_eaten <- as.numeric(leedt$protein_eaten)

# Adjusting other numeric vars
leedt$lifespan <- as.numeric(leedt$lifespan) 
leedt$lifetimeegg <- as.numeric(leedt$lifetimeegg)
leedt$dailyeggs <- as.numeric(leedt$dailyeggs)



# --- Rescale x and y columns ---
leedt <- leedt %>%
  mutate(
    P_scaled = scale(P_fixed, scale = TRUE)[,1],
    C_scaled = scale(C_fixed, scale = TRUE)[,1]
  )


# --- Run ---
fit_lifespan <- fit_quadratic_surface(leedt, "P_scaled", "C_scaled", "lifespan")
fit_dailyeggs <- fit_quadratic_surface(leedt, "P_scaled", "C_scaled", "dailyeggs")
fit_lifetimeegg <- fit_quadratic_surface(leedt, "P_scaled", "C_scaled", "lifetimeegg")

### ---- Peak
# * Lifespan * #
predictonML_lm_lifespan_fixed <- predict_region_GF_lm(leedt,
                                                      Pcol = 'P_scaled',
                                                      Ccol = 'C_scaled', 
                                                      Zcol = 'lifespan',
                                                      peak = 'yes')


## Lifetime eggs
predictonML_lm_lifetimeegg_fixed <- predict_region_GF_lm(leedt,
                                                         Pcol = 'P_scaled',
                                                         Ccol = 'C_scaled', 
                                                         Zcol = 'lifetimeegg',
                                                         peak = 'yes')


# * Reproductive rate * #
predictonML_lm_dailyeggs_fixed <- predict_region_GF_lm(leedt,
                                                       Pcol = 'P_scaled',
                                                       Ccol = 'C_scaled', 
                                                       Zcol = 'dailyeggs',
                                                       peak = 'yes')


## Extracting parameters

### Creating predicted convex hulld ####
# * Lifespan * #
hullLM_lifespan_fixed <- design_hull(predictonML_lm_lifespan_fixed)
hullLM_lifespan_fixed
##Centroid
## MIP centroid variation
Lifespan_centroid <- as.data.frame(terra::crds(terra::centroids(convHull(vect(cbind(x = hullLM_lifespan_fixed$PROT, 
                                                                                    y = hullLM_lifespan_fixed$CARB))))))
colnames(Lifespan_centroid) <- c("PROT", "CARB")


## * Lifetime egg
hullLM_lifetimeegg_fixed <- design_hull(predictonML_lm_lifetimeegg_fixed)
hullLM_lifetimeegg_fixed
##Centroid
## MIP centroid variation
lifetimeegg_centroid <- as.data.frame(terra::crds(terra::centroids(convHull(vect(cbind(x = hullLM_lifetimeegg_fixed$PROT, 
                                                                                       y = hullLM_lifetimeegg_fixed$CARB))))))
colnames(lifetimeegg_centroid) <- c("PROT", "CARB")



# * Reproductive Rate * #
hullLM_RR_fixed <- design_hull(predictonML_lm_dailyeggs_fixed)
## MIP centroid variation
RR_centroid <- as.data.frame(terra::crds(terra::centroids(convHull(vect(cbind(x = hullLM_RR_fixed$PROT, 
                                                                              y = hullLM_RR_fixed$CARB))))))
colnames(RR_centroid) <- c("PROT", "CARB")




### scaled surfaces
## Predicting landscapes for traits ##
# * Source the function scripts * #
Leedt_Survsurface_fixed_scaled <- predict_surface(leedt, 
                                                  x1 = 'P_scaled',
                                                  x2 = 'C_scaled', 
                                                  y = 'lifespan', lambda = 0.5)
Leedt_Lifetimeeggsurface_fixed_scaled <- predict_surface(leedt, 
                                                         x1 = 'P_scaled',
                                                         x2 = 'C_scaled', 
                                                         y = 'lifetimeegg', lambda = 0.5)

Leedt_RRsurface_fixed_scaled <- predict_surface(leedt, 
                                                x1 = 'P_scaled',
                                                x2 = 'C_scaled', 
                                                y = 'dailyeggs', lambda = 0.5)




#  * Lifespan * #
Survlandscape_plot_fixed_scaled <- ggplot(data = Leedt_Survsurface_fixed_scaled, 
                                          aes(x = PROT, 
                                              y = CARB, 
                                              z = tps, 
                                              fill = tps))+
  geom_tile(alpha = 0.9) + 
  geom_contour(data = Leedt_Survsurface_fixed_scaled, 
               aes(x = PROT, y = CARB, z = tps), 
               color = "black", 
               breaks = seq(from = 0, to = 30, by = 2), 
               size = 0.08) +
  ylab("Carbohydrates (mg)") + 
  xlab("Protein (mg)") + 
  ggtitle('') + 
  theme_bw() + 
  theme(panel.grid=element_blank(),
        aspect.ratio = 1.2,
        plot.title = element_text(hjust = 0.5, face = 'bold', size = 11),
        axis.text = element_text(size = 9, color = 'black'),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        legend.justification=c(1,1),
        legend.title = element_text(size=6), 
        legend.text = element_text(size = 6)) +
  scale_fill_distiller(palette = "Spectral", 
                       name = 'Lifespan (days)',   
                       limits = c(min(Leedt_Survsurface_fixed_scaled$tps), 
                                  max(Leedt_Survsurface_fixed_scaled$tps))) +
  geom_polygon(hullLM_lifespan_fixed, mapping = aes(x = PROT, y = CARB),
               fill = 'black',
               col = 'black', 
               inherit.aes = FALSE, 
               alpha = 0.1,
               size = 0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 11),
        axis.text = element_text(size = 7, color = 'black'),
        axis.title = element_text(size = 10)) + 
  geom_point(data = Lifespan_centroid, mapping = aes(x = PROT, y = CARB), 
             pch = 21, 
             fill = "steelblue2", 
             col = "black",
             inherit.aes = FALSE,
             size = 3)
Survlandscape_plot_fixed_scaled




###* Lifetime Egg plot

Lifetimeegglandscape_plot_fixed_scaled <- ggplot(data = Leedt_Lifetimeeggsurface_fixed_scaled, 
                                                 aes(x = PROT, 
                                                     y = CARB, 
                                                     z = tps, 
                                                     fill = tps))+
  geom_tile(alpha = 0.9) + 
  geom_contour(data = Leedt_Lifetimeeggsurface_fixed_scaled, 
               aes(x = PROT, y = CARB, z = tps), 
               color = "black", 
               breaks = seq(from = 0, to =100, by = 8), 
               size = 0.08) +
  ylab("Carbohydrates (mg)") + 
  xlab("Protein (mg)") + 
  ggtitle('') + 
  theme_bw() + 
  theme(panel.grid=element_blank(),
        aspect.ratio = 1.2,
        plot.title = element_text(hjust = 0.5, face = 'bold', size = 11),
        axis.text = element_text(size = 9, color = 'black'),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        legend.justification=c(1,1),
        legend.title = element_text(size=6), 
        legend.text = element_text(size = 6)) +
  scale_fill_distiller(palette = "Spectral", 
                       name = 'Lifespan (days)',   
                       limits = c(min(Leedt_Lifetimeeggsurface_fixed_scaled$tps), 
                                  max(Leedt_Lifetimeeggsurface_fixed_scaled$tps))) +
  geom_polygon(hullLM_lifetimeegg_fixed, mapping = aes(x = PROT, y = CARB),
               fill = 'black',
               col = 'black', 
               inherit.aes = FALSE, 
               alpha = 0.1,
               size = 0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 11),
        axis.text = element_text(size = 7, color = 'black'),
        axis.title = element_text(size = 10)) + 
  geom_point(data = lifetimeegg_centroid, mapping = aes(x = PROT, y = CARB), 
             pch = 21, 
             fill = "orange", 
             col = "black",
             inherit.aes = FALSE,
             size = 3)
Lifetimeegglandscape_plot_fixed_scaled

# * Reproductive Rate * #
RRlandscape_plot_fixed_scaled <- ggplot(data = Leedt_RRsurface_fixed_scaled, 
                                        aes(x = PROT, 
                                            y = CARB, 
                                            z = tps, 
                                            fill = tps))+
  geom_tile(alpha = 0.9) + 
  geom_contour(data = Leedt_RRsurface_fixed_scaled, 
               aes(x = PROT, y = CARB, z = tps), 
               color = "black", 
               breaks = seq(from = 0, to = 15, by = 0.5), 
               size = 0.08) +
  ylab("Carbohydrates (mg)") + 
  xlab("Protein (mg)") + 
  ggtitle('') + 
  theme_bw() + 
  theme(panel.grid=element_blank(),
        aspect.ratio = 1.2,
        plot.title = element_text(hjust = 0.5, face = 'bold', size = 11),
        axis.text = element_text(size = 9, color = 'black'),
        axis.title = element_text(size = 10),
        legend.position = 'none',
        legend.justification=c(1,1),
        legend.title = element_text(size=6), 
        legend.text = element_text(size = 6)) +
  scale_fill_distiller(palette = "Spectral", 
                       name = 'Lifespan (days)',   
                       limits = c(min(Leedt_RRsurface_fixed_scaled$tps), 
                                  max(Leedt_RRsurface_fixed_scaled$tps))) + 
  geom_polygon(hullLM_RR_fixed, mapping = aes(x = PROT, y = CARB),
               fill = 'black',
               col = 'black', 
               inherit.aes = FALSE, 
               alpha = 0.1,
               size = 0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 11),
        axis.text = element_text(size = 7, color = 'black'),
        axis.title = element_text(size = 10)) + 
  geom_point(data = RR_centroid, mapping = aes(x = PROT, y = CARB), 
             pch = 21, 
             fill = "firebrick2", 
             col = "black",
             inherit.aes = FALSE,
             size = 3)
RRlandscape_plot_fixed_scaled



##geodesics

##lifespan
a <- fit_lifespan$coefs[["a"]]
b <- fit_lifespan$coefs[["b"]]
c <- fit_lifespan$coefs[["c"]]
d <- fit_lifespan$coefs[["d"]]
e <- fit_lifespan$coefs[["e"]]

best_vT_lifespan <- robust_optimize_geodesic(Lifespan_centroid[[1]], Lifespan_centroid[[2]], 
                                             RR_centroid[[1]], RR_centroid[[2]]) ## may take a while


sol_df_lifespan  <- integrate_geodesic(best_vT_lifespan, Lifespan_centroid[[1]], Lifespan_centroid[[2]])



## dailyeggs
a <- fit_dailyeggs$coefs[["a"]]
b <- fit_dailyeggs$coefs[["b"]]
c <- fit_dailyeggs$coefs[["c"]]
d <- fit_dailyeggs$coefs[["d"]]
e <- fit_dailyeggs$coefs[["e"]]

best_vT_dailyeggs <- robust_optimize_geodesic(Lifespan_centroid[[1]], Lifespan_centroid[[2]], 
                                              RR_centroid[[1]], RR_centroid[[2]]) ## may take a while


sol_df_dailyeggs <- integrate_geodesic(best_vT_dailyeggs, Lifespan_centroid[[1]], Lifespan_centroid[[2]])



##lifetime eggs
a <- fit_lifetimeegg$coefs[["a"]]
b <- fit_lifetimeegg$coefs[["b"]]
c <- fit_lifetimeegg$coefs[["c"]]
d <- fit_lifetimeegg$coefs[["d"]]
e <- fit_lifetimeegg$coefs[["e"]]

best_vT_lifetimeegg_lifespan <- robust_optimize_geodesic(lifetimeegg_centroid[[1]], lifetimeegg_centroid[[2]],
                                                         Lifespan_centroid[[1]], Lifespan_centroid[[2]]) ## may take a while


sol_df_lifetimeegg_lifespan  <- integrate_geodesic(best_vT_lifetimeegg_lifespan , 
                                                   lifetimeegg_centroid[[1]], 
                                                   lifetimeegg_centroid[[2]])




### Calculating Area in between
### Area in between (Geo and Euclidean; Lifespan surface)
area_between_lifespan_geo <- area_geodesic_with_peak(
  surface_df  = Leedt_Survsurface_fixed_scaled,
  x_col_surf  = "PROT",
  y_col_surf  = "CARB",
  z_col_surf  = "tps",
  path_df     = sol_df_lifespan,
  standardize = TRUE,
  use_metric  = TRUE,
  signed_area = FALSE  # below-area uses |P|
)

area_between_lifespan_euc <- area_geodesic_with_peak(
  surface_df  = Leedt_Survsurface_fixed_scaled,
  x_col_surf  = "PROT",
  y_col_surf  = "CARB",
  z_col_surf  = "tps",
  path_df     = sol_df_lifespan,
  standardize = TRUE,
  use_metric  = FALSE,
  signed_area = FALSE  # below-area uses |P|
)

### Area in between (Geo and Euclidean; RR surface)
area_between_RR_geo <- area_geodesic_with_peak(
  surface_df  = Leedt_RRsurface_fixed_scaled,
  x_col_surf  = "PROT",
  y_col_surf  = "CARB",
  z_col_surf  = "tps",
  path_df     = sol_df_lifespan,
  standardize = TRUE,
  use_metric  = TRUE,
  signed_area = FALSE  # below-area uses |P|
)



area_between_RR_euc <- area_geodesic_with_peak(
  surface_df  = Leedt_RRsurface_fixed_scaled,
  x_col_surf  = "PROT",
  y_col_surf  = "CARB",
  z_col_surf  = "tps",
  path_df     = sol_df_lifespan,
  standardize = TRUE,
  use_metric  = FALSE,
  signed_area = FALSE  # below-area uses |P|
)





## 
### 

### Area in between (Geo and Euclidean; Lifetime surface)
area_between_lifetime_geo <- area_geodesic_with_peak(
  surface_df  = Leedt_Lifetimeeggsurface_fixed_scaled,
  x_col_surf  = "PROT",
  y_col_surf  = "CARB",
  z_col_surf  = "tps",
  path_df     = sol_df_lifetimeegg_lifespan,
  standardize = TRUE,
  use_metric  = TRUE,
  signed_area = FALSE  # below-area uses |P|
)

area_between_lifetime_euc <- area_geodesic_with_peak(
  surface_df  = Leedt_Lifetimeeggsurface_fixed_scaled,
  x_col_surf  = "PROT",
  y_col_surf  = "CARB",
  z_col_surf  = "tps",
  path_df     = sol_df_lifetimeegg_lifespan,
  standardize = TRUE,
  use_metric  = FALSE,
  signed_area = FALSE  # below-area uses |P|
)


### Area in between (Geo and Euclidean; RR surface)

data.frame(lifespan_geo = area_between_lifespan_geo$area_between_geo_peak,
           lifespan_euc =area_between_lifespan_euc$area_between_geo_peak,
           RR_geo = area_between_RR_geo$area_between_geo_peak,
           RR_euc =area_between_RR_euc$area_between_geo_peak,
           lifetime_geo = area_between_lifetime_geo$area_between_geo_peak,
           lifetime_euc =area_between_lifetime_euc$area_between_geo_peak)








### running multiple times
##
### Lifespan
### running multiple times
##
### Lifespan
N <- 50

a <- fit_lifespan$coefs[["a"]]
b <- fit_lifespan$coefs[["b"]]
c <- fit_lifespan$coefs[["c"]]
d <- fit_lifespan$coefs[["d"]]
e <- fit_lifespan$coefs[["e"]]


lifespan_peak_small <- predictonML_lm_lifespan_fixed$predicted_df %>% sample_n(N)
lifetime_peak_small <- predictonML_lm_lifetimeegg_fixed$predicted_df %>% sample_n(N)
rr_peak_small       <- predictonML_lm_dailyeggs_fixed$predicted_df %>% sample_n(N)


### Lifespan surface: Lifespan ↔ RR
all_paths_lifespan <- list()

for (i in 1:N) {
  x0 <- lifespan_peak_small$PROT[i]
  y0 <- lifespan_peak_small$CARB[i]
  x1 <- rr_peak_small$PROT[i]
  y1 <- rr_peak_small$CARB[i]
  
  cat(sprintf("Sample %d of %d...\n", i, N))
  
  path <- tryCatch({
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    df <- data.frame(integrate_geodesic(vT, x0, y0))
    df$sample <- i
    all_paths_lifespan[[i]] <- df
  }, error = function(e) NULL)
}

paths_df_lifespan <- bind_rows(all_paths_lifespan) %>%
  mutate(trait   = "Lifespan-RR",
         surface = "Lifespan",
         type    = "geodesic")



### RR surface: Lifespan ↔ RR
a <- fit_dailyeggs$coefs[["a"]]
b <- fit_dailyeggs$coefs[["b"]]
c <- fit_dailyeggs$coefs[["c"]]
d <- fit_dailyeggs$coefs[["d"]]
e <- fit_dailyeggs$coefs[["e"]]

all_paths_RR <- list()

for (i in 1:N) {
  x0 <- lifespan_peak_small$PROT[i]
  y0 <- lifespan_peak_small$CARB[i]
  x1 <- rr_peak_small$PROT[i]
  y1 <- rr_peak_small$CARB[i]
  
  cat(sprintf("Sample %d of %d...\n", i, N))
  
  path <- tryCatch({
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    df <- data.frame(integrate_geodesic(vT, x0, y0))
    df$sample <- i
    all_paths_RR[[i]] <- df
  }, error = function(e) NULL)
}

paths_df_RR <- bind_rows(all_paths_RR) %>%
  mutate(trait   = "Lifespan-RR",
         surface = "RR",
         type    = "geodesic")



### Lifetime egg surface – centroids (Lifespan vs Lifetime, RR vs Lifetime)
a <- fit_lifetimeegg$coefs[["a"]]
b <- fit_lifetimeegg$coefs[["b"]]
c <- fit_lifetimeegg$coefs[["c"]]
d <- fit_lifetimeegg$coefs[["d"]]
e <- fit_lifetimeegg$coefs[["e"]]

best_vT_lifetimeegg_lifespan <- robust_optimize_geodesic(
  lifetimeegg_centroid[[1]], lifetimeegg_centroid[[2]],
  Lifespan_centroid[[1]],    Lifespan_centroid[[2]]
) ## may take a while

sol_df_lifetimeegg_lifespan  <- integrate_geodesic(
  best_vT_lifetimeegg_lifespan , 
  lifetimeegg_centroid[[1]], 
  lifetimeegg_centroid[[2]]
)

best_vT_lifetimeegg_dailyeggs <- robust_optimize_geodesic(
  RR_centroid[[1]],          RR_centroid[[2]],
  lifetimeegg_centroid[[1]], lifetimeegg_centroid[[2]]
) ## may take a while

sol_df_lifetimeegg_dailyeggs  <- integrate_geodesic(
  best_vT_lifetimeegg_dailyeggs, 
  RR_centroid[[1]], 
  RR_centroid[[2]]
)



### Lifetime eggs surface: Lifetime ↔ Lifespan
all_paths_lifetime <- list()

for (i in 1:N) {
  x0 <- lifetime_peak_small$PROT[i]
  y0 <- lifetime_peak_small$CARB[i]
  x1 <- lifespan_peak_small$PROT[i]
  y1 <- lifespan_peak_small$CARB[i]
  
  cat(sprintf("Sample %d of %d...\n", i, N))
  
  path <- tryCatch({
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    df <- data.frame(integrate_geodesic(vT, x0, y0))
    df$sample <- i
    all_paths_lifetime[[i]] <- df
  }, error = function(e) NULL)
}

paths_df_lifetime <- bind_rows(all_paths_lifetime) %>%
  mutate(trait   = "Lifespan-Lifetime",
         surface = "Lifetime eggs",
         type    = "geodesic")



### Lifetime eggs surface: Lifetime ↔ RR
all_paths_lifetime_dailyeggs <- list()

for (i in 1:N) {
  x0 <- lifetime_peak_small$PROT[i]
  y0 <- lifetime_peak_small$CARB[i]
  x1 <- rr_peak_small$PROT[i]
  y1 <- rr_peak_small$CARB[i]
  
  cat(sprintf("Sample %d of %d...\n", i, N))
  
  path <- tryCatch({
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    df <- data.frame(integrate_geodesic(vT, x0, y0))
    df$sample <- i
    all_paths_lifetime_dailyeggs[[i]] <- df
  }, error = function(e) NULL)
}

paths_df_lifetime_dailyeggs <- bind_rows(all_paths_lifetime_dailyeggs ) %>%
  mutate(trait   = "Lifetime-RR",
         surface = "Lifetime eggs",
         type    = "geodesic")


## Lifetime-Lifespan with Lifespan surface
a <- fit_lifespan$coefs[["a"]]
b <- fit_lifespan$coefs[["b"]]
c <- fit_lifespan$coefs[["c"]]
d <- fit_lifespan$coefs[["d"]]
e <- fit_lifespan$coefs[["e"]]


all_paths_lifetime_lifespan_surflifespan <- list()

for (i in 1:N) {
  x0 <- lifetime_peak_small$PROT[i]
  y0 <- lifetime_peak_small$CARB[i]
  x1 <- lifespan_peak_small$PROT[i]
  y1 <- lifespan_peak_small$CARB[i]
  
  cat(sprintf("Sample %d of %d...\n", i, N))
  
  path <- tryCatch({
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    df <- data.frame(integrate_geodesic(vT, x0, y0))
    df$sample <- i
    all_paths_lifetime_lifespan_surflifespan[[i]] <- df
  }, error = function(e) NULL)
}

paths_df_lifetime_lifespan_surflifespan <- bind_rows(all_paths_lifetime_lifespan_surflifespan) %>%
  mutate(trait   = "Lifespan-Lifetime", 
         surface = "Lifespan",
         type    = "geodesic")



best_vT_lifetimeegg_lifespan_surflifespan <- robust_optimize_geodesic(
  lifetimeegg_centroid[[1]], lifetimeegg_centroid[[2]],
  Lifespan_centroid[[1]],    Lifespan_centroid[[2]]
) ## may take a while

sol_df_lifetimeegg_lifespan_surflifespan  <- integrate_geodesic(
  best_vT_lifetimeegg_lifespan_surflifespan, 
  lifetimeegg_centroid[[1]], 
  lifetimeegg_centroid[[2]]
)



### Lifetime-RR with RR surface
a <- fit_dailyeggs$coefs[["a"]]
b <- fit_dailyeggs$coefs[["b"]]
c <- fit_dailyeggs$coefs[["c"]]
d <- fit_dailyeggs$coefs[["d"]]
e <- fit_dailyeggs$coefs[["e"]]

all_paths_lifetime_dailyeggs_surfdailyeggs <- list()

for (i in 1:N) {
  x0 <- lifetime_peak_small$PROT[i]
  y0 <- lifetime_peak_small$CARB[i]
  x1 <- rr_peak_small$PROT[i]
  y1 <- rr_peak_small$CARB[i]
  
  cat(sprintf("Sample %d of %d...\n", i, N))
  
  path <- tryCatch({
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    df <- data.frame(integrate_geodesic(vT, x0, y0))
    df$sample <- i
    all_paths_lifetime_dailyeggs_surfdailyeggs[[i]] <- df
  }, error = function(e) NULL)
}

paths_df_lifetime_dailyeggs_surfdailyeggs <- bind_rows(all_paths_lifetime_dailyeggs_surfdailyeggs) %>%
  mutate(trait   = "Lifetime-RR", 
         surface = "RR",
         type    = "geodesic")



best_vT_lifetimeegg_dailyeggs_surfdailyeggs <- robust_optimize_geodesic(
  RR_centroid[[1]],          RR_centroid[[2]],
  lifetimeegg_centroid[[1]], lifetimeegg_centroid[[2]]
) ## may take a while

sol_df_lifetimeegg_dailyeggs_surfdailyeggs  <- integrate_geodesic(
  best_vT_lifetimeegg_dailyeggs_surfdailyeggs, 
  RR_centroid[[1]], 
  RR_centroid[[2]]
)



## Combine all paths
paths_df_all <- bind_rows(
  paths_df_lifespan,
  paths_df_RR,
  paths_df_lifetime, 
  paths_df_lifetime_dailyeggs,
  paths_df_lifetime_lifespan_surflifespan,
  paths_df_lifetime_dailyeggs_surfdailyeggs
)
unique(paths_df_all$surface)

## saving as dataframe so we dont need to run all the time
write.csv(paths_df_all, "Lee_et_al_Geodesic_paths_Lifespan_RR_Lifetime_N_50.csv")

#paths_df_all  <- read.csv("Lee_et_al_Geodesic_paths_Lifespan_RR_Lifetime_N_50.csv")



#### AREA / DISTANCE ESTIMATES
#### Now we compute areas/distances for EACH trait × surface combination
#### matching the paths exactly.

area_df_trait_surface <- c()

for (j in 1:N) {
  cat(sprintf("Area sample %d of %d...\n", j, N))
  
  ## 1) Lifespan-RR on Lifespan surface
  path_df_lr_Lifespan <- subset(paths_df_all,
                                trait   == "Lifespan-RR" &
                                  surface == "Lifespan" &
                                  sample  == j)
  if (nrow(path_df_lr_Lifespan) > 0) {
    inter_geo  <- area_geodesic_with_peak(
      surface_df  = Leedt_Survsurface_fixed_scaled,
      x_col_surf  = "PROT",
      y_col_surf  = "CARB",
      z_col_surf  = "tps",
      path_df     = path_df_lr_Lifespan,
      standardize = TRUE,
      use_metric  = TRUE,
      signed_area = FALSE
    )
    
    inter_euc  <- area_geodesic_with_peak(
      surface_df  = Leedt_Survsurface_fixed_scaled,
      x_col_surf  = "PROT",
      y_col_surf  = "CARB",
      z_col_surf  = "tps",
      path_df     = path_df_lr_Lifespan,
      standardize = TRUE,
      use_metric  = FALSE,
      signed_area = FALSE
    )
    
    area_df_trait_surface <- bind_rows(
      area_df_trait_surface,
      data.frame(area_below_geo   = inter_geo$area_below_geo_baseline,
                 area_between_geo = inter_geo$area_between_geo_peak,
                 path_length_geo  = inter_geo$path_length,
                 area_below_euc   = inter_euc$area_below_geo_baseline,
                 area_between_euc = inter_euc$area_between_geo_peak,
                 path_length_euc  = inter_euc$path_length,
                 surface          = "Lifespan",
                 trait            = "Lifespan-RR",
                 sample           = j)
    )
  }
  
  ## 2) Lifespan-RR on RR surface
  path_df_lr_RR <- subset(paths_df_all,
                          trait   == "Lifespan-RR" &
                            surface == "RR" &
                            sample  == j)
  if (nrow(path_df_lr_RR) > 0) {
    inter_geo  <- area_geodesic_with_peak(
      surface_df  = Leedt_RRsurface_fixed_scaled,
      x_col_surf  = "PROT",
      y_col_surf  = "CARB",
      z_col_surf  = "tps",
      path_df     = path_df_lr_RR,
      standardize = TRUE,
      use_metric  = TRUE,
      signed_area = FALSE
    )
    
    inter_euc  <- area_geodesic_with_peak(
      surface_df  = Leedt_RRsurface_fixed_scaled,
      x_col_surf  = "PROT",
      y_col_surf  = "CARB",
      z_col_surf  = "tps",
      path_df     = path_df_lr_RR,
      standardize = TRUE,
      use_metric  = FALSE,
      signed_area = FALSE
    )
    
    area_df_trait_surface <- bind_rows(
      area_df_trait_surface,
      data.frame(area_below_geo   = inter_geo$area_below_geo_baseline,
                 area_between_geo = inter_geo$area_between_geo_peak,
                 path_length_geo  = inter_geo$path_length,
                 area_below_euc   = inter_euc$area_below_geo_baseline,
                 area_between_euc = inter_euc$area_between_geo_peak,
                 path_length_euc  = inter_euc$path_length,
                 surface          = "RR",
                 trait            = "Lifespan-RR",
                 sample           = j)
    )
  }
  
  ## 3) Lifespan-Lifetime on Lifespan surface
  path_df_ll_Lifespan <- subset(paths_df_all,
                                trait   == "Lifespan-Lifetime" &
                                  surface == "Lifespan" &
                                  sample  == j)
  if (nrow(path_df_ll_Lifespan) > 0) {
    inter_geo  <- area_geodesic_with_peak(
      surface_df  = Leedt_Survsurface_fixed_scaled,
      x_col_surf  = "PROT",
      y_col_surf  = "CARB",
      z_col_surf  = "tps",
      path_df     = path_df_ll_Lifespan,
      standardize = TRUE,
      use_metric  = TRUE,
      signed_area = FALSE
    )
    
    inter_euc  <- area_geodesic_with_peak(
      surface_df  = Leedt_Survsurface_fixed_scaled,
      x_col_surf  = "PROT",
      y_col_surf  = "CARB",
      z_col_surf  = "tps",
      path_df     = path_df_ll_Lifespan,
      standardize = TRUE,
      use_metric  = FALSE,
      signed_area = FALSE
    )
    
    area_df_trait_surface <- bind_rows(
      area_df_trait_surface,
      data.frame(area_below_geo   = inter_geo$area_below_geo_baseline,
                 area_between_geo = inter_geo$area_between_geo_peak,
                 path_length_geo  = inter_geo$path_length,
                 area_below_euc   = inter_euc$area_below_geo_baseline,
                 area_between_euc = inter_euc$area_between_geo_peak,
                 path_length_euc  = inter_euc$path_length,
                 surface          = "Lifespan",
                 trait            = "Lifespan-Lifetime",
                 sample           = j)
    )
  }
  
  ## 4) Lifespan-Lifetime on Lifetime eggs surface
  path_df_ll_Lifetime <- subset(paths_df_all,
                                trait   == "Lifespan-Lifetime" &
                                  surface == "Lifetime eggs" &
                                  sample  == j)
  if (nrow(path_df_ll_Lifetime) > 0) {
    inter_geo  <- area_geodesic_with_peak(
      surface_df  = Leedt_Lifetimeeggsurface_fixed_scaled,
      x_col_surf  = "PROT",
      y_col_surf  = "CARB",
      z_col_surf  = "tps",
      path_df     = path_df_ll_Lifetime,
      standardize = TRUE,
      use_metric  = TRUE,
      signed_area = FALSE
    )
    
    inter_euc  <- area_geodesic_with_peak(
      surface_df  = Leedt_Lifetimeeggsurface_fixed_scaled,
      x_col_surf  = "PROT",
      y_col_surf  = "CARB",
      z_col_surf  = "tps",
      path_df     = path_df_ll_Lifetime,
      standardize = TRUE,
      use_metric  = FALSE,
      signed_area = FALSE
    )
    
    area_df_trait_surface <- bind_rows(
      area_df_trait_surface,
      data.frame(area_below_geo   = inter_geo$area_below_geo_baseline,
                 area_between_geo = inter_geo$area_between_geo_peak,
                 path_length_geo  = inter_geo$path_length,
                 area_below_euc   = inter_euc$area_below_geo_baseline,
                 area_between_euc = inter_euc$area_between_geo_peak,
                 path_length_euc  = inter_euc$path_length,
                 surface          = "Lifetime eggs",
                 trait            = "Lifespan-Lifetime",
                 sample           = j)
    )
  }
  
  ## 5) Lifetime-RR on Lifetime eggs surface
  path_df_lr2_Lifetime <- subset(paths_df_all,
                                 trait   == "Lifetime-RR" &
                                   surface == "Lifetime eggs" &
                                   sample  == j)
  if (nrow(path_df_lr2_Lifetime) > 0) {
    inter_geo  <- area_geodesic_with_peak(
      surface_df  = Leedt_Lifetimeeggsurface_fixed_scaled,
      x_col_surf  = "PROT",
      y_col_surf  = "CARB",
      z_col_surf  = "tps",
      path_df     = path_df_lr2_Lifetime,
      standardize = TRUE,
      use_metric  = TRUE,
      signed_area = FALSE
    )
    
    inter_euc  <- area_geodesic_with_peak(
      surface_df  = Leedt_Lifetimeeggsurface_fixed_scaled,
      x_col_surf  = "PROT",
      y_col_surf  = "CARB",
      z_col_surf  = "tps",
      path_df     = path_df_lr2_Lifetime,
      standardize = TRUE,
      use_metric  = FALSE,
      signed_area = FALSE
    )
    
    area_df_trait_surface <- bind_rows(
      area_df_trait_surface,
      data.frame(area_below_geo   = inter_geo$area_below_geo_baseline,
                 area_between_geo = inter_geo$area_between_geo_peak,
                 path_length_geo  = inter_geo$path_length,
                 area_below_euc   = inter_euc$area_below_geo_baseline,
                 area_between_euc = inter_euc$area_between_geo_peak,
                 path_length_euc  = inter_euc$path_length,
                 surface          = "Lifetime eggs",
                 trait            = "Lifetime-RR",
                 sample           = j)
    )
  }
  
  ## 6) Lifetime-RR on RR surface
  path_df_lr2_RR <- subset(paths_df_all,
                           trait   == "Lifetime-RR" &
                             surface == "RR" &
                             sample  == j)
  if (nrow(path_df_lr2_RR) > 0) {
    inter_geo  <- area_geodesic_with_peak(
      surface_df  = Leedt_RRsurface_fixed_scaled,
      x_col_surf  = "PROT",
      y_col_surf  = "CARB",
      z_col_surf  = "tps",
      path_df     = path_df_lr2_RR,
      standardize = TRUE,
      use_metric  = TRUE,
      signed_area = FALSE
    )
    
    inter_euc  <- area_geodesic_with_peak(
      surface_df  = Leedt_RRsurface_fixed_scaled,
      x_col_surf  = "PROT",
      y_col_surf  = "CARB",
      z_col_surf  = "tps",
      path_df     = path_df_lr2_RR,
      standardize = TRUE,
      use_metric  = FALSE,
      signed_area = FALSE
    )
    
    area_df_trait_surface <- bind_rows(
      area_df_trait_surface,
      data.frame(area_below_geo   = inter_geo$area_below_geo_baseline,
                 area_between_geo = inter_geo$area_between_geo_peak,
                 path_length_geo  = inter_geo$path_length,
                 area_below_euc   = inter_euc$area_below_geo_baseline,
                 area_between_euc = inter_euc$area_between_geo_peak,
                 path_length_euc  = inter_euc$path_length,
                 surface          = "RR",
                 trait            = "Lifetime-RR",
                 sample           = j)
    )
  }
}



## ------------------------------------------------------------------
## 1) "Per surface" dataset similar to your original area_df_lee_all
## ------------------------------------------------------------------

area_df_lee_all <- area_df_trait_surface %>%
  mutate(surface = ifelse(surface == "Lifespan", "lifespan",
                          ifelse(surface == "Lifetime eggs", "lifetime", "RR")))

write.csv(area_df_lee_all, "Lee_et_al_Geodesic_area_between_all_traits.csv")

### standardising
# 1. Compute means from leedt
mean_lifespan     <- mean(leedt$lifespan, na.rm = TRUE)
mean_lifetimeegg  <- mean(leedt$lifetimeegg, na.rm = TRUE)
mean_dailyeggs    <- mean(leedt$dailyeggs, na.rm = TRUE)

# 2. Create a lookup table telling which mean to use for each "surface" type
scale_lookup <- tibble(
  surface = c("lifespan", "lifetime", "RR"),
  scale   = c(mean_lifespan, mean_lifetimeegg, mean_dailyeggs)
)

# 3. Join the lookup table and scale numeric columns
area_df_scaled <- area_df_lee_all %>%
  left_join(scale_lookup, by = "surface") %>%
  mutate(across(
    .cols = where(is.numeric) & !c(sample),   # choose numeric columns except "sample"
    .fns = ~ .x / scale                        # divide by correct mean
  )) %>%
  select(-scale)  # remove helper column


### Wilcoxon tests geodesic vs Euclidean by surface
wilcox.test(
  area_df_scaled$area_between_geo[area_df_lee_all$surface == "lifespan"], 
  area_df_scaled$area_between_euc[area_df_lee_all$surface == "lifespan"],
  paired = TRUE
)

wilcox.test(
  area_df_scaled$area_between_geo[area_df_lee_all$surface == "RR"], 
  area_df_scaled$area_between_euc[area_df_lee_all$surface == "RR"],
  paired = TRUE
)

wilcox.test(
  area_df_scaled$area_between_geo[area_df_lee_all$surface == "lifetime"], 
  area_df_scaled$area_between_euc[area_df_lee_all$surface == "lifetime"],
  paired = TRUE
)



### plotting area between (per surface)
area_df_long_lee <- area_df_scaled %>%
  group_by(surface) %>%
  tidyr::pivot_longer(
    cols = c(area_between_geo, area_between_euc),
    names_to = "Type",
    values_to = "Area"
  ) %>%
  mutate(
    Type    = ifelse(Type == "area_between_geo", "Geodesic", "Euclidean"),
    surface = ifelse(surface == "RR",        "Surface:\nRepr. rate", 
                     ifelse(surface == "lifespan", "Surface:\nLifespan",
                            "Surface:\nLifetime eggs"))
  )

area_df_long_summary_lee <- area_df_long_lee %>%
  group_by(Type, surface) %>%
  summarise(
    avg = mean(Area),
    std = sd(Area) / sqrt(n()),
    .groups = "drop"
  )

## barplot distance and area between
##### Same traits, different landscapes (by surface)


## ------------------------------------------------------------------
## 2) TRAIT-RELATIVE COMPARISONS: directly from area_df_trait_surface
## ------------------------------------------------------------------

## Clean labels for trait-relative plotting
area_trait_relative_all <- area_df_trait_surface %>%
  mutate(surface = ifelse(surface == "Lifespan", "Lifespan", 
                          ifelse(surface == "RR", "Repr. rate", "Lifetime eggs")),
         trait   = dplyr::recode(trait,
                                 "Lifetime-RR"       = "Lifetime-Repr. rate",
                                 "Lifespan-RR"       = "Lifespan-Repr. rate",
                                 "Lifespan-Lifetime" = "Lifespan-Lifetime eggs"))


area_trait_relative_all_std <- area_trait_relative_all %>%
  mutate(
    across(
      c(area_below_geo:path_length_euc),
      ~ case_when(
        surface == "Lifespan"      ~ .x / mean_lifespan,
        surface == "Lifetime eggs" ~ .x / mean_lifetimeegg,
        surface == "Repr. rate"    ~ .x / mean_dailyeggs,
        TRUE                       ~ .x          # fallback (no change)
      )
    )
  )

### barplot for trait-relative comparisons (same trait, different surfaces)
area_trait_relative_all_long <- area_trait_relative_all_std %>%
  group_by(surface, trait) %>%
  tidyr::pivot_longer(
    cols = c(area_between_geo, area_between_euc),
    names_to = "Type",
    values_to = "Area"
  ) %>%
  mutate(
    Type = ifelse(Type == "area_between_geo", "Geodesic", "Euclidean")
  )

area_trait_relative_all_long_summary <- area_trait_relative_all_long %>%
  group_by(Type, trait, surface) %>%
  summarise(
    avg = mean(Area),
    std = sd(Area) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Comparison = ifelse(trait == "Lifespan-Lifetime eggs", "Lifespan\nvs\nLifetime eggs",
                             ifelse(trait == "Lifespan-Repr. rate", "Lifespan\nvs\nRepr. rate", "Lifetime eggs\nvs\nRepr. rate")))

## Example: table for paper (Table 1)
area_trait_relative_all_long_summary %>% 
  group_by(Type)




## graphs
costs_plot_surflifetime <- ggplot(subset(area_trait_relative_all_long_summary,  surface == "Lifetime eggs"), aes(x = Type, y = avg, fill = Type)) +
  facet_grid(~Comparison) + 
  geom_col(col = "black", width = 0.9) +
  geom_errorbar(aes(ymin = avg - std, ymax = avg + std), width = 0.1) +
  ylab(expression(italic(C[xy]))) + 
  theme_linedraw() + 
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5),
        axis.title  = element_text(size = 16),
        panel.grid  = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "lavenderblush2"),
        strip.text      = element_text(size = 13, color = "black")) + 
  scale_fill_manual('Type', values = c("lightskyblue2", "indianred2"))
costs_plot_surflifetime


### plot of surface across all landscapes
euclidean_path_lifespan_dailyeggs <- euclidean_path_2d(p1 = Lifespan_centroid,
                                                       p2 = RR_centroid)
euclidean_path_lifetime_dailyeggs <- euclidean_path_2d(p1 = lifetimeegg_centroid,
                                                       p2 = RR_centroid)
euclidean_path_lifetime_lifespan <- euclidean_path_2d(p1 = lifetimeegg_centroid,
                                                      p2 = Lifespan_centroid)




## ------------------------------------------------------------------
## Plots: same trait across different landscapes
## ------------------------------------------------------------------

### Lifespan vs Repr. rate – Lifespan surface
Survlandscape_plot_fixed_scaled_geo2 <- Survlandscape_plot_fixed_scaled +
  geom_point(data = RR_centroid,
             mapping = aes(x = PROT, y = CARB), 
             pch = 21, 
             fill = "firebrick2", 
             col  = "black",
             inherit.aes = FALSE,
             size = 3) +
  geom_path(data = sol_df_lifespan,           ## Geodesic
            mapping = aes(x = x, y = y),
            color = "indianred2", 
            linewidth = 1, 
            inherit.aes = FALSE) +
  geom_path(data = euclidean_path_lifespan_dailyeggs,  ## Euclidean
            mapping = aes(x = x, y = y), 
            color = "lightskyblue2",
            linetype = "dashed", 
            linewidth = 0.7,
            inherit.aes = FALSE) + 
  ggtitle("Lifespan\nvs\nRepr. rate", subtitle = "Landscape: Lifespan") +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title    = element_text(size = 14),
        axis.text     = element_text(size = 12))

Survlandscape_plot_fixed_scaled_geo2


### Lifespan vs Repr. rate – RR surface
Survlandscape_plot_fixed_scaled_geo3 <- RRlandscape_plot_fixed_scaled +
  geom_point(data = Lifespan_centroid,
             mapping = aes(x = PROT, y = CARB), 
             pch = 21, 
             fill = "steelblue2", 
             col  = "black",
             inherit.aes = FALSE,
             size = 3) +
  geom_path(data = sol_df_dailyeggs,          ## Geodesic
            mapping = aes(x = x, y = y),
            color = "indianred2", 
            linewidth = 1, 
            inherit.aes = FALSE) +
  geom_path(data = euclidean_path_lifespan_dailyeggs,  ## Euclidean
            mapping = aes(x = x, y = y), 
            color = "lightskyblue2",
            linetype = "dashed", 
            linewidth = 0.7,
            inherit.aes = FALSE) + 
  ggtitle("Lifespan\nvs\nRepr. rate", subtitle = "Landscape: Repr. rate") +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title    = element_text(size = 14),
        axis.text     = element_text(size = 12))

Survlandscape_plot_fixed_scaled_geo3



## ------------------------------------------------------------------
## Same landscape, different traits
## ------------------------------------------------------------------

### Lifetime eggs vs Lifespan – Lifetime eggs surface
Lifetimeggslandscape_plot_fixed_scaled_geo_1 <- Lifetimeegglandscape_plot_fixed_scaled +
  geom_point(data = Lifespan_centroid,
             mapping = aes(x = PROT, y = CARB), 
             pch = 21, 
             fill = "steelblue2", 
             col  = "black",
             inherit.aes = FALSE,
             size = 3) +
  geom_point(data = lifetimeegg_centroid,
             mapping = aes(x = PROT, y = CARB), 
             pch = 21, 
             fill = "orange", 
             col  = "black",
             inherit.aes = FALSE,
             size = 3) +
  geom_path(data = sol_df_lifetimeegg_lifespan,        ## Geodesic
            mapping = aes(x = x, y = y),
            color = "indianred2", 
            linewidth = 1, 
            inherit.aes = FALSE) +
  geom_path(data = euclidean_path_lifetime_lifespan,    ## Euclidean
            mapping = aes(x = x, y = y), 
            color = "lightskyblue2",
            linetype = "dashed", 
            linewidth = 0.7,
            inherit.aes = FALSE) + 
  ggtitle("Lifetime eggs\nvs\nLifespan", subtitle = "Landscape: Lifetime eggs") +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title    = element_text(size = 14),
        axis.text     = element_text(size = 12))

Lifetimeggslandscape_plot_fixed_scaled_geo_1


### Lifetime eggs vs Repr. rate – Lifetime eggs surface
Lifetimeggslandscape_plot_fixed_scaled_geo_2 <- Lifetimeegglandscape_plot_fixed_scaled +
  geom_point(data = RR_centroid,
             mapping = aes(x = PROT, y = CARB), 
             pch = 21, 
             fill = "firebrick2", 
             col  = "black",
             inherit.aes = FALSE,
             size = 3) +
  geom_path(data = sol_df_lifetimeegg_dailyeggs,        ## Geodesic
            mapping = aes(x = x, y = y),
            color = "indianred2", 
            linewidth = 1, 
            inherit.aes = FALSE) +
  geom_path(data = euclidean_path_lifetime_dailyeggs,   ## Euclidean
            mapping = aes(x = x, y = y), 
            color = "lightskyblue2",
            linetype = "dashed", 
            linewidth = 0.7,
            inherit.aes = FALSE) + 
  ggtitle("Lifetime eggs\nvs\nRepr. rate", subtitle = "Landscape: Lifetime eggs") +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title    = element_text(size = 14),
        axis.text     = element_text(size = 12))


### Lifetime eggs vs Repr. rate – RR surface
Lifetimeggslandscape_plot_fixed_scaled_geo_3 <- RRlandscape_plot_fixed_scaled +
  geom_point(data = lifetimeegg_centroid,
             mapping = aes(x = PROT, y = CARB), 
             pch = 21, 
             fill = "orange", 
             col  = "black",
             inherit.aes = FALSE,
             size = 3) +
  geom_point(data = RR_centroid,
             mapping = aes(x = PROT, y = CARB), 
             pch = 21, 
             fill = "firebrick2", 
             col  = "black",
             inherit.aes = FALSE,
             size = 3) +
  geom_path(data = sol_df_lifetimeegg_dailyeggs_surfdailyeggs,  ## Geodesic
            mapping = aes(x = x, y = y),
            color = "indianred2", 
            linewidth = 1, 
            inherit.aes = FALSE) +
  geom_path(data = euclidean_path_lifetime_dailyeggs,           ## Euclidean
            mapping = aes(x = x, y = y), 
            color = "lightskyblue2",
            linetype = "dashed", 
            linewidth = 0.7,
            inherit.aes = FALSE) + 
  ggtitle("Lifetime eggs\nvs\nRepr. rate", subtitle = "Landscape: Repr. rate") +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title    = element_text(size = 14),
        axis.text     = element_text(size = 12))

Lifetimeggslandscape_plot_fixed_scaled_geo_3


### Lifetime eggs vs Lifespan – Lifespan surface
Lifetimeggslandscape_plot_fixed_scaled_geo_4 <- Survlandscape_plot_fixed_scaled +
  geom_point(data = lifetimeegg_centroid,
             mapping = aes(x = PROT, y = CARB), 
             pch = 21, 
             fill = "orange", 
             col  = "black",
             inherit.aes = FALSE,
             size = 3) +
  geom_path(data = sol_df_lifetimeegg_lifespan_surflifespan,   ## Geodesic
            mapping = aes(x = x, y = y),
            color = "indianred2", 
            linewidth = 1, 
            inherit.aes = FALSE) +
  geom_path(data = euclidean_path_lifetime_lifespan,           ## Euclidean
            mapping = aes(x = x, y = y), 
            color = "lightskyblue2",
            linetype = "dashed", 
            linewidth = 0.7,
            inherit.aes = FALSE) + 
  ggtitle("Lifetime eggs\nvs\nLifespan", subtitle = "Landscape: Lifespan") +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title    = element_text(size = 14),
        axis.text     = element_text(size = 12))

Lifetimeggslandscape_plot_fixed_scaled_geo_4


(Survlandscape_plot_fixed_scaled_geo2 +
    Survlandscape_plot_fixed_scaled_geo3 +
    Lifetimeggslandscape_plot_fixed_scaled_geo_1) /
  (Lifetimeggslandscape_plot_fixed_scaled_geo_4 +
     Lifetimeggslandscape_plot_fixed_scaled_geo_3 +
     Lifetimeggslandscape_plot_fixed_scaled_geo_2)



