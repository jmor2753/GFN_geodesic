
# --- Load packages ---
library(ggplot2)
library(dplyr)
library(optimx)
library(deSolve)
library(patchwork)
library(broom)

# --- Load geodesic functions ---
source('...Geodesic_consolidated_Functions_6_December_2025.R')

set.seed(321908)

############################################################
#### 1. Parameters and landscape design #####################
############################################################

# --- Nutrient space ---
x_vals <- seq(0, 10, by = 0.05)
y_vals <- seq(0, 10, by = 0.05)

# --- Centre of the landscape ---
xc <- 5
yc <- 5

# --- Minimal curvature grid ---
curv_values <- c(-0.20, -0.05, 0, 0.05, 0.20)

landscape_grid <- expand.grid(
  A = curv_values,
  B = curv_values
) %>%
  mutate(
    landscape_id = row_number(),
    topology = case_when(
      A == 0 & B == 0 ~ "Flat",
      A < 0 & B < 0 ~ "Dome",
      A > 0 & B > 0 ~ "Bowl",
      A * B < 0 ~ "Saddle",
      A == 0 | B == 0 ~ "One-dimensional curvature"
    ),
    curvature_strength = abs(A) + abs(B)
  )

# --- Endpoint scenarios ---
endpoint_scenarios <- data.frame(
  scenario = c("Close points", "Far points"),
  x0 = c(4.5, 2),
  y0 = c(4.5, 2),
  x1 = c(5.5, 8),
  y1 = c(5.5, 8)
)

# --- Number of numerical repeats ---
n_repeats <- 10


############################################################
#### 2. Helper functions ####################################
############################################################

# Converts centred quadratic:
# z = A(x - xc)^2 + B(y - yc)^2
# into your existing form:
# z = ax^2 + by^2 + cx + dy + exy
#
# Constant term is omitted because it does not affect gradients,
# metric tensor, geodesics, or distance calculations.

set_landscape_parameters <- function(A, B, xc = 5, yc = 5) {
  
  assign("a", A, envir = .GlobalEnv)
  assign("b", B, envir = .GlobalEnv)
  assign("c", -2 * A * xc, envir = .GlobalEnv)
  assign("d", -2 * B * yc, envir = .GlobalEnv)
  assign("e", 0, envir = .GlobalEnv)
  
}


make_surface_grid <- function(A, B, xc = 5, yc = 5) {
  
  set_landscape_parameters(A, B, xc, yc)
  
  grid <- expand.grid(x = x_vals, y = y_vals)
  grid$P <- with(grid, P(x, y))
  
  return(grid)
}


make_euclidean_path <- function(x0, y0, x1, y1, n = 300) {
  data.frame(
    x = seq(x0, x1, length.out = n),
    y = seq(y0, y1, length.out = n)
  )
}


euclidean_distance <- function(x0, y0, x1, y1) {
  sqrt((x1 - x0)^2 + (y1 - y0)^2)
}


# --- Optimizer with jittered initial conditions ---
# This is only for numerical stability checks.
# The geodesic itself is deterministic for fixed surface and endpoints.

robust_optimize_geodesic_jitter <- function(x0, y0, x1, y1,
                                            jitter_sd = 0.05,
                                            lower = -20,
                                            upper = 20,
                                            tol = 1e-6) {
  
  guess <- generate_initial_guess(x0, y0, x1, y1)
  guess <- guess + rnorm(length(guess), mean = 0, sd = jitter_sd)
  guess[3] <- abs(guess[3])
  
  d_xy <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
  
  T_lower <- max(1e-4, 0.1 * d_xy)
  T_upper <- min(50, max(1, 10 * d_xy + 1))
  
  res <- tryCatch({
    optimx::optimx(
      par     = guess,
      fn      = function(vT) geodesic_objective(vT, x0, y0, x1, y1),
      method  = "L-BFGS-B",
      lower   = c(lower, lower, T_lower),
      upper   = c(upper, upper, T_upper),
      control = list(maxit = 100)
    )
  }, error = function(e) NULL)
  
  if (is.null(res)) return(NULL)
  if (any(is.na(res[1, c("p1", "p2", "p3")]))) return(NULL)
  
  vT <- as.numeric(res[1, c("p1", "p2", "p3")])
  err <- as.numeric(res[1, "value"])
  
  if (!is.finite(err)) return(NULL)
  
  if (err > tol) {
    warning(sprintf("Endpoint mismatch %.3e larger than tolerance %.3e", err, tol))
  }
  
  return(vT)
}


run_single_sensitivity <- function(A, B, topology,
                                   scenario, x0, y0, x1, y1,
                                   repeat_id = 1) {
  
  set_landscape_parameters(A, B, xc, yc)
  
  vT <- tryCatch({
    robust_optimize_geodesic_jitter(x0, y0, x1, y1)
  }, error = function(e) NULL)
  
  if (is.null(vT)) {
    return(data.frame(
      A = A, B = B,
      topology = topology,
      scenario = scenario,
      repeat_id = repeat_id,
      geo_dist = NA,
      euc_dist = euclidean_distance(x0, y0, x1, y1),
      ratio = NA,
      endpoint_error = NA,
      curvature_strength = abs(A) + abs(B)
    ))
  }
  
  sol_df <- integrate_geodesic(vT, x0, y0)
  
  endpoint_error <- sqrt(
    (tail(sol_df$x, 1) - x1)^2 +
      (tail(sol_df$y, 1) - y1)^2
  )
  
  geo_dist <- geodesic_distance(sol_df)
  euc_dist <- euclidean_distance(x0, y0, x1, y1)
  
  data.frame(
    A = A,
    B = B,
    topology = topology,
    scenario = scenario,
    repeat_id = repeat_id,
    geo_dist = geo_dist,
    euc_dist = euc_dist,
    ratio = geo_dist / euc_dist,
    endpoint_error = endpoint_error,
    curvature_strength = abs(A) + abs(B)
  )
}


############################################################
#### 3. Run sensitivity analysis ############################
############################################################

sensitivity_results <- list()
counter <- 1

for (i in 1:nrow(landscape_grid)) {
  
  for (j in 1:nrow(endpoint_scenarios)) {
    
    for (r in 1:n_repeats) {
      
      cat("Running landscape", i,
          "scenario", endpoint_scenarios$scenario[j],
          "repeat", r, "\n")
      
      sensitivity_results[[counter]] <- run_single_sensitivity(
        A = landscape_grid$A[i],
        B = landscape_grid$B[i],
        topology = landscape_grid$topology[i],
        scenario = endpoint_scenarios$scenario[j],
        x0 = endpoint_scenarios$x0[j],
        y0 = endpoint_scenarios$y0[j],
        x1 = endpoint_scenarios$x1[j],
        y1 = endpoint_scenarios$y1[j],
        repeat_id = r
      )
      
      counter <- counter + 1
    }
  }
}

sensitivity_results <- bind_rows(sensitivity_results)



############################################################
#### 4. Summarise results ###################################
############################################################

sensitivity_summary <- sensitivity_results %>%
  filter(!is.na(ratio)) %>%
  group_by(A, B, topology, scenario, curvature_strength) %>%
  summarise(
    geo_mean = mean(geo_dist),
    geo_sd = sd(geo_dist),
    euc_mean = mean(euc_dist),
    ratio_mean = mean(ratio),
    ratio_sd = sd(ratio),
    endpoint_error_mean = mean(endpoint_error),
    endpoint_error_max = max(endpoint_error),
    .groups = "drop"
  )


# Save raw results
saveRDS(sensitivity_results, file = "/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Nutritional Geometry/Nutrigonometry/23. Geodesics/Final Scripts/Final codes/Revisions 1/sensitivity_geodesic_euclidean_results.rds")




write.csv(
  sensitivity_summary,
  file = "/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Nutritional Geometry/Nutrigonometry/23. Geodesics/Final Scripts/Final codes/Revisions 1/sensitivity_geodesic_euclidean_summary.csv",
  row.names = FALSE
)


############################################################
#### 5. Example landscapes for figure panel A ###############
############################################################

example_landscapes <- data.frame(
  example = c("1D curvature", "Dome", "Bowl", "Saddle"),
  A = c(0.20, -0.20, 0.20, -0.20),
  B = c(0,    -0.20, 0.20,  0.20)
)

example_surface_list <- list()

for (i in 1:nrow(example_landscapes)) {
  
  tmp <- make_surface_grid(
    A = example_landscapes$A[i],
    B = example_landscapes$B[i],
    xc = xc,
    yc = yc
  )
  
  tmp$example <- example_landscapes$example[i]
  tmp$A <- example_landscapes$A[i]
  tmp$B <- example_landscapes$B[i]
  
  example_surface_list[[i]] <- tmp
}

example_surfaces <- bind_rows(example_surface_list)


############################################################
#### Add geodesic paths for close and far scenarios #########
############################################################

example_paths <- list()
counter <- 1

for (i in 1:nrow(example_landscapes)) {
  
  set_landscape_parameters(
    A = example_landscapes$A[i],
    B = example_landscapes$B[i],
    xc = xc,
    yc = yc
  )
  
  for (j in 1:nrow(endpoint_scenarios)) {
    
    x0 <- endpoint_scenarios$x0[j]
    y0 <- endpoint_scenarios$y0[j]
    x1 <- endpoint_scenarios$x1[j]
    y1 <- endpoint_scenarios$y1[j]
    
    scenario_now <- endpoint_scenarios$scenario[j]
    
    vT <- robust_optimize_geodesic(x0, y0, x1, y1)
    sol_df <- integrate_geodesic(vT, x0, y0)
    
    sol_df$type <- "Geodesic"
    sol_df$example <- example_landscapes$example[i]
    sol_df$scenario <- scenario_now
    
    example_paths[[counter]] <- sol_df[, c("x", "y", "type", "example", "scenario")]
    counter <- counter + 1
  }
}

example_paths <- bind_rows(example_paths)

example_paths$scenario <- factor(
  example_paths$scenario,
  levels = c("Close points", "Far points")
)

example_paths$example <- factor(
  example_paths$example,
  levels = c("1D curvature", "Dome", "Bowl", "Saddle")
)





# --- Panel A: examples ---
############################################################
#### Endpoint points for close and far scenarios ############
############################################################

endpoint_points <- endpoint_scenarios %>%
  tidyr::pivot_longer(
    cols = c(x0, y0, x1, y1),
    names_to = "coord",
    values_to = "value"
  ) %>%
  mutate(
    point = ifelse(coord %in% c("x0", "y0"), "Initial", "Final"),
    axis  = ifelse(coord %in% c("x0", "x1"), "x", "y")
  ) %>%
  select(scenario, point, axis, value) %>%
  tidyr::pivot_wider(
    names_from = axis,
    values_from = value
  ) %>%
  tidyr::crossing(
    example = unique(example_surfaces$example)
  )

endpoint_points$scenario <- factor(
  endpoint_points$scenario,
  levels = c("Close points", "Far points")
)



panel_A <- ggplot() +
  geom_tile(
    data = example_surfaces,
    aes(x = x, y = y, fill = P),
    colour = NA
  ) +
  geom_contour(
    data = example_surfaces,
    aes(x = x, y = y, z = P),
    colour = "grey40",
    linewidth = 0.15
  ) +
  geom_path(
    data = example_paths,
    aes(x = x, y = y, colour = scenario),
    linewidth = 0.8, alpha = 1
  ) +
  geom_point(
    data = endpoint_points,
    aes(x = x, y = y, colour = scenario),
    size = 2.5,
    alpha = 1
  ) +
  facet_wrap(~ example, nrow = 1) +
  scale_colour_manual(
    values = c(
      "Close points" = "royalblue2",
      "Far points" = "tomato2"
    )
  ) +
  scale_fill_viridis_c() +
  coord_equal() +
  labs(
    x = "Nutrient 1",
    y = "Nutrient 2",
    fill = "Performance",
    colour = "Endpoint scenario"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(size = 13),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),  # Increase label size
    legend.title = element_text(size = 12),  # Increase title size
    legend.position = "bottom"
  ) 
panel_A



# --- Panel B: sensitivity by curvature strength ---

############################################################
#### Add curvature per dimension ###########################
############################################################

sensitivity_summary2 <- sensitivity_summary %>%
  mutate(
    n_curved_dimensions = (A != 0) + (B != 0),
    curvature_per_dimension = ifelse(
      n_curved_dimensions == 0,
      0,
      (abs(A) + abs(B)) / n_curved_dimensions
    )
  )


############################################################
#### Duplicate Flat into all topology panels ###############
############################################################

topology_levels <- c(
  "One-dimensional curvature",
  "Dome",
  "Bowl",
  "Saddle"
)

flat_expanded <- sensitivity_summary2 %>%
  filter(A == 0, B == 0) %>%
  tidyr::crossing(topology_plot = topology_levels) %>%
  mutate(topology = topology_plot) %>%
  select(-topology_plot)


sensitivity_summary_plot <- sensitivity_summary2 %>%
  filter(!(A == 0 & B == 0)) %>%
  bind_rows(flat_expanded) %>%
  mutate(
    topology = factor(topology, levels = topology_levels),
    scenario = factor(scenario, levels = c("Close points", "Far points"))
  )


############################################################
#### Replicate sensitivity plot ############################
############################################################
sensitivity_summary_plot <- sensitivity_summary_plot %>%
  mutate(
    topology = recode(
      as.character(topology),
      "One-dimensional curvature" = "1D curvature"
    )
  )


## plot

panel_sensitivity_line <- ggplot(
  sensitivity_summary_plot,
  aes(
    x = curvature_per_dimension,
    y = ratio_mean,
    colour = scenario,
    group = scenario
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(size = 3, alpha = 0.9, pch = 18) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_errorbar(
    aes(
      ymin = ratio_mean - ratio_sd,
      ymax = ratio_mean + ratio_sd
    ),
    width = 0.005,
    linewidth = 0.3,
    alpha = 0.7
  ) +
  facet_grid(~ topology) +
  scale_colour_manual(values = c("royalblue2", "tomato2")) +
  scale_x_continuous(
    breaks = c(0, 0.05, 0.10, 0.15, 0.20),
    limits = c(0, 0.20)
  ) +
  labs(
    x = "Curvature magnitude per curved dimension",
    y = "Geodesic / Euclidean distance",
    colour = "Landscape topology",
    shape = "Landscape topology"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(size = 13),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),  # Increase label size
    legend.title = element_text(size = 12),  # Increase title size
    legend.position = "bottom", aspect.ratio = 1
  ) 

panel_sensitivity_line



# --- Panel C: When to use geodesics ---



# --- Panel C: heatmap across A and B values ---
panel_C <- ggplot(
  sensitivity_summary,
  aes(x = factor(A), y = factor(B), fill = ratio_mean)
) +
  geom_tile(color = "white") +
  facet_wrap(~ scenario) +
  scale_fill_viridis_c(name = "Geo / Euc") +
  labs(
    x = "Curvature in nutrient 1 (A)",
    y = "Curvature in nutrient 2 (B)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(size = 13),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),  # Increase label size
    legend.title = element_text(size = 12),  # Increase title size
    legend.position = "bottom", aspect.ratio = 1
  ) 

panel_C


panel_A / panel_sensitivity_line / panel_C



save.image(
  file = "Sensitivity_Geodesic_May2026_Revision1.RData"
)
