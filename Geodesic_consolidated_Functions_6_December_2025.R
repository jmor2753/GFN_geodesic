#### All functions Geodesic updated. - 06 Dec 2025


#### ~~~~ Scirpt to calculate Geodesics in performance landscapes ~~~~~ ######
# --- PACKAGES ---
library(deSolve)
library(ggplot2)
library(dplyr)
library(optimx)


#### -- functions
# --- SCALAR FIELD ---
P <- function(x, y) {
  a * x^2 + b * y^2 + c * x + d * y + e * x * y
}


### functions
P_grad <- function(x, y) {
  c(2 * a * x + c + e * y,
    2 * b * y + d + e * x)
}

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

# --- AUTO-GENERATE INITIAL GUESS ---
generate_initial_guess <- function(x0, y0, x1, y1, min_norm = 1e-4, default_dir = c(0.5, 0)) {
  dir <- c(x1 - x0, y1 - y0)
  norm <- sqrt(sum(dir^2))
  
  # Fallback if norm is too small
  if (norm < min_norm) {
    warning("Start and end points are too close. Using default direction for initial guess.")
    dir <- default_dir
    norm <- 1
  }
  
  v_unit <- dir / norm
  T0 <- max(norm * 2, 1)  # Minimum total time enforced
  return(c(v_unit * T0, T0))
} ### using this function instead of the other because this appears to accomodate more cases instead of throwing errors, and is equally inacurate.




# --- GEODESIC EQUATION ---
geodesic_eqs <- function(t, state, parameters) {
  x <- state[1]; dx <- state[2]
  y <- state[3]; dy <- state[4]
  
  Px <- 2 * a * x + c + e * y
  Py <- 2 * b * y + d + e * x
  
  g <- metric_tensor(x, y)
  speed2 <- dx * (g[1,1] * dx + g[1,2] * dy) + dy * (g[2,1] * dx + g[2,2] * dy)
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
  
  return(list(c(dx, ddx, dy, ddy)))
}




# --- INTEGRATE GEODESIC WITH PROGRESS BAR ---
## Shared objective: same as your current geodesic_objective
geodesic_objective <- function(vT, x0, y0, x_target, y_target) {
  df <- tryCatch({
    integrate_geodesic(vT, x0, y0)
  }, error = function(e) return(NULL))
  
  if (is.null(df) || nrow(df) == 0L) return(Inf)
  
  endpoint <- df[nrow(df), c("x", "y")]
  sum((endpoint - c(x_target, y_target))^2)
}


robust_optimize_geodesic <- function(x0, y0, x1, y1,
                                     lower = -20, upper = 20,
                                     tol = 1e-6) {
  guess <- generate_initial_guess(x0, y0, x1, y1)
  
  d <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
  
  # Version 1 bounds (fixed T)
  fixed_lower <- c(-20, -20, 1)
  fixed_upper <- c( 20,  20, 50)
  
  # Version 2 bounds (distance-scaled T, capped at 50)
  T_lower <- max(1e-4, 0.1 * d)
  T_upper <- min(50, max(1, 10 * d + 1))
  dyn_lower <- c(lower, lower, T_lower)
  dyn_upper <- c(upper, upper, T_upper)
  
  run_optim <- function(lower_vec, upper_vec) {
    res <- tryCatch({
      optimx::optimx(
        par     = guess,
        fn      = function(vT) geodesic_objective(vT, x0, y0, x1, y1),
        method  = "L-BFGS-B",
        lower   = lower_vec,
        upper   = upper_vec,
        control = list(maxit = 100)
      )
    }, error = function(e) NULL)
    
    if (is.null(res)) return(NULL)
    if (any(is.na(res[1, c("p1", "p2", "p3")]))) return(NULL)
    
    vT  <- as.numeric(res[1, c("p1", "p2", "p3")])
    err <- as.numeric(res[1, "value"])   # use optimx objective value
    
    if (!is.finite(err)) return(NULL)
    
    list(vT = vT, err = err)
  }
  
  cand_fixed <- run_optim(fixed_lower, fixed_upper)
  cand_dyn   <- run_optim(dyn_lower,   dyn_upper)
  
  if (is.null(cand_fixed) && is.null(cand_dyn)) {
    stop("Optimization failed for both fixed and distance-scaled bounds. ",
         "Try adjusting the initial guess or bounds.")
  }
  
  chosen <- NULL
  if (!is.null(cand_fixed) && is.null(cand_dyn)) {
    chosen <- cand_fixed
  } else if (is.null(cand_fixed) && !is.null(cand_dyn)) {
    chosen <- cand_dyn
  } else {
    chosen <- if (cand_fixed$err <= cand_dyn$err) cand_fixed else cand_dyn
  }
  
  if (chosen$err > tol) {
    warning(sprintf(
      "Geodesic endpoint mismatch is %.3e (> tol = %.3e). ",
      chosen$err, tol
    ))
  }
  
  return(chosen$vT)
}



## Your distance function can stay exactly as in Version 2
geodesic_distance <- function(df) {
  total <- 0
  for (i in 2:nrow(df)) {
    dx_vec <- c(df$x[i] - df$x[i-1], df$y[i] - df$y[i-1])
    mid_x  <- (df$x[i] + df$x[i-1]) / 2
    mid_y  <- (df$y[i] + df$y[i-1]) / 2
    g      <- metric_tensor(mid_x, mid_y)
    ds2    <- t(dx_vec) %*% g %*% dx_vec
    total  <- total + sqrt(ds2)
  }
  return(as.numeric(total))
}

integrate_geodesic <- function(vT, x0, y0) {
  v0 <- vT[1:2]
  T  <- abs(vT[3])  # Ensure T is positive
  
  dir_metric <- normalize_direction_signed(x0, y0, v0)
  state0     <- c(x = x0, dx = dir_metric[1], y = y0, dy = dir_metric[2])
  
  # Ensure we always have at least two time points (0 and T)
  if (T == 0) {
    T <- 1e-4  # or some small positive fallback
  }
  
  # Use a maximum step of 0.001, but shrink it if T is smaller
  dt    <- min(0.001, T / 10)  # e.g. ~10 steps even for tiny T
  times <- seq(0, T, by = dt)
  if (length(times) < 2L) {
    times <- c(0, T)
  }
  
  sol <- ode(
    y      = state0,
    times  = times,
    func   = geodesic_eqs,
    parms  = NULL,
    method = "lsoda",
    atol   = 1e-8,
    rtol   = 1e-8
  )
  
  as.data.frame(sol)
}



# --- Fit quadratic surface ---
fit_quadratic_surface <- function(data, x_col, y_col, z_col) {
  formula_str <- paste0(
    z_col, " ~ ",
    x_col, " + I(", x_col, "^2) + ",
    y_col, " + I(", y_col, "^2) + ",
    x_col, "*", y_col
  )
  model <- lm(as.formula(formula_str), data = data)
  coefs <- coef(model)
  
  a <- coefs[[paste0("I(", x_col, "^2)")]]
  b <- coefs[[paste0("I(", y_col, "^2)")]]
  c <- coefs[[x_col]]
  d <- coefs[[y_col]]
  e <- coefs[[paste0(x_col, ":", y_col)]]
  
  return(list(model = model, coefs = c(a = a, b = b, c = c, d = d, e = e)))
}


# ----------------------------------------------------------
# Function: test_geodesic_vs_euclidean
# ----------------------------------------------------------
# Compares geodesic and Euclidean distances using:
# - Permutation test of paired differences
# - Optional bootstrap of the mean difference
# - Returns bootstrap distribution for plotting
# ----------------------------------------------------------

test_geodesic_vs_euclidean <- function(geo, euc, n_perm = 10000, n_boot = 10000, do_bootstrap = TRUE, store_distributions = FALSE) {
  if(length(geo) != length(euc)) {
    stop("Geodesic and Euclidean vectors must be the same length.")
  }
  
  diff_vec <- geo - euc
  observed_diff <- mean(diff_vec)
  
  # --- Permutation test ---
  perm_diffs <- replicate(n_perm, {
    signs <- sample(c(1, -1), length(diff_vec), replace = TRUE)
    mean(signs * diff_vec)
  })
  
  p_val <- mean(abs(perm_diffs) >= abs(observed_diff))
  
  result <- list(
    observed_difference = observed_diff,
    permutation_p_value = p_val,
    bootstrap_mean = NULL,
    bootstrap_CI = NULL,
    bootstrap_distribution = NULL
  )
  
  # --- Optional bootstrap ---
  if (do_bootstrap) {
    boot_diffs <- replicate(n_boot, {
      sample_diffs <- sample(diff_vec, replace = TRUE)
      mean(sample_diffs)
    })
    
    ci <- quantile(boot_diffs, probs = c(0.025, 0.975))
    
    result$bootstrap_mean <- mean(boot_diffs)
    result$bootstrap_CI <- ci
    if (store_distributions){
      result$bootstrap_distribution <- boot_diffs
    } else {
      result$bootstrap_distribution <- NULL
    }
    
  }
  
  return(result)
}






####

# ---- AREA BELOW & ABOVE GEODESIC RELATIVE TO A PEAK CAP ----
# surface_df  : surface grid (x, y, z)
# path_df     : geodesic path (x, y)
# peak_df     : optional dataframe with x, y coordinates of peak region
#               (used only to define a cap height via max P over peak_df)
# x_col_path  : name of x column in path_df
# y_col_path  : name of y column in path_df
# x_col_surf  : name of x column in surface_df
# y_col_surf  : name of y column in surface_df
# z_col_surf  : name of z column in surface_df
# x_col_peak  : name of x column in peak_df (if provided)
# y_col_peak  : name of y column in peak_df (if provided)
# standardize : if TRUE, z is rescaled to [0,1] using surface_df
# use_metric  : if TRUE, ds uses metric_tensor; otherwise Euclidean
# signed_area : if TRUE, "below" area integrates signed z;
#               if FALSE, integrates |z| (total magnitude below baseline)
# cap_value   : optional peak height on the ORIGINAL P scale (e.g. 10);
#               if NULL, it is inferred from peak_df or surface_df



















area_geodesic_with_peak <- function(surface_df,
                                    path_df,
                                    peak_df     = NULL,
                                    x_col_path  = "x",
                                    y_col_path  = "y",
                                    x_col_surf  = "x",
                                    y_col_surf  = "y",
                                    z_col_surf  = "P",
                                    x_col_peak  = "x",
                                    y_col_peak  = "y",
                                    standardize = FALSE,
                                    use_metric  = FALSE,
                                    signed_area = TRUE,
                                    cap_value   = NULL) {
  
  
  
  
  ## --- 0. Basic checks ---
  if (!all(c(x_col_path, y_col_path) %in% names(path_df))) {
    stop("path_df must contain columns '", x_col_path,
         "' and '", y_col_path, "'.")
  }
  if (!all(c(x_col_surf, y_col_surf, z_col_surf) %in% names(surface_df))) {
    stop("surface_df must contain columns '", x_col_surf, "', '",
         y_col_surf, "', and '", z_col_surf, "'.")
  }
  if (!is.null(peak_df) &&
      !all(c(x_col_peak, y_col_peak) %in% names(peak_df))) {
    stop("peak_df must contain columns '", x_col_peak, "' and '",
         y_col_peak, "'.")
  }
  
  ## --- 1. Extract path coordinates ---
  x_path <- path_df[[x_col_path]]
  y_path <- path_df[[y_col_path]]
  
  ## --- 1b. Fit quadratic surface and define internal P(x, y) ---
  quad_fit <- fit_quadratic_surface(
    data  = surface_df,
    x_col = x_col_surf,
    y_col = y_col_surf,
    z_col = z_col_surf
  )
  coefs <- quad_fit$coefs
  
  a <- as.numeric(coefs["a"])
  b <- as.numeric(coefs["b"])
  c <- as.numeric(coefs["c"])
  d <- as.numeric(coefs["d"])
  e <- as.numeric(coefs["e"])
  
  # Internal, self-contained quadratic surface
  P_local <- function(x, y) {
    a * x^2 + b * y^2 + c * x + d * y + e * x * y
  }
  
  ## internal gradient
  P_grad_local <- function(x, y, coefs = quad_fit$coefs) {
    a <- coefs["a"]; b <- coefs["b"]; c <- coefs["c"]; d <- coefs["d"]; e <- coefs["e"]
    c(2 * a * x + c + e * y,
      2 * b * y + d + e * x)
  }
  
  metric_tensor_local <- function(x, y) {
    grad <- P_grad_local(x, y)
    Px <- grad[1]; Py <- grad[2]
    matrix(c(1 + Px^2, Px * Py, Px * Py, 1 + Py^2), nrow = 2)
  }
  
  ## --- 2. Evaluate P along the path using the internal quadratic ---
  z_raw <- P_local(x_path, y_path)
  
  ## --- 3. Possibly standardize z to [0,1] based on the whole surface ---
  z_surf_raw <- surface_df[[z_col_surf]]
  z_min <- min(z_surf_raw, na.rm = TRUE)
  z_max <- max(z_surf_raw, na.rm = TRUE)
  
  if (standardize) {
    if (z_max == z_min) {
      warning("Surface z-range is zero; standardized values set to 0.")
      z <- rep(0, length(z_raw))
    } else {
      z <- (z_raw - z_min) / (z_max - z_min)
    }
  } else {
    z <- z_raw
  }
  
  ## --- 4. Define "below" integrand: signed or absolute ---
  if (signed_area) {
    z_below_int <- z          # signed area below baseline
  } else {
    z_below_int <- abs(z)     # total magnitude below baseline
  }
  
  ## --- 5. Compute ds for each segment along the path ---
  dx <- diff(x_path)
  dy <- diff(y_path)
  
  if (use_metric) {
    if (!exists("metric_tensor_local") || !is.function(metric_tensor_local)) {
      stop("metric_tensor_local() must exist if use_metric = TRUE.")
    }
    ds <- numeric(length(dx))
    for (i in seq_along(dx)) {
      mid_x <- (x_path[i] + x_path[i + 1]) / 2
      mid_y <- (y_path[i] + y_path[i + 1]) / 2
      g <- metric_tensor_local(mid_x, mid_y)
      v <- c(dx[i], dy[i])
      ds[i] <- sqrt(t(v) %*% g %*% v)
    }
  } else {
    ds <- sqrt(dx^2 + dy^2)
  }
  
  path_length <- sum(ds)
  
  ## --- 6. Trapezoidal integration for area BELOW the geodesic ---
  z_below_mid <- (z_below_int[-1] + z_below_int[-length(z_below_int)]) / 2
  area_below  <- sum(z_below_mid * ds)
  
  ## ====== BRANCH: NO peak_df (old behavior: single list) ====== ##
  if (is.null(peak_df)) {
    
    ## --- 7a. Single cap_raw (from cap_value or global max) ---
    cap_raw <- NA_real_
    
    if (!is.null(cap_value)) {
      cap_raw <- cap_value
    } else {
      # Use quadratic surface at the *maximum* of the original z as before
      cap_raw <- max(z_surf_raw, na.rm = TRUE)
    }
    
    # Scale cap to match z
    if (standardize) {
      if (z_max == z_min) {
        cap_z <- 0
      } else {
        cap_z <- (cap_raw - z_min) / (z_max - z_min)
      }
    } else {
      cap_z <- cap_raw
    }
    
    ## --- 8a. Area ABOVE the geodesic up to the cap (single cap) ---
    delta <- cap_z - z
    delta_mid <- (delta[-1] + delta[-length(delta)]) / 2
    delta_mid <- pmax(delta_mid, 0)
    area_above <- sum(delta_mid * ds)
    
    return(list(
      area_below_geo_baseline   = area_below,
      area_between_geo_peak     = area_above,
      total_area_baseline_peak  = area_below + area_above,
      path_length               = path_length,
      cap_value_raw             = cap_raw,
      cap_value_scaled          = cap_z,
      standardized              = standardize,
      metric_based              = use_metric,
      signed_area               = signed_area
    ))
  }
  
  ## ====== BRANCH: peak_df PROVIDED (row-wise results) ====== ##
  
  # If both peak_df and cap_value are provided, ignore cap_value (or warn)
  if (!is.null(cap_value)) {
    warning("Both peak_df and cap_value were provided. cap_value is ignored; ",
            "using quadratic surface evaluated at each peak_df row instead.")
  }
  
  # Number of peak candidates
  n_peak <- nrow(peak_df)
  
  ## --- 7b. Compute cap_raw for EACH row of peak_df using internal P_local ---
  cap_raw_vec <- P_local(peak_df[[x_col_peak]], peak_df[[y_col_peak]])
  
  # Scale each cap_raw to match z
  if (standardize) {
    if (z_max == z_min) {
      cap_z_vec <- rep(0, length(cap_raw_vec))
    } else {
      cap_z_vec <- (cap_raw_vec - z_min) / (z_max - z_min)
    }
  } else {
    cap_z_vec <- cap_raw_vec
  }
  
  ## --- 8b. For each peak candidate, compute area_above and totals ---
  area_above_vec <- numeric(n_peak)
  total_area_vec <- numeric(n_peak)
  
  for (j in seq_len(n_peak)) {
    cap_z_j <- cap_z_vec[j]
    delta_j <- cap_z_j - z           # vertical gap at each path point
    delta_mid_j <- (delta_j[-1] + delta_j[-length(delta_j)]) / 2
    delta_mid_j <- pmax(delta_mid_j, 0)
    area_above_j <- sum(delta_mid_j * ds)
    area_above_vec[j] <- area_above_j
    total_area_vec[j] <- area_below + area_above_j
  }
  
  ## --- 9b. Return peak_df augmented with the new columns ---
  out_df <- peak_df
  out_df$area_below_geo_baseline   <- area_below          # same for all rows
  out_df$area_between_geo_peak     <- area_above_vec
  out_df$total_area_baseline_peak  <- total_area_vec
  out_df$path_length               <- path_length
  out_df$cap_value_raw             <- cap_raw_vec
  out_df$cap_value_scaled          <- cap_z_vec
  out_df$standardized              <- standardize
  out_df$metric_based              <- use_metric
  out_df$signed_area               <- signed_area
  
  return(out_df)
}



# --- Run ---
euclidean_path_2d <- function(p1, p2, n = 100) {
  # p1 and p2 must each have columns x and y
  stopifnot(all(c("PROT", "CARB") %in% names(p1)))
  stopifnot(all(c("PROT", "CARB") %in% names(p2)))
  
  t <- seq(0, 1, length.out = n)
  
  data.frame(
    x = p1$PROT + t * (p2$PROT - p1$PROT),
    y = p1$CARB + t * (p2$CARB - p1$CARB)
  )
}
