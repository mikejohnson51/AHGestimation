#' @title Approximate channel coeffiecent
#' @description Approximate the hydraulic values from AHG fit
#' @param fit output of ahg_estimate
#' @return numeric 
#' @family hydraulics 
#' @export

compute_hydraulic_params <- function(fit){
  
  n <- names(fit)
  fit$r <- fit$Y_exp / fit$TW_exp
  fit$p <- fit$V_exp / fit$Y_exp
  
  fit$d <- (1 + fit$r + (fit$r * fit$p))
  
  fit$R <- (1 + fit$r) / fit$r
  
  fit$bd <- 1 / fit$d
  
  fit$fd <- fit$r / fit$d
  
  fit$md <- (fit$p * fit$r) / fit$d
  
  select(fit, -any_of(n))
}

#' @title Approximate Roughness
#' @description Approximate median roughness using Manning Equation
#' @param df a data.frame with at least Y and V. 
#' @param S reach scale longitudinal slope (m/m). Default mean of the nhdplusV2
#' @return numeric
#' @family hydraulics 
#' @export

compute_n <- function(df, S = .02){
  # Method from: https://zenodo.org/record/7868764
  if(all(c("Y", "V") %in% names(df))){
    # median streamdepth as approximation of R
    r <- (median(df$Y)) ^ (2/3)
    # Median V for V
    v <- median(df$V)
    # S from NHDPlus or other
    s <- sqrt(S)
    
    (r * s) / v
  } else { 
    stop("Can't compute n without V and Y")
  }
}

#' @title Approximate channel shape
#' @description Get a list of points from x axis of a cross section 
#' and max depth and produce depth values for 
#' those points based on channel shape
#' @param r The corresponding Dingman's r coefficient 
#' @param TW width of the channel at bankfull 
#' @param Ymax maximum depth of the channel at bankfull 
#' @param n the number of points to construct in the XS
#' @return depth values every 1m along the cross section  
#' @family hydraulics 
#' @export
#' 
cross_section <- function(r, TW = 30, Ymax = 2, n = 30, interpolate=FALSE){
  
  TW <- as.integer(TW)
  half <- ceiling(n/2)
  
  # Compute half of a symmetrical channel
  x_list <- 1:half

  z_list <- Ymax * (2/TW)^r * (x_list^r)
  
  # Handle Rectangular shaped channels that results in NaN or Inf values in z_list 
  z_list[is.nan(z_list) | is.infinite(z_list)] <- 0
  # Add banks relative elevation
  z_list[length(z_list)] <- Ymax
  
  if (interpolate) {
    # Perform interpolation
    interpolated <- interpolate_values(x_list, z_list)
    # Build a dataframe
    df <- data.frame(
      ind = seq_along(interpolated$x),
      x = interpolated$x, # add in the last point
      Y = interpolated$y, # add in the last point
      A = NA
    )
  } else {
    # Build a dataframe
    df <- data.frame(
      ind = seq_along(z_list),
      x = seq(0, TW, length.out = 2*length(z_list)),
      Y = c(rev(z_list), z_list),
      A = NA
    )
    
  }
  
  # Calculate area per stage
  df$A <- sapply(df$Y, function(y) .findCA(df, depth = y))
  return(df)
}

.findCA <- function(df, depth) {
  # Filter df to keep only rows below given stage
  t <- df[df$Y < depth, ]
  
  # Ensure t is not empty before calculating AUC
  if (nrow(t) == 0) {
    return(0)
  }
  
  # Calculate AUC for the given depth and the actual Y values
  auc_depth <- AUC(x = t$x, y = rep(depth, nrow(t)), absolutearea = FALSE)
  auc_y <- AUC(x = t$x, y = t$Y, absolutearea = FALSE)
  
  # Ensure non-negative values
  x <- pmax(0, auc_depth - auc_y)
  
  # Return x, replacing NA with 0
  return(ifelse(is.na(x), 0, x))
}

#' @title interpolate between points
#' @description given a array of coordinates of points populate in between them
#' @param x The x coordinates 
#' @param y The y coordinates 
#' @param spacing euclidean distance between points 
#' @return a list of new x and y coordinates   
#' @family hydraulics 
#' @export
#' 
interpolate_values <- function(x, y, spacing = 0.1) {
  # Calculate the differences between consecutive points
  x <- c(x, x[[length(x)]]+x)
  y <- c(rev(y),y)
  dx <- diff(x)
  dy <- diff(y)
  
  # Calculate the Euclidean distances between consecutive points
  distances <- sqrt(dx^2 + dy^2)
  
  # Calculate the cumulative distance along the points
  cumulative_distances <- c(0, cumsum(distances))
  
  # Create a new distance array at every `spacing` units
  new_distances <- seq(0, max(cumulative_distances), by = spacing)
  
  # Interpolate x and y coordinates based on the cumulative distances
  interp_x <- approx(cumulative_distances, x, xout = new_distances, method = "linear")$y
  interp_y <- approx(cumulative_distances, y, xout = new_distances, method = "linear")$y
  
  return(list(x = c(interp_x, x[length(x)]), y = c(interp_y, y[length(y)]))) # Add in the right bank
}
