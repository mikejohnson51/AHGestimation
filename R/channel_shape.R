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
#' @return depth values every 1m along the cross section  
#' @family hydraulics 
#' @export

cross_section <- function(r, TW = 30, Ymax = 2){
  
  TW <- as.integer(TW)
  half <- TW/2
  x_list <- numeric()
  z_list <- numeric()

  for(i in 1:half){
    Z <- Ymax * (2/TW)**r * (i**r)
    x_list <- append(x_list, i)
    z_list <- append(z_list, Z)
  }
  
  df <- data.frame(
    ind = seq_along(z_list),
    x = seq(0, TW, length.out = 2*length(z_list)),
    Y = c(rev(z_list), z_list),
    A = NA
  )
  
  for(i in seq(nrow(df))){
    df$A[i] <- .findCA(df, depth = df$Y[i])
  }
  
  df

}

.findCA <- function(df, depth){
  
  Y <-  NULL
  t <- filter(df, Y < depth)
  
  x <- pmax(0, AUC(x = t$x, 
      y = rep(depth, nrow(t)), 
      absolutearea = FALSE) - 
  AUC(x = t$x, 
      y = t$Y, 
      absolutearea = FALSE))

  ifelse(is.na(x), 0, x)
}
