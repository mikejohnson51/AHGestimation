check_validity = function(r, type = "nls", allowance = .05){
  
  method <- NULL
  is_between = function(x, allowance) {
    between(x, 1 - allowance, 1 + allowance)
  }
  
  tw = filter(r$TW, method == !!type)
  y = filter(r$Y, method == !!type)
  v = filter(r$V, method == !!type)
  
  c1 = is_between(prod(y$coef, tw$coef, v$coef), allowance)
  c2 = is_between(sum(y$exp, tw$exp,v$exp), allowance)
  
  sum(c1, c2) == 2
  
}


#' @title Report best optimal
#' @param best best performing method (character string)
#' @param check values to check against
#' @param verbose should messages be emitted
#' @return vector
#' @family FHG
#' @export

best_optimal = function(best, check, verbose = TRUE) {
  if (verbose) {
    message(
      paste0(
        "The best performing method (",
        best,
        ") ",
        ifelse(check, "is", "is not"),
        " physically valid ... ",
        "\n"
      )
    )
  }
  
  return(check)
}

#' @title Properly estimate FHG values
#' @param df hydraulic data.frame
#' @param allowance allowed deviation from continuity
#' @param verbose should messages be emitted?
#' @return list
#' @details
#' If filter_data is TRUE. An best fit curve is fit to the raw data relationships and is used to predict
#' A value for each Q. Data is only retained if the provided value is within +/- the `filter_threshold` (as a percent) of 
#' the estimated value
#' @family FHG
#' @export

fhg_estimate = function(df,
                        allowance = .05,
                        verbose = FALSE) {
  
  type <- NULL
  
  if(!"Q" %in% names(df)){ stop("Q must be present in df") }
  
  if(!any(grepl("TW|V|Y", names(df)))){ stop("At least one of TW, V, or Y must be present in df") }

  df = select(df, any_of(c("date", "TW", "V", "Y", "Q")))
  df = df[df > 0,]
  df = df[is.finite(rowSums(select(df, -any_of('date')))), ]
  df = df[complete.cases(df), ]
  
  fhg_y  = if ("Y" %in% names(df)) {
    compute_fhg(df$Q, df$Y, "Y")
  }
  
  fhg_tw = if ("TW" %in% names(df)) {
    compute_fhg(df$Q, df$TW, "TW")
  }
  
  fhg_v  = if ("V" %in% names(df)) {
    compute_fhg(df$Q, df$V, "V")
  }
  
  r = list(fhg_y, fhg_tw, fhg_v)
  r = Filter(Negate(is.null), r)
  n = vector()
  for(i in 1:length(r)){ n = append(n, r[[i]]$type[1])}
  names(r)  = n
  
  best = NULL
  
  for (i in 1:length(r)) {
    best[i] = r[[i]]$method[which.min(r[[i]]$nrmse)]
  }
  
  if (verbose) {
    message(paste0(
      unique(toupper(best)),
      " performs best for the Q-",
      names(r),
      ' realtionship\n'
    ))
  }
  
  if (length(r) == 3) {
    
    nls_viable = check_validity(r, "nls", allowance)
    ols_viable = check_validity(r, "ols", allowance)
    
    if (verbose) {
      message(paste(
        "NLS",
        ifelse(nls_viable, "meets", "does not meet"),
        "continuity ... "
      ))
      message(paste(
        "OLS",
        ifelse(ols_viable, "meets", "does not meet"),
        "continuity ... ",
        "\n"
      ))
    }
    
    best = unique(best)
    
    cond = best_optimal(best,
                        ifelse(best[1] == "nls", nls_viable, ols_viable),
                        verbose = verbose)
    
    if (any(!cond)) {
      if (verbose) {
        message(
          "Launching Evolutionary Algorithm... allowance (",
          allowance,
          ") "
        )
      }
      
      r = calc_nsga(
        df,
        allowance,
        r
      )
      
      m = mismash(v = c("ols", "nls", "nsga2"),
                  V = df$V,
                  TW = df$TW,
                  Y = df$Y,
                  Q = df$Q,
                  r,
                  allowance)
      
    } else{
      m = mismash(v = c("ols", "nls"), df$V, df$TW, df$Y, df$Q, r, allowance)
    }
    
  } else if (length(r) == 2) {
    m = mismash(v = c("ols", "nls"), df$V, df$TW, df$Y, df$Q, r, allowance)
  } else {
    return(arrange(bind_rows(r), type, nrmse))
}
  
  return(m$summary)

}

