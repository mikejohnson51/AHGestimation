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
#' @family AHG
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

#' @title Properly estimate AHG values
#' @param df hydraulic data.frame with columns named (Q, V, TW, Y). Q and at least one other are required.
#' @param allowance allowed deviation from continuity
#' @param gen Number of generations to breed. 
#' @param pop Size of population
#' @inheritParams mco::nsga2
#' @param times how many times (seeds) should nsga2 be run
#' @param scale should a scale factor be applied to data pre NSGA-2 fitting
#' @param full_fitting should all fits be returned?
#' @param verbose should messages be emitted?
#' @return list
#' @family AHG
#' @export

ahg_estimate = function(df,
                        allowance = .05,
                        gen = 192,
                        pop = 200,
                        cprob = .4,
                        mprob = .4, 
                        times = 1,
                        scale = 1.5,
                        full_fitting = FALSE,
                        verbose = FALSE) {
  
  type <- NULL
  
  if(!"Q" %in% names(df)){ stop("Q must be present in df") }
  
  if(!any(grepl("TW|V|Y", names(df)))){ stop("At least one of TW, V, or Y must be present in df") }

  df = as.data.frame(select(df, any_of(c("date", "TW", "V", "Y", "Q"))))
  df = df[df > 0,]
  df = df[is.finite(rowSums(select(df, -any_of('date')))), ]
  df = df[complete.cases(df), ]
  
  ahg_y  = if ("Y" %in% names(df)) {
    compute_ahg(df$Q, df$Y, "Y")
  }
  
  ahg_tw = if ("TW" %in% names(df)) {
    compute_ahg(df$Q, df$TW, "TW")
  }
  
  ahg_v  = if ("V" %in% names(df)) {
    compute_ahg(df$Q, df$V, "V")
  }
  
  r = list(ahg_y, ahg_tw, ahg_v)
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
        df = df,
        allowance = allowance,
        r = r,
        scale = scale, 
        gen = gen,
        pop = pop,
        cprob = cprob,
        mprob = mprob, 
        times = times
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
    m = mismash(v = c("ols", "nls"), V = df$V, TW = df$TW, Y = df$Y, Q = df$Q, r, allowance)
  } else {
    return(arrange(bind_rows(r), type, nrmse))
}
  
  if(full_fitting){
    return(m)
  } else {
    return(m$summary)
  }
  

}

