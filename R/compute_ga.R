min_max = function(df, s = 2) {
  
  method <- NULL
  
  if (is.null(scale)) {
    list(
      value = c('k', 'm', 'a', 'b', 'c', 'f'),
      low   = c(0, 0, 0, 0, 0, 0),
      high  = c(3.5, 1, 642, 1, 20, 1)
    )
  } else {
    o = filter(compute_fhg(df$Q, df$Y, "Y"), method == "nls")
      Y_e_l =  (1 / s) * o$exp
      Y_e_h =  s     * o$exp
      Y_c_l = (1 / s) * o$coef
      Y_c_h = s     * o$coef
      Y_e = o$nrmse
      Y_e_nls = o$exp
      Y_c_nls = o$coef
    
    o = filter(compute_fhg(df$Q, df$V, "V"), method == "nls")
      V_e_l =  (1 / s) * o$exp
      V_e_h =  s     * o$exp
      V_c_l = (1 / s) * o$coef
      V_c_h = s     * o$coef
      V_e = o$nrmse 
      V_e_nls = o$exp
      V_c_nls = o$coef
    
    o = filter(compute_fhg(df$Q, df$TW, "TW"), method == "nls")
      TW_e_l =  (1 / s) * o$exp
      TW_e_h =  s     * o$exp
      TW_c_l = (1 / s) * o$coef
      TW_c_h = s     * o$coef
      TW_e     = o$nrmse
      TW_e_nls = o$exp
      TW_c_nls = o$coef
    
    list(
      value = c('k', 'm', 'a', 'b', 'c', 'f'),
      error = c(Y_e, NA, V_e, NA, TW_e),
      base  = c(V_c_nls, V_e_nls, TW_c_nls, TW_e_nls, Y_c_nls, Y_e_nls),
      low   = c(V_c_l, V_e_l, TW_c_l, TW_e_l, Y_c_l, Y_e_l),
      high  = c(V_c_h, V_e_h, TW_c_h, TW_e_h, Y_c_h, Y_e_h)
    )
  }
}

run_it = function(fitness, 
                  constraints, 
                  seed, 
                  cprob, 
                  mprob,  
                  gen, 
                  pop,
                  lower.bounds  = c(0, 0, 0, 0, 0, 0),
                  upper.bounds  = c(3.5, 1, 642, 1, 20, 1)) {
  
  set.seed(seed)
  
  cont = function(par) {
    c1 <- c2 <- c1_e <- c2_e <- c3 <- NULL
    data.frame(
      id = 1:nrow(par),
      r = par[, 6] / par[, 4],
      c1 = par[, 1] * par[, 3] * par[, 5],
      c2 = par[, 2] + par[, 4] + par[, 6]
    ) %>%
      mutate(
        c1_e = abs(1 - c1),
        c2_e = abs(1 - c2),
        c3 = c1_e + c2_e,
        rank = rank(c3)
      )
  }
  
  errors = function(vals) {
    error <- NULL
    vals %>% data.frame() %>%
      mutate(id = 1:nrow(vals), error = rowSums(vals)) %>%
      mutate(rank = rank(error))
  }
  
  error <-  floor_error <- rank.y <-  NULL
  
  nsga = nsga2(
    fitness,
    idim = 6,
    odim = 3,
    lower.bounds  = lower.bounds,
    upper.bounds  = upper.bounds,
    generations   = gen,
    popsize       = pop,
    cprob = cprob,
    mprob = mprob,
    constraints   = constraints,
    cdim = 4
  )
  
  vals = nsga$value[nsga$pareto.optimal, ]
  
  if(is.null(nrow(vals))){
    vals = matrix(vals, byrow = 1, nrow = 1)
  }
  vals = vals[!duplicated(vals), ]
  par  = nsga$par[nsga$pareto.optimal, ]
  
  if(is.null(nrow(par))){
    par = matrix(par, byrow = 1, nrow = 1)
  }
  par  = par[!duplicated(par), ]
  
  
  
  if (!is.null(dim(vals))) {
    err  = left_join(errors(vals), cont(par), by = "id") %>%
      mutate(floor_error = floor(error)) %>%
      filter(floor_error == min(floor_error)) %>%
      slice_min(rank.y, n = 1)
    
    min = par[err[1, "id"],]
    
  } else {
    min = vals
  }
  min
}



#' Percent Bias
#' @description Percent Bias between sim and obs, with treatment of missing values.
#' @param sim numeric vector simulated values
#' @param obs numeric vector observed values
#' @return numeric
#' @family evaluation
#' @export

pbias = function(sim, obs){
  vi <- valindex(sim, obs)
  
  if (length(vi) > 0) {	 
    obs <- obs[vi]
    sim <- sim[vi]
    
    # length of the data sets that will be considered for the computations
    n <- length(obs)
    
    denominator <- sum( obs )
    
    if (denominator != 0) {      
      pbias <- 100 * ( sum( sim - obs ) / denominator )
      pbias <- round(pbias, 2)     
    } else {
      pbias <- NA
      warning("'sum((obs)=0' -> it is not possible to compute 'pbias' !")  
    } 
  } else {
    pbias <- NA
    warning("There are no pairs of 'sim' and 'obs' without missing values !")
  } 
  
  return( pbias )
}

#' Normalized Root Mean Square Error
#' @description Normalized root mean square error (NRMSE) between sim and obs, with treatment of missing values
#' @param sim numeric vector simulated values
#' @param obs numeric vector observed values
#' @return numeric
#' @family evaluation
#' @export
#' 
nrmse <- function(sim, obs) {
  
  # index of those elements that are present both in 'sim' and 'obs' (NON- NA values)
  vi <- valindex(sim, obs)
  
  if (length(vi) > 0) {	 
    obs <- obs[vi]
    sim <- sim[vi]
    
    cte <- ( max(obs, na.rm= TRUE) - min(obs, na.rm =TRUE) )
    
    rmse <- rmse(sim, obs, TRUE) 
    
    if (max(obs, na.rm= TRUE) - min(obs, na.rm= TRUE) != 0) {     
      nrmse <- rmse / cte     
    } else {
      nrmse <- NA
      warning("'obs' is constant -> it is not possible to compute 'nrmse' !")  
    } 
  } else {
    nrmse <- NA
    warning("There are no pairs of 'sim' and 'obs' without missing values !")
  } # ELSE end
  
  return( round( 100*nrmse, 2) )
  
}

valindex <- function(sim, obs) {  
    index <- which(!is.na(sim) & !is.na(obs))
    if (length(index)==0) warning("'sim' and 'obs' are empty or they do not have any common pair of elements with data !!")
    return( index  )
} 

rmse <- function (sim, obs, na.rm=TRUE) {
  
  if ( length(obs) != length(sim) ) 
    stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")     
  
  rmse <- sqrt( mean( (sim - obs)^2, na.rm = na.rm) )
  
  return(rmse)
  
}

#' @title Calculate GA FHG
#' @param df hydraulic data.frame
#' @param allowance allowable deviation from continuity
#' @param r fit list
#' @inheritParams fhg_estimate
#' @return data.frame
#' @family FHG
#' @export

calc_nsga = function(df, 
                     allowance = .05, 
                     r, 
                     scale = 2, 
                     gen = 96,
                     pop = 500,
                     cprob = .8,
                     mprob = .05, 
                     times = 1) {
  
    fitness <- function(x) {
      ## order: k,m,a,b,c,f
      V.tmp = x[1]*df$Q^x[2]
      T.tmp = x[3]*df$Q^x[4]
      D.tmp = x[5]*df$Q^x[6]
  
      return(c(nrmse(V.tmp, df$V),
               nrmse(T.tmp, df$TW),
               nrmse(D.tmp, df$Y)))
    }

  # Define Constraints
  constraints = function(x) {
    c1 = (x[1] * x[3] * x[5])
    c2 = (x[2] + x[4] + x[6])
    
    return(c((1 + allowance) - c1,
             c1 - (1 - allowance),
             (1 + allowance) - c2,
             c2 - (1 - allowance)
    ))
  }
  
  mm = min_max(df, scale)
  
  o2 = list()
  
  for(i in seq_along(times)){
    o2[[i]] = run_it(
           fitness, 
           constraints, 
           i,
           cprob, 
           mprob,  
           gen, 
           pop,
           lower.bounds = mm$low, 
           upper.bounds = mm$high)
  }
  
  o2 = do.call('rbind', o2)
  #   # https://ai.stackexchange.com/questions/12019/how-to-find-optimal-mutation-probability-and-crossover-probability
  #   o2 = rbind(# high crossover, low mutation
  #              run_it(fitness, c4, 1, .8, .05, lower.bounds = mm$low, upper.bounds = mm$high),
  #              run_it(fitness, c4, 4, .8, .05, lower.bounds = mm$low, upper.bounds = mm$high),
  #              run_it(fitness, c4, 7, .8, .05, lower.bounds = mm$low, upper.bounds = mm$high),
  #              
  #              # high-ish crossover, low-ish mutation
  #              run_it(fitness, c4, 2, .6, .2, lower.bounds = mm$low, upper.bounds = mm$high),
  #              run_it(fitness, c4, 5, .6, .2, lower.bounds = mm$low, upper.bounds = mm$high),
  #              run_it(fitness, c4, 8, .6, .2, lower.bounds = mm$low, upper.bounds = mm$high),
  #              
  #              # moderate crossover, moderate mutation
  #              run_it(fitness, c4, 3, .4, .4, lower.bounds = mm$low, upper.bounds = mm$high),
  #              run_it(fitness, c4, 6, .4, .4, lower.bounds = mm$low, upper.bounds = mm$high),
  #              run_it(fitness, c4, 9, .4, .4, lower.bounds = mm$low, upper.bounds = mm$high)
  #   )
  # })
  
  err = apply(o2, 1, fitness) %>% 
    apply(2, sum)

  min = o2[which.min(err), ]

  vp = (min[1] * df$Q ^ min[2])
  tp = (min[3] * df$Q ^ min[4])
  dp = (min[5] * df$Q ^ min[6])
  
  nrmse_v = nrmse(vp, df$V) 
  nrmse_t = nrmse(tp, df$TW) 
  nrmse_d = nrmse(dp, df$Y) 
  
  pbias_v = pbias(vp,df$V)
  pbias_t = pbias(tp, df$TW)
  pbias_d = pbias(dp, df$Y)
  
  r$Y[3, ] = data.frame(
    type = "Y",
    exp = min[6],
    coef = min[5],
    nrmse = nrmse_d,
    pb = pbias_d,
    method = "nsga2")
  
  r$TW[3, ] = data.frame(
    type = "TW",
    exp = min[4] ,
    coef = min[3],
    nrmse = nrmse_t,
    pb = pbias_t,
    method = "nsga2",
    stringsAsFactors = FALSE
  )
  
  r$V[3, ]  = data.frame(
    type = "V",
    exp = min[2] ,
    coef = min[1],
    nrmse = nrmse_v,
    pb = pbias_v,
    method = "nsga2",
    stringsAsFactors = FALSE
  )

  return(r)
}
