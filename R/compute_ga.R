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
#' @param type metric to evaluate
#' @return data.frame
#' @family FHG
#' @export

calc_nsga = function(df, allowance = .05, r, type = "nrmse") {
  
  if (type == "nrmse") {
    
    fitness <- function(x) {
      ## order: k,m,a,b,c,f
      ## order: k,m,a,b,c,f
      V.tmp = x[1]*df$Q^x[2]
      v = nrmse(V.tmp, df$V) 
      
      T.tmp = x[3]*df$Q^x[4]
      t =nrmse(T.tmp, df$TW) 
      
      D.tmp = x[5]*df$Q^x[6]
      d = nrmse(D.tmp, df$Y) 
      
      return(c(v,t,d))
    }
    
  } else {
    
    fitness <- function(x) {
      
      ## order: k,m,a,b,c,f
      v = mean(abs((df$V - (x[1] * df$Q ^ x[2])) / df$V))
      t = mean(abs((df$TW - (x[3] * df$Q ^ x[4])) / df$TW))
      d = mean(abs((df$Y - (x[5] * df$Q ^ x[6])) / df$Y))
      return(c(v, t, d))
    }
    
  }
  
  # Define Constraints
  c4 = function(x) {
    c1 = (x[1] * x[3] * x[5])
    c2 = (x[2] + x[4] + x[6])
    
    return(c((1 + allowance) - c1,
             c1 - (1 - allowance),
             (1 + allowance) - c2,
             c2 - (1 - allowance)
    ))
  }
  
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
  
  
  run_it = function(fitness, c4, seed, cprob, mprob) {
    set.seed(seed)
    error <-  floor_error <- rank.y <-  NULL
    nsga = nsga2(
      fitness,
      idim = 6,
      odim = 3,
      lower.bounds  = c(0, 0, 0, 0, 0, 0),
      upper.bounds  = c(3.5, 1, 642, 1, 20, 1),
      generations   = 200,
      popsize       = 32,
      cprob = cprob,
      mprob = mprob,
      constraints   = c4,
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
  
  m = c(.1, .2, .3)
  c = c(.5, .75, 1)
  
  g = expand.grid(m, c)
  
  o2 = rbind(run_it(fitness, c4, 1, g$Var2[1], g$Var1[1]),
             run_it(fitness, c4, 2, g$Var2[2], g$Var1[2]),
             run_it(fitness, c4, 3, g$Var2[3], g$Var1[3]),
             run_it(fitness, c4, 4, g$Var2[4], g$Var1[4]),
             run_it(fitness, c4, 5, g$Var2[5], g$Var1[5]),
             run_it(fitness, c4, 6, g$Var2[6], g$Var1[6]),
             run_it(fitness, c4, 7, g$Var2[7], g$Var1[7]),
             run_it(fitness, c4, 8, cprob = g$Var2[8], mprob = g$Var1[8]),
             run_it(fitness, c4, 9, cprob = g$Var2[9], mprob = g$Var1[9]))
  
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
