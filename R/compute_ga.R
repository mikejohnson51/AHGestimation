#' @title Calculate GA FHG
#' @param Q streamflow time series
#' @param Y river depth time series
#' @param V river velocity time series
#' @param TW river top width time series
#' @param allowance allowable deviation from continuity
#' @param r rrr
#' @param type metric to evaluate
#' @return data.frame
#' @export

calc_ga = function(Q,
                   Y,
                   V,
                   TW,
                   allowance = .05,
                   r,
                   type = "nrmse") {
  
  if (type == "nrmse") {
    
    fitness <- function(x) {
      ## order: k,m,a,b,c,f
      v = rmse(x[1] * Q ^ x[2], V) / mean(V)
      t = rmse(x[3] * Q ^ x[4], TW) / mean(TW)
      d = rmse(x[5] * Q ^ x[6], Y) / mean(Y)
      return(c(v, t, d))
    }
    
  } else {
    
    fitness <- function(x) {
      
      ## order: k,m,a,b,c,f
      v = mean(abs((V - (x[1] * Q ^ x[2])) / V))
      t = mean(abs((TW - (x[3] * Q ^ x[4])) / TW))
      d = mean(abs((Y - (x[5] * Q ^ x[6])) / Y))
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
    vals %>% data.frame() %>%
      mutate(id = 1:nrow(vals), error = rowSums(vals)) %>%
      mutate(rank = rank(error))
  }
  
  
  set.seed(1234)
  
  run_it = function(seed, fitness, c4) {
    #set.seed(seed)
    nsga = nsga2(
      fitness,
      6,
      3,
      lower.bounds  = c(0, 0, 0, 0, 0, 0),
      upper.bounds  = c(3.5, 1, 642, 1, 20, 1),
      generations   = 250,
      popsize       = 100,
      constraints   = c4,
      cdim = 4
    )
    
    vals = nsga$value[nsga$pareto.optimal, ]
    vals = vals[!duplicated(vals), ]
    par  = nsga$par[nsga$pareto.optimal, ]
    par  =  par[!duplicated(par), ]
    
    
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
  
  o = rbind(run_it(10291991, fitness, c4),
            run_it(1991, fitness, c4) ,
            run_it(1234, fitness, c4))
  
  err = apply(o, 1, fitness) %>% apply(2, sum)
  
  min = o[which.min(err), ]
  
  vp = (min[1] * Q ^ min[2])
  tp = (min[3] * Q ^ min[4])
  dp = (min[5] * Q ^ min[6])
  
  v = rmse(vp, V) / mean(V)
  t = rmse(tp, TW) / mean(TW)
  d = rmse(dp, Y) / mean(Y)
  
  v2 = pbias(vp, V)
  t2 = pbias(tp, TW)
  d2 = pbias(dp, Y)
  
  v3 = mean(abs((V - vp) / V))
  t3 = mean(abs((TW - tp) / TW))
  d3 = mean(abs((Y - dp) / Y))
  
  v4 = min(abs((V - vp) / V))
  t4 = min(abs((TW - tp) / TW))
  d4 = min(abs((Y - dp) / Y))
  
  v5 = max(abs((V - vp) / V))
  t5 = max(abs((TW - tp) / TW))
  d5 = max(abs((Y - dp) / Y))
  
  
  r$Y[3, ] = data.frame(
    type = "Y",
    exp = min[6],
    coef = min[5],
    nrmse = d,
    pb = d2,
    mean_ape = d3,
    min_ape = d4,
    max_ape = d5,
    method = "GA",
    stringsAsFactors = FALSE
  )
  r$TW[3, ] = data.frame(
    type = "TW",
    exp = min[4] ,
    coef = min[3],
    nrmse = t,
    pb = t2,
    mean_ape = t3,
    min_ape = t4,
    max_ape = t5,
    method = "GA",
    stringsAsFactors = FALSE
  )
  r$V[3, ]  = data.frame(
    type = "V",
    exp = min[2] ,
    coef = min[1],
    nrmse = v,
    pb = v2,
    mean_ape = v3,
    min_ape = v4,
    max_ape = v5,
    method = "GA",
    stringsAsFactors = FALSE
  )

  return(r)
}
