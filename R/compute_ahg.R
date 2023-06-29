#' @title Approximate AHG realtionships
#' @description Approxiamte AHG realtionships using both OLS and NLS methods
#' @param Q a streamflow timeseries
#' @param P a corresponding timeseries of a second hydraulic variable
#' @return data.frame
#' @export
#' @importFrom dplyr arrange mutate filter
#' @importFrom stats nls nls.control
#' @importFrom hydroGOF nrmse pbias

compute_ahg = function(Q,P, type = NA){
  
  df = data.frame(Q,P) %>% 
    arrange(Q) %>% 
    mutate(Qlog = log(Q), Ylog = log(P)) %>% 
    filter(is.finite(Qlog))
  
  fit = lm(df$Ylog ~ df$Qlog)
  
  coef = exp(fit$coefficients[1])
  exp  = fit$coefficients[2]
  
  model <-  suppressWarnings({
    nls(P ~ alpha * Q ^ x,  start = list(alpha = coef, x = exp), trace = F,
        control = nls.control(maxiter = 50,
                              tol=1e-09, 
                              warnOnly=TRUE))
  })
  
  s = model %>% summary()
  coef2 = s$coefficients[1,1]
  exp2  = s$coefficients[2,1]
  
  OLS = coef  * df$Q ^ exp
  NLS = coef2 * df$Q ^ exp2
  
  res  = list()
  
  res[[type]] = data.frame(
    type = type,
    exp = exp,
    coef = coef,
    nrmse = hydroGOF::rmse(OLS, df$P) / mean(df$P),
    pb = hydroGOF::pbias(OLS, df$P),
    mean_ape = mean(abs((df$P - OLS) / df$P)),
    min_ape = min(abs((df$P - OLS) / df$P)),
    max_ape = max(abs((df$P - OLS) / df$P)),
    method = "ols",
    row.names = NULL,
    stringsAsFactors = FALSE) %>% 
  rbind(
    data.frame(
    type = type,
    exp  = exp2,
    coef = coef2,
    nrmse  = hydroGOF::rmse(NLS, df$P) / mean(df$P),
    pb = hydroGOF::pbias(NLS, df$P),
    mean_ape = mean(abs((df$P - NLS) / df$P)),
    min_ape = min(abs((df$P - NLS) / df$P)),
    max_ape = max(abs((df$P - NLS) / df$P)),
    method = "nls",
    row.names = NULL,
    stringsAsFactors = FALSE))
      

  res
}
