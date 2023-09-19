#' @title Approximate FHG relationships
#' @description Approximate FHG relationships using both OLS and NLS methods
#' @param Q a stream flow time series
#' @param P a corresponding time series of a second hydraulic variable
#' @param type relationship being tested
#' @return data.frame
#' @family FHG
#' @export

compute_fhg = function(Q, P, type = "relation") {
  Qlog <- Ylog <- NULL
  
  df = data.frame(Q, P) %>%
    arrange(Q) %>%
    mutate(Qlog = log(Q), Ylog = log(P)) %>%
    filter(is.finite(Qlog), is.finite(Ylog))
  
  fit = lm(df$Ylog ~ df$Qlog)
  
  coef = exp(fit$coefficients[1])
  exp  = fit$coefficients[2]
  p2 = df$P
  q = df$Q
  model <-  suppressWarnings({
    nls(
      p2 ~ alpha * q ^ x,
      start = list(alpha = coef, x = exp),
      trace = FALSE,
      control = nls.control(
        maxiter = 50,
        tol = 1e-09,
        warnOnly = TRUE
      )
    )
  })
  
  s = model %>% summary()
  
  coef2 = s$coefficients[1, 1]
  exp2  = s$coefficients[2, 1]
  
  OLS = coef  * df$Q ^ exp
  NLS = coef2 * df$Q ^ exp2
  
  res  = list()
  
  res = data.frame(
    type = type,
    exp = exp,
    coef = coef,
    nrmse = nrmse(OLS, df$P),
    pb = pbias(OLS, df$P),
    method = "ols",
    row.names = NULL
  ) %>%
    rbind(
      data.frame(
        type = type,
        exp  = exp2,
        coef = coef2,
        nrmse = nrmse(NLS, df$P),
        pb = pbias(NLS, df$P),
        method = "nls",
        row.names = NULL
      )
    )
  
  
  arrange(res, nrmse)
}
