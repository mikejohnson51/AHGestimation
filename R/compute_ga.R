#' @title Calculate GA AHG
#' @param Q streamflow timeseries
#' @param Y streamflow timeseries
#' @param V streamflow timeseries
#' @param TW streamflow timeseries
#' @param allowance streamflow timeseries
#' @param r streamflow timeseries
#' @return data.frame
#' @export

calc_ga = function(Q, Y, V, TW, allowance = .05, r){
  
  set.seed(10291991)
  Q = Q
  V = V
  Y =Y
  TW = TW
  allowance = allowance

fitness <- function(x) {
  ## order: k,m,a,b,c,f
  v = (x[1]*Q^x[2]) %>% hydroGOF::nrmse(V,  norm = "maxmin")
  t = (x[3]*Q^x[4]) %>% hydroGOF::nrmse(TW, norm = "maxmin")
  d = (x[5]*Q^x[6]) %>% hydroGOF::nrmse(Y,  norm = "maxmin")
  
  return(c(v,t,d))
}

# Define Constraints
c4 = function(x){
  c1 = (x[1] * x[3] * x[5])
  c2 = (x[2] + x[4] + x[6])
  
  return(c((1+allowance) - c1, c1 - (1-allowance),
           (1+allowance) - c2, c2 - (1-allowance)))
}


out = mco::nsga2(fitness, 6, 3,
                 lower.bounds  = c(0,0, 0,0, 0,0),
                 upper.bounds =  c(3.5, 1, 642, 1, 20, 1),
                 generations= 150,
                 popsize = 160,
                 constraints = c4, cdim = 4)

vals = out$value[out$pareto.optimal,]
par = out$par[out$pareto.optimal,]
min_index = which.min(rowSums(vals)) 
min = par[min_index,]

v = (min[1]*Q^min[2]) %>% hydroGOF::nrmse(V,  norm = "maxmin")
t = (min[3]*Q^min[4]) %>% hydroGOF::nrmse(TW, norm = "maxmin")
d = (min[5]*Q^min[6]) %>% hydroGOF::nrmse(Y,  norm = "maxmin")

v2 = (min[1]*Q^min[2]) %>% hydroGOF::pbias(V)
t2 = (min[3]*Q^min[4]) %>% hydroGOF::pbias(TW)
d2 = (min[5]*Q^min[6]) %>% hydroGOF::pbias(Y)

r$Y[3,] = data.frame(exp = min[6] , coef = min[5], nrmse = d, pb = d2, method = "GA", 
                     stringsAsFactors = FALSE)
r$TW[3,] = data.frame(exp = min[4] , coef = min[3], nrmse = t, pb = t2, method = "GA", 
                     stringsAsFactors = FALSE)
r$V[3,]  = data.frame(exp = min[2] , coef = min[1], nrmse = v , pb = v2, method = "GA", 
                     stringsAsFactors = FALSE)

return(r)
}
