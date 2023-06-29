#' @title Calculate GA AHG
#' @param Q streamflow timeseries
#' @param Y streamflow timeseries
#' @param V streamflow timeseries
#' @param TW streamflow timeseries
#' @param allowance streamflow timeseries
#' @param r streamflow timeseries
#' @return data.frame
#' @export

calc_ga_2 = function(Q, Y, V, TW, allowance = .05, r){

f2 <- function(x) {
  ## order: k,m,a,b,c,f
  V.tmp = x[1]*Q^x[2]
  v = hydroGOF::nrmse(V.tmp, V, norm = "maxmin") 
  
  T.tmp = x[3]*Q^x[4]
  t = hydroGOF::nrmse(T.tmp, TW, norm = "maxmin") 
  
  D.tmp = x[5]*Q^x[6]
  d = hydroGOF::nrmse(D.tmp, Y, norm = "maxmin") 
  
  return(sum(v,t,d))
}

f3 <- function(x) {
  ## order: k,m,a,b,c,f
  V.tmp = x[1]*Q^x[2]
  v  = hydroGOF::nrmse(V.tmp, V, norm = "maxmin")
  v2 = hydroGOF::pbias(V.tmp, V)

  T.tmp = x[3]*Q^x[4]
  t  = hydroGOF::nrmse(T.tmp, TW, norm = "maxmin")
  t2 = hydroGOF::pbias(T.tmp, TW)

  D.tmp = x[5]*Q^x[6]
  d  = hydroGOF::nrmse(D.tmp, Y, norm = "maxmin")
  d2 = hydroGOF::pbias(D.tmp, Y)

  return(c(v,t,d))
}


c1 <- function(x){ 
  !between(x[1]*x[3]*x[5],(1 - allowance), (1 + allowance))
  # t1 = min(0, x[1]*x[3]*x[5] - (1 - allowance))
  # t2 = max(0, x[1]*x[3]*x[5] - (1 + allowance))
  # abs(t1+t2)
}

c2 <- function(x){ 
  !between(x[2]+x[4]+x[6],(1 - allowance), (1 + allowance))
  # t1 = min(0, x[2]+x[4]+x[6] - (1 - allowance))
  # t2 = max(0, x[2]+x[4]+x[6] - (1 + allowance))
  # abs(t1+t2)
}

fitness <- function(x){
  f <- -f2(x)
  
  pen <- 1e6
  
  if(any(c1(x),c2(x))){ 
    f - pen 
  } else {
    f
  }
  # penalty1 <- c1(x)*pen
  # penalty2 <- c2(x)*pen
  # 
  # f - penalty1 - penalty2
}

GA <- GA::ga(type = "real-valued",
             fitness =  fitness,
             lower = c(0,0, 0,0, 0,0),
             upper = c(3.5, 1, 642, 1, 20, 1), 
             # suggestions = c(r$V$coef[1],  r$V$exp[1],
             #                 r$TW$coef[1], r$TW$exp[1],
             #                 r$Y$coef[1],  r$Y$exp[1]),
             names = c('k','m','a','b','c','f'),
             seed = 10291991, 
             maxiter = 2500, run = 250,
             monitor = T, optim = TRUE)

sol = GA@solution[1,]

out = f3(GA@solution[1,])
c1(sol)

r$Y[3,] = data.frame(exp = sol[6] , 
                     coef = sol[5], 
                     nrmse = out[3], 
                     pb = out[6], 
                     method = "GA", 
                     stringsAsFactors = FALSE)

r$TW[3,] = data.frame(exp = sol[4] , 
                      coef = sol[3], 
                      nrmse = out[2],
                      pb = out[5], method = "GA", 
                      stringsAsFactors = FALSE)

r$V[3,]  = data.frame(exp = sol[2], 
                      coef = sol[1], 
                      nrmse = out[1] , 
                      pb = out[4], 
                      method = "GA", 
                      stringsAsFactors = FALSE)

return(r)

}
