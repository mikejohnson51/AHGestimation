#' @title Compute all combos!
#' @param v 
#' @param params 
#' @param V 
#' @param TW 
#' @param Y 
#' @param Q 
#' @param r
#' @param allowance
#' @return list
#' @export

mismash = function(v, V, TW, Y, Q, r, allowance){
  
  fit =  function(x,V, TW, Y, Q) {
    ## order: k,m,a,b,c,f
    v = (x[1]*Q^x[2]) %>% hydroGOF::rmse(V) /mean(V)
    t = (x[3]*Q^x[4]) %>% hydroGOF::rmse(TW) /mean(TW)
    d = (x[5]*Q^x[6]) %>% hydroGOF::rmse(Y)/mean(Y)

    return(c(v,t,d))
  }


  l = list()
  
  for(i in 1:length(v)){
    l[[i]] = v
  }
  
  g = expand.grid(l, stringsAsFactors = F) %>% 
    mutate(c1 = NA, c2 = NA, viable = NA,
           Y_error = NA, TW_error = NA, V_error = NA, tot_error = NA)
  
  for(i in 1:nrow(g)){
    Ytmp = c(filter(r$Y, method == g[i,1])$coef, filter(r$Y, method == g[i,1])$exp)
    Ttmp = c(filter(r$TW, method == g[i,2])$coef, filter(r$TW, method == g[i,2])$exp)
    Vtmp = c(filter(r$V, method == g[i,3])$coef, filter(r$V, method == g[i,3])$exp)
    
    ttt = c(Vtmp, Ttmp, Ytmp)
    
    g$c1[i] = round(prod(ttt[c(1,3,5)]), 3)
    g$c2[i] = round(sum(ttt[c(2,4,6)]), 3)
    
    g$viable[i] = sum(
      between(g$c1[i], 1-allowance, 1+allowance),
      between(g$c2[i], 1-allowance, 1+allowance)) == 2
    
    f = fit(x = ttt,V, TW, Y, Q)
    
    g$tot_error[i] = sum(f)
    g$V_error[i] = f[1]
    g$TW_error[i] = f[2]
    g$Y_error[i] = f[3]
    
  }

  g = g %>% 
    arrange(tot_error) %>% 
    rename(Y = Var1, TW = Var2, V = Var3)
 
  physical = g %>% 
    filter(viable == TRUE) %>% 
    mutate(improve = tot_error[nrow(.)] - tot_error[1])  %>% 
    slice(1) %>% 
    mutate(condition = "bestValid")
  
  combo =  g %>% 
    filter(apply(g,1,function(r){length(unique(r[1:3]))}) != 1) %>% 
    filter(viable) %>% 
    slice(1) %>%
    mutate(condition = "combo")
  
  ols = g %>% 
    filter(Y == "ols", TW == "ols", V == "ols") %>% 
    mutate(condition = "ols")
  
  nls = g %>% 
    filter(Y == "nls", TW == "nls", V == "nls") %>% 
    mutate(condition = "nls")
  
  if("GA" %in% v){
  ga = g %>% 
    filter(Y == "GA", TW == "GA", V == "GA") %>% 
    mutate(condition = "GA")
  } else {
    ga = NULL
  }
  
  return(list(g = g %>% 
                arrange(!viable, tot_error), 
                summary = bind_rows(combo, ols, nls, ga) %>%  
                dplyr::select(-c1, -c2) %>% 
                arrange(!viable, tot_error)))
}
