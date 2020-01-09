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

mismash = function(v, params, V, TW, Y, Q, r, allowance){
  
  fit =  function(x,V, TW, Y, Q) {
    ## order: k,m,a,b,c,f
    v = (x[1]*Q^x[2]) %>% hydroGOF::nrmse(V,  norm = "maxmin")
    t = (x[3]*Q^x[4]) %>% hydroGOF::nrmse(TW, norm = "maxmin")
    d = (x[5]*Q^x[6]) %>% hydroGOF::nrmse(Y,  norm = "maxmin")
    
    return(c(v,t,d))
  }
  
  l = list()
  
  for(i in 1:params){
    l[[i]] = v
  }
  
  g = expand.grid(l, stringsAsFactors = F) %>% mutate(c1 = NA, c2 = NA, viable = NA,
                                                      Verror = NA,
                                                      TWerror = NA, Yerror = NA, tot_error = NA)
  
  for(i in 1:nrow(g)){
  
    Ytmp = c(r$Y$coef[which(r$Y$method ==   g[i,1])],
             r$Y$exp [which(r$Y$method ==   g[i,1])])
    Ttmp = c(r$TW$coef[which(r$TW$method == g[i,2])],
             r$TW$exp [which(r$TW$method == g[i,2])])
    Vtmp = c(r$V$coef[which(r$V$method ==   g[i,3])],
             r$V$exp [which(r$V$method ==   g[i,3])])
    
    ttt = c(Vtmp, Ttmp, Ytmp)
    
    g$c1[i] = prod(ttt[c(1,3,5)])
    g$c2[i] = sum(ttt[c(2,4,6)])
    
    g$viable[i] = sum(
      between(g$c1[i], 1-allowance, 1+allowance),
      between(g$c2[i], 1-allowance, 1+allowance)) == 2
    f = fit(ttt,V, TW, Y, Q)
    g$tot_error[i] = sum(f)
    g$Verror[i] = f[1]
    g$Yerror[i] = f[3]
    g$TWerror[i] = f[2]
    
  }

  g = g %>% arrange(tot_error)
  best = g %>% slice(n = 1) %>% mutate(condition = "best")
  physical = g %>% filter(viable == TRUE) %>% slice(n = 1) %>% mutate(condition = "physical")
  ols = g %>% filter(Var1 == "ols", Var2 == "ols", Var3 == "ols") %>% mutate(condition = "ols")
  nls = g %>% filter(Var1 == "nls", Var2 == "nls", Var3 == "nls") %>% mutate(condition = "nls")
  
  return(list(g = g %>% arrange(!viable, tot_error), summary = bind_rows(best, physical, ols, nls) %>%  arrange(!viable, tot_error)))
}