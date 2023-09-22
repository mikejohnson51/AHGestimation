#' @title Compute all permutations of fitting types
#' @param v values
#' @param V Velocity time series
#' @param TW Top width time series
#' @param Y Depth time series
#' @param Q Discharge time series
#' @param r list of fitting methods and parameters
#' @param allowance Allowable deviation from continuity
#' @return list
#' @family FHG
#' @export

mismash = function(v, V, TW, Y, Q, r, allowance){
  
  viable <- tot_error <-  Y_method <-  TW_method <- V_method <-  c1 <-  c2 <- NULL
  
  fit =  function(g, ind, V, TW, Y, Q) {

    x = g[ind,]
    
    if(!is.null(V)){
    g$V_error[ind] = tryCatch({
      (x$V_coef*Q^x$V_exp) %>% rmse(V) /mean(V)},
      error = function(e){
        NA
      })
    }
      
    if(!is.null(TW)){
    g$TW_error[ind] = tryCatch({
      (x$TW_coef*Q^x$TW_exp) %>% rmse(TW) /mean(TW)},
      error = function(e){
        NA
      })
    }
    
    if(!is.null(Y)){
    g$Y_error[ind] = tryCatch({
      (x$Y_coef*Q^x$Y_exp) %>% rmse(Y) /mean(T)},
      error = function(e){
        NA
      })
    }

    g
  }

  l = list()
  
  num = sum(!is.null(V), !is.null(TW), !is.null(Y))
  
  names = c("V", "TW", "Y")[c(!is.null(V), !is.null(TW), !is.null(Y))]
  
  
  for(i in 1:num){
    l[[i]] = v
  }
  
  g = expand.grid(l, stringsAsFactors = F) %>% 
    setNames(paste0(names, "_method")) %>% 
    mutate(c1 = NA, c2 = NA, viable = NA, tot_error = NA) 
  
  g = bind_cols(g, 
                setNames(data.frame(matrix(NA, ncol = 3*num, nrow =nrow(g))), c(paste0(names, "_error"), 
                                                                                paste0(names, "_coef"), 
                                                                                paste0(names, "_exp"))))

  
  
  types = c("TW", "Y", "V")
  
  for(t in 1:3){
    x = g[[paste0(types[t], "_method")]]
    
    ind = match(x, r[[types[t]]]$method)
    
    g[[paste0(types[t], '_exp')]] = r[[types[t]]]$exp[ind]
    g[[paste0(types[t], '_coef')]] = r[[types[t]]]$coef[ind]
  }
  
  for(i in 1:nrow(g)){
    g = fit(g, i, V, TW, Y, Q)
  }
  
    
    if(num == 3){
      g$c1 = round(g$V_coef * g$Y_coef * g$TW_coef, 3)
      g$c2 = round(g$V_exp + g$Y_exp + g$TW_exp, 3)
      
      g$viable =  (between(g$c1, 1-allowance, 1+allowance) +  between(g$c2, 1-allowance, 1+allowance)) == 2
    }
    
    g$tot_error = rowSums(g[, grepl("error", names(g))], na.rm = TRUE)

  
  if(num == 3){
    combo = g %>% 
      filter(viable == TRUE) %>% 
      arrange(tot_error) %>% 
      slice(1) %>% 
      mutate(condition = "bestValid")
    # 
    # combo =  g %>% 
    #   filter(apply(g,1,function(r){length(unique(r[1:3]))}) != 1) %>% 
    #   #filter(viable) %>% 
    #   arrange(tot_error) %>% 
    #   slice(1) %>%
    #   mutate(condition = "combo")
    
    ols = g %>% 
      filter(Y_method == "ols", TW_method == "ols", V_method == "ols") %>% 
      mutate(condition = "ols")
    
    nls = g %>% 
      filter(Y_method == "nls", TW_method == "nls", V_method == "nls") %>% 
      mutate(condition = "nls")
    
    if("nsga2" %in% v){
      ga = g %>% 
        filter(Y_method == "nsga2", TW_method == "nsga2", V_method == "nsga2") %>% 
        mutate(condition = "nsga2")
      
    } else {
      ga = NULL
  }
 
    return(list(g = g %>% 
                  arrange(!viable, tot_error), 
                summary = bind_rows(combo, ols, nls, ga) %>%  
                  dplyr::select(-c1, -c2) %>% 
                  arrange(!viable, tot_error)))
    
  } else {
    list(g = g %>% 
          arrange(tot_error), 
        summary = g %>%  
          dplyr::select(-c1, -c2, -viable) %>% 
          arrange(tot_error) %>% 
          slice(1))
  }
}
