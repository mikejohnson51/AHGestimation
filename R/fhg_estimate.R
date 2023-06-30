#' @title Report best optimal
#' @param best best performing method (character string)
#' @param check values to check against
#' @return vector
#' @export

best_optimal = function(best, check){
  
  message(paste0("The best performing method (", best,") ", 
                 ifelse(check, "is", "is not"), 
                 " physically valid ... ", 
                 ifelse(check, emo::ji("joy"), emo::ji("sad")), "\n")
  )
  
  return(check)
}

#' @title Properly estimate FHG values
#' @param Q streamflow time series
#' @param Y depth time series
#' @param V velocity time series
#' @param TW top width time series
#' @param allowance allowed deviation from continuity
#' @param quiet should messages be emitted?
#' @param forceGA use GA over nsga2
#' @param nwis NWIS gage ID
#' @param location provide a locations
#' @return list
#' @export

fhg_estimate = function(Q, Y = NULL, V = NULL, TW = NULL, 
                        allowance = .05, quiet = FALSE, forceGA= FALSE,
                        nwis = NULL, location = NULL){
  
  fhg_y  = if(!is.null(Y)){ compute_fhg(Q,Y, "Y") }
  fhg_tw = if(!is.null(TW)){ compute_fhg(Q,TW, "TW") }
  fhg_v  = if(!is.null(V)){ compute_fhg(Q,V, "V") }
  
  r = c(fhg_y, fhg_tw, fhg_v) 
  
  best = NULL
  
  for(i in 1:length(r)){ best[i] = r[[i]]$method[which.min(r[[i]]$nrmse)]}

  if(!quiet){ message(paste0(unique(toupper(best)), " performs best for the Q-", names(r), ' realtionship\n')) }
    
  if(length(r) == 3){
    
    is_between = function(x, allowance) { between(x, 1-allowance, 1+allowance)}
    
    c1_nls = is_between(prod(r$Y$coef[2], r$TW$coef[2], r$V$coef[2] ), allowance)
    c1_ols = is_between(prod(r$Y$coef[1], r$TW$coef[1], r$V$coef[1] ), allowance)
    
    c2_nls = is_between(sum(r$Y$exp[2], r$TW$exp[2], r$V$exp[2] ), allowance)
    c2_ols = is_between(sum(r$Y$exp[1], r$TW$exp[1], r$V$exp[1] ), allowance)
    
    nls_viable = sum(c1_nls, c2_nls) == 2
    ols_viable = sum(c1_ols, c2_ols) == 2
    
    message(paste("NLS", ifelse(nls_viable, "meets", "does not meet"), "continuity ... ",
                  ifelse(nls_viable, emo::ji("+1"), emo::ji("-1"))
                  ))
    message(paste("OLS", ifelse(ols_viable, "meets", "does not meet"), "continuity ... ",
                  ifelse(ols_viable, emo::ji("+1"), emo::ji("-1")), "\n"))

    best = unique(best)
    
    cond = best_optimal(best, ifelse(best[1] == "nls", nls_viable, ols_viable))
  
    if(any(!cond) | forceGA){ 
      message("Launching Evolutionary Algorithm... allowance (", allowance,") ", emo::ji("biceps"))
      
      r = calc_ga(Q, Y, V, TW, allowance, r)
      m = mismash(v = c("ols", "nls", "GA"), V, TW, Y, Q, r, allowance)
      
    } else{
      m = mismash(v = c("ols", "nls"), V, TW, Y, Q, r, allowance)
    }
    
    } else if(length(r) == 2){
      m = mismash(v = c("ols", "nls"), V, TW, Y, Q, r, allowance)
    } else {
      m = NULL
    }
  
  r = bind_rows(r) %>%  arrange(nrmse)

  if(!is.null(nwis)) {
    r$comid = dataRetrieval::findNLDI(nwis = nwis)$comid
  } 
  
  if(!is.null(location)) {
    r$comid = dataRetrieval::findNLDI(location = location)$comid
  } 
  
    
  # if(!is.null(m)){
  #   tmp = filter(r, type == "Y")
  #   m$summary$c =  tmp$coef[match(m$summary$Y, tmp$method)]
  #   m$summary$f =  tmp$exp[match(m$summary$Y, tmp$method)]
  # 
  #   tmp = filter(r, type == "TW")
  #   m$summary$a =  tmp$coef[match(m$summary$TW, tmp$method)]
  #   m$summary$b =  tmp$exp[match(m$summary$TW, tmp$method)]
  #   
  #   tmp = filter(r, type == "V")
  #   m$summary$k =  tmp$coef[match(m$summary$V, tmp$method)]
  #   m$summary$m =  tmp$exp[match(m$summary$V, tmp$method)]
  #   m$summary$comid = r$comid[1]
  #   m$summary$r = m$summary$f / m$summary$b
  # }
      
  out = list(summary = r, output = m$summary)
  out = Filter(Negate(is.null), out)
  if(length(out) == 1){ out = out[[1]]}
  
  return(out) 

}

  