#' @title Report best optimal
#' @param best best performing method (character string)
#' @param check values to check against
#' @return
#' @export
#' @examples

best_optimal = function(best, check){
  
  message(paste0("The best performing method (", best,") ", ifelse(check, "is", "is not"), " physically valid ... ", ifelse(check, emo::ji("joy"), emo::ji("sad")), "\n")
  )
  
  return(check)
}

#' @title Properly estimate AHG values
#' @param Q streamflow timeseries
#' @param Y streamflow timeseries
#' @param V streamflow timeseries
#' @param TW streamflow timeseries
#' @param allowance streamflow timeseries
#' @param quiet streamflow timeseries
#' @param GA use GA or NSGA2?
#' @return
#' @export

ahg_estimate = function(Q, Y = NULL, V = NULL, TW = NULL, allowance = .05, quiet = FALSE, GA = T){
  
  ahg_y  = if(!is.null(Y)){ compute_ahg(Q,Y, "Y") }
  ahg_tw = if(!is.null(TW)){ compute_ahg(Q,TW, "TW") }
  ahg_v  = if(!is.null(V)){ compute_ahg(Q,V, "V") }
  
  r = c(ahg_y, ahg_tw, ahg_v) 
  
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

    cond = best_optimal(best, ifelse(best[1] == "nls", nls_viable, ols_viable))
    
    if(!cond){ 
      message("Launching Evolutionary Algorithm... allownace (", allowance,") ", emo::ji("biceps"))
      if(GA){
        r = calc_ga_2(Q, Y, V, TW, allowance, r)
      } else {
        r = calc_ga(Q, Y, V, TW, allowance, r)
      }
      
      m = mismash(v = c("ols", "nls", "GA"), params= length(r),V, TW, Y, Q, r, allowance )
    } else{
      m = mismash(v = c("ols", "nls"), params = length(r),V, TW, Y, Q, r, allowance)
    }
    
    }else if(length(r) == 2){
      m = mismash(v = c("ols", "nls"), params = length(r),V, TW, Y, Q,r, allowance)
    } else {
      m = NULL
    }
  
    
    return(c(r, m)) 

}

  