#' Implements filtering by date
#' @description Data is filtered when it is beyond a specified year threshold (e.g. 5 years old). The relative date is based on the newest 
#' observation in the data set. Optionally, the maximum flow (Q) record can be retained.
#' @param df a data.frame with at least a date and Q field. 
#' @param years the number of allowed history
#' @param keep_max Should the largest flow record be kept, even if older then "years"
#' @return data.frame
#' @family filters
#' @export

date_filter = function(df, years, keep_max = FALSE){

  Q <-  NULL
  if(!is.null(df$date)){
    y = format(max(df$date), "%Y")
    y = as.numeric(y) - years
    fdate = as.Date(paste0(y, "-01-01"))
    
    if(keep_max){
      df2 = filter(df, date >= fdate | Q == max(df$Q))
    } else {
      df2 = filter(df, date >= fdate)
    }
  
    if(nrow(df2) < 10){
      warning("Too much data was filtered based on date. Input data returned", call. = FALSE)
    } else {
      df = df2
    }
  } else {
    warning("date column not provided", call. = FALSE)
  }
  
  df
}

#' Implements filtering by continuity
#' @description The function tests if the measured Q is outside of the expected range based on the product 
#' of measured velocity, top-width, and depth (e.g. Q≠vA) 
#' @param df a data.frame with a Q, Y, TW, V and field. 
#' @param allowance how much deviation from equality should be allowed (default = .05)
#' @return data.frame
#' @family filters
#' @export


qva_filter = function(df, allowance = .05) {
  
  Qt <- NULL
  if (!all(c("TW", "V", "Q", "Y") %in% names(df))) {
    warning("Q, TW, V, Y are all needed to run this filter.")
    df
  } else {
    df$Qt = df$TW * df$Y * df$V
    df2 = df[between(df$Qt, (1 - allowance) * df$Q, (1 + allowance) * df$Q),]
    
    if (nrow(df2) < 10) {
      warning("Too much data was filtered using Q = vA metric. Input data returned",
              call. = FALSE)
    } else {
      df = df2
    }
    
    select(df, -Qt)
  }
}

#' Implements filtering by median absolute deviation
#' @description An iterative outlier detection procedure is run based on to the linear regression residuals. 
#' Values of log-transformed TW, V, and Y residuals falling outside a specified median absolute deviation (MAD) envelope are excluded. 
#' Regression coefficients were recalculated and the outlier detection procedure was reapplied until no outliers are detected. This method
#'  was identified in \href{https://zenodo.org/record/7868764}{HyG}
#' @param df a data.frame with at least a Q and one other AHG field (Y. TW, V). 
#' @param envelope MAD envelope
#' @return data.frame
#' @family filters
#' @export

mad_filter = function(df, envelope = 3){
  
  mad_f = function(df, relation, envelope){
    
    if(relation %in% names(df)){
      num = 100
      
      while(num != 0){
        m = lm(log(df$Q) ~ log(df[[relation]]))
        s = summary(m)
        ind = s$residuals > (envelope*mad(s$residuals))
        df = df[!ind, ]
        
        num = sum(ind)
      }
      
      df
    } else {
      df
    }
   
  }
  
  df2 = mad_f(df, "Y", envelope)
  df2 = mad_f(df2, "TW", envelope)
  df2 = mad_f(df2, "V", envelope)
  
  if(nrow(df2) < 10){
    warning("Too much data was filtered using MAD. Less then 10 obs remain. Keeping complete history.", call. = FALSE)
  } else {
    df = df2
  }
  df
}

#' Implements NLS filtering
#' @description An NLS fit provides the best relation by relation fit. For each provided relationship, an NLS fit is computed and used to 
#' estimate the predicted {V,TW,Y} for a given Q. If the actual value is outside the specified allowance it is removed.
#' @param df a data.frame with at least a Q and one other AHG field (Y. TW, V). 
#' @param allowance how much deviation from observed should be allowed (default = .5)
#' @return data.frame
#' @family filters
#' @export

nls_filter = function(df, allowance = .5){
  
  method <- NULL
  
  predict_nls_fhg = function(Q, Y) {
    if(is.null(Y)){
      NULL
    } else {
      o = compute_fhg(Q, Y, "zzz") %>% 
        filter(method == "nls")
      o$coef * (Q ^ o$exp)
    }
  }
  
  n = nrow(df)
  
  df2 = df
  df2$pY  = predict_nls_fhg(Q = df2$Q, Y = df2$Y)
  df2$pTW = predict_nls_fhg(df2$Q, df2$TW)
  df2$pV  = predict_nls_fhg(df2$Q, df2$V)
  
  if(!is.null(df2$pY)){
    df2 = df2 %>%
      filter(between(
        df2$Y,
        (1 - allowance) * df2$pY,
        (1 + allowance) * df2$pY
      ))
  }
  
  if(!is.null(df2$pV)){
    df2 = df2 %>%
      filter(between(
        df2$V,
        (1 - allowance) * df2$pV,
        (1 + allowance) * df2$pV
      ))
  }
  
  if(!is.null(df2$pTW)){
    df2 = df2 %>%
      filter(between(
        df2$TW,
        (1 - allowance) * df2$pTW,
        (1 + allowance) * df2$pTW
      ))
  }
  
  
  if(nrow(df2) < 10){
    warning("Too much data was filtered using NLS window. Less then 10 obs remain. Keeping all data.", call. = FALSE)
  } else {
    df = df2
  }
  
  select(df, -any_of(c("pTW", "pV", "pY")))
}

#' Implements significance check
#' @description The relationship between all supplied log transformed variables are computed. 
#' If the p-value of any of these is less then the supplied p-value an error message is emitted.
#' @param df a data.frame with at least a Q and one other AHG field (Y. TW, V). 
#' @param pvalue Significant p-value (default = .05)
#' @return data.frame
#' @family filters
#' @export

significance_check = function(df, pvalue = .05){
  
  sig_f = function(df, relation, pvalue){
    
    if(relation %in% names(df)){
      m = lm(log(df$Q) ~ log(df[[relation]]))
      f = summary(m)$fstatistic
      p <- pf(f[1],f[2],f[3],lower.tail=F)
      attributes(p) <- NULL
      p < pvalue
    } else {
      TRUE
    }
  }

  x1 = sig_f(df, "Y", pvalue)
  x2 = sig_f(df, "TW", pvalue)
  x3 = sig_f(df, "V", pvalue)
  
  if(any(!x1,!x2,!x3)){
    stop("Significance missing for relation: ", c("Y ", "TW ", "V ")[which(c(!x1,!x2,!x3))], call. = FALSE)
  }
  
  df
  
}
