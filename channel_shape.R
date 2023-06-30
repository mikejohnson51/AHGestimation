#' @title Approximate channel shape
#' @description Get a list of points from x axis of a cross section and max depth and produce depth values for those points based on channel shape
#' @param r The corresponding dingman's r coeficent 
#' @param TW width of the channel at bankfull 
#' @param Ymax maximum depth of the channel at bankfull 
#' @return depth values every 1m along the cross section   
#' @export

cross_section = function(r, TW = 30, Ymax = 2){
  TW = as.integer(TW)
  half = TW/2
  x_list <- list()
  z_list <- list()
  print(TW)
  for(point in -half:half){
    print(point)
    Z = Ymax * (2/TW)**r * (point**r)
    x_list <- append(x_list, point)
    z_list <- append(z_list, Z)
  }
  nested_list <- list(X=x_list, Z=z_list)
  df <- as.data.frame(do.call(cbind, nested_list))
  return (df)
}