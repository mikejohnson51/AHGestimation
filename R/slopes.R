#' Extract Thalweg
#' From a data.frame of cross sections, a classified thalweg can be extracted as the connected LINESTRING
#' @param xs a data.frame containing cross sectional data. Required columns are hf_id, cs_id, X, Y, Z
#' @param crs the CRS of the XY coordinates
#' @return XYZ LINESTRING object
#' @export

extract_thalweg = function(xs, crs = 5070){
  xs |> 
    group_by(hf_id, cs_id) |> 
    slice_min(Z, with_ties = FALSE) |> 
    ungroup() |> 
    st_as_sf(coords = c("X", "Y", "Z"), crs = crs) |> 
    group_by(hf_id) |> 
    summarize(geometry = st_combine(geometry)) |> 
    st_cast("LINESTRING") 
}


#' @title Calculate the slope of 3D linestring
#' @description Given a sf object with 'XYZ' coordinates, return a vector
#' of numeric values representing the average slope of each linestring in the
#' sf data frame input.
#'
#' The default calculates the  slope using `slope_weighted()`.
#' You can also use `slope_mean()` or any other
#' function that takes the same inputs as these functions.
#' @param path an XYZ LINESTRING representing the path of travel
#' @param fun The slope function to calculate per element, `slope_weighted` is the default.
#' @param directed Should the value be directed? `FALSE` by default.
#'   If `TRUE` the result will be negative when it represents a downslope
#'   (when the end point is lower than the start point).
#' @return A vector of slopes associated with each linear element
#'   The value is a proportion representing the change in elevation
#'   for a given change in horizontal distance.
#' @export

compute_channel_slope = function(path, 
                                 fun = slope_weighted, 
                                 directed = FALSE) {
  
  if (inherits(path, "sf") | inherits(path, "sfc")) {
    lonlat = st_is_longlat(path)
    xyz    = as.data.frame(st_coordinates(path))
  }
  
  data.frame(
    hf_id = path[[1]],
    slope = slope_matrices(split(x = xyz, 
                      f = xyz[, "L1"]), 
                      lonlat = lonlat, 
                      fun = fun, 
                      directed = directed)
  )
}


#' @title Calculate the gradient of line segments from a 3D matrix of coordinates
#' @param mat  Matrix containing coordinates and elevations.
#'   The matrix should have three columns: X, Y, and Z. 
#'   In data with geographic coordinates, Z values are assumed to be in
#'   meters. In data with projected coordinates, Z values are assumed to have
#'   the same units as the X and Y coordinates.
#' @param lonlat Are the elements provided in longitude/latitude coordinates?
#'   By default, value is from the CRS of the routes (`sf::st_is_longlat(...)`).
#' @return A vector of slopes associated with each LINE element
#'   The output value is a proportion representing the change in elevation
#'   for a given change in horizontal distance.
#' @export

slope_matrix = function(mat, 
                        lonlat = TRUE) {
  
  if(lonlat) {
    d = geodist(mat[, c("X", "Y")], sequential = TRUE)
  } else {
    d = sqrt(diff(mat[, "X"])^2 + diff(mat[, "Y"])^2)
  }
  
  diff(mat[,"Z"]) / d
}

#' @rdname slope_matrix
#' @export

slope_weighted = function(mat,
                          lonlat = TRUE,
                          directed = FALSE) {
  if(lonlat) {
    d = geodist(mat[, c("X", "Y")], sequential = TRUE)
  } else {
    d = sqrt(diff(mat[, "X"])^2 + diff(mat[, "Y"])^2)
  }
  
  x = weighted.mean(abs(slope_matrix(mat, lonlat = lonlat)), d, na.rm = TRUE)
  
  if (directed) {
    x * sign(tail(mat[,"Z"], 1) - head(mat[,"Z"], 1))
  } else {
    x
  }
}

#' @rdname slope_matrix
#' @export

slope_mean = function(mat, 
                      lonlat = TRUE, 
                      directed = FALSE) {
  
  x = mean(abs(slope_matrix(mat, lonlat = lonlat)), na.rm = TRUE)
  
  if(directed) {
     x * sign(tail(mat[,"Z"], 1) - head(mat[, "Z"], 1))
  } else {
    x
  }
}

slope_matrices = function(mat_split, fun = slope_matrix_weighted, ...) {
  unlist(pblapply(mat_split, fun, ...))
}
