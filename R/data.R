#' Sample gage data
#' Manual measurements made at NWIS site 01096500
#' Q_cms is a mandatory argument and at least one of TW_m, V_ms, or Y_m. 
#' @format
#' A data frame with 245 rows and 6 columns:
#' \describe{
#'   \item{siteID}{NWIS ID}
#'   \item{date}{date of measurement}
#'   \item{Q_cms}{Steamflow (cubic meters per second)}
#'   \item{Y_m}{Depth (meters)}
#'   \item{V_ms}{Velocity (meters per second)}
#'   \item{TW_m}{Top width (meters)}
#' }
#' @family data

"nwis"
