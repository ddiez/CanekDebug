#' GetTotalTime
#'
#' @param x CanekDebug object.
#'
#' @export
#'
GetTotalTime <- function(x) {
  x <- GetDebugData(x)
  x[["Total_Correction_Time"]]
}
