#' GetMNNPairsData
#'
#' @param x CanekDebug object.
#' @param batch batch correction number/name.
#'
#' @export
#'
GetMNNPairsData <- function(x, batch = 1) {
  tmp <- GetDebugData(x)
  #tmp[[batch]][["Correction Data"]][["MNN Pairs"]]
  tmp[[batch]][["debug"]][["pairs"]]
}
