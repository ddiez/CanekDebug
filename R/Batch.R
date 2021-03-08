#' GetNumberOfBatches
#'
#' @param x CanekDebug object.
#'
#' @export
#'
GetNumberOfBatches <- function(x) {
  length(GetDebugData(x)) - 1
}

#' GetRefBatchName
#'
#' @param x CanekDebug object.
#'
#' @export
#'
GetRefBatchName <- function(x) {
  tmp <- GetDebugData(x)
  sub("/.*", "", names(tmp)[[1]])
}


#' GetBatchNames
#'
#' @param x CanekDebug object.
#'
#' @export
#'
GetBatchNames <- function(x) {
  rev(rev(names(GetDebugData(x)))[-1])
}

#' GetBatchCellNames
#'
#' @param x CanekDebug object.
#' @param batch batch correction number/name.
#'
#' @export
#'
GetBatchCellNames <- function(x, batch = 1) {
  tmp <- GetDebugData(x)
  #ref <- GetRefBatchName(x)
  if (is.numeric(batch))
    query <- GetBatchNames(x)[batch]
  colnames(tmp[[batch]][["Query Batch (B2)"]])
}

#' GetRefCellNames
#'
#' @param x CanekDebug object.
#'
#' @export
#'
GetRefCellNames <- function(x) {
  tmp <- GetDebugData(x)
  colnames(tmp[[1]][["Reference Batch (B1)"]])
}
