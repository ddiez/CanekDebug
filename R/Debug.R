#' GetDebugData
#'
#' @param x CanekDebug object.
#'
#' @export
#'
GetDebugData <- function(x) {
  if (inherits(x, "Seurat"))
    Tool(x, "RunCanekSeurat")
  else
    x[-2]
}

#' GetDebugInfo
#'
#' @param x CanekDebug object.
#'
#' @export
#'
GetDebugInfo <- function(x) {
  x <- GetDebugData(x)
  tt <- x[["Total_Correction_Time"]]
  x[["Total_Correction_Time"]] <- NULL
  tmp <- enframe(sapply(x, "[[", "Correction_Time"), name = "integration", value = "time")
  tmp$units <- sapply(lapply(x, "[[", "Correction_Time"), units)
  tmp$memberships <- sapply(x, function(xx) xx[["Correction Data"]][["Membership Data"]][["Cluster Membership"]])
  tmp$reference <- unique(sub("/.*", "", tmp$integration))
  tmp$batch <- factor(basename(tmp$integration), levels = basename(tmp$integration))
  tmp <- tmp %>% select(.data[["integration"]], .data[["reference"]], .data[["batch"]], everything())
  tmp
}
