#' GetCorrectionMatrix
#'
#' @param x CanekDebug object.
#' @param batch batch correction number/name.
#'
#' @export
#'
GetCorrectionMatrix <- function(x, batch = 1) {
  tmp <- GetDebugData(x)[[batch]][["debug"]][["matrix"]][["corGene"]]
  rownames(tmp) <- GetDebugData(x)[[batch]][["debug"]][["matrix"]][["features"]]
  colnames(tmp) <- names(GetDebugData(x)[[batch]][["debug"]][["membership"]])
  tmp
}

#' SummarizeTopGenes
#'
#' @param x CanekDebug object.
#' @param batch batch correction number/name.
#'
#' @export
#'
SummarizeTopGenes <- function(x, batch = NULL) {
  tmp <- GetDebugData(x)
  tmp[["Total_Correction_Time"]] <- NULL
  if (!is.null(batch)) tmp <- tmp[batch]
  lapply(names(tmp), function(n) {
    enframe(rowSums(GetCorrectionMatrix(x, n) ^ 2), name = "gene", value = "correction") %>% add_column(batch = basename(n))
  }) %>% bind_rows() %>% group_by(.data[["gene"]]) %>% summarize(correction = sum(.data[["correction"]])) %>% arrange(desc(.data[["correction"]]))
}
