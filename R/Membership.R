#' CheckUncorrectedMemberships
#'
#' @param x CanekDebug object.
#'
#' @export
#'
CheckUncorrectedMemberships <- function(x) {
  x <- GetDebugData(x)
  x[["Total_Correction_Time"]] <- NULL
  bind_rows(lapply(x, function(xx) bind_rows(xx$debug$membership)), .id = "batch") %>% count(batch, membership, corrected) %>% filter(!corrected) %>% select(-corrected)
}
