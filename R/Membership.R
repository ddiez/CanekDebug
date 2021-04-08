#' GetMembershipData
#'
#' @param x CanekDebug object.
#' @param batch batch correction number/name.
#'
#' @export
#'
GetMembershipData <- function(x, batch = 1) {
  tmp <- GetDebugData(x)
  #tmp[[batch]][["Correction Data"]][["Membership Data"]][["Membership Correction Data"]]
  tmp[[batch]][["debug"]][["membership"]]
}

#' CheckUncorrectedMemberships
#'
#' @param x CanekDebug object.
#'
#' @export
#'
CheckUncorrectedMemberships <- function(x) {
  x <- GetDebugData(x)
  x[["Total_Correction_Time"]] <- NULL
  bind_rows(lapply(x, function(xx) bind_rows(xx$debug$membership)), .id = "batch") %>% count(.data[["batch"]], .data[["membership"]], .data[["corrected"]]) %>% filter(!.data[["corrected"]]) %>% select(-.data[["corrected"]])
}

#' GetMembershipIndex
#'
#' @param x CanekDebug object.
#' @param batch batch correction number/name.
#'
#' @export
#'
GetMembershipIndex <- function(x, batch = 1) {
  tmp <- GetMembershipData(x, batch)

  # tmp <- lapply(names(tmp), function(membership) {
  #   data.frame(membership = membership, index = tmp[[membership]][["Cells Index"]])
  # }) %>% bind_rows() %>% arrange(index) %>% select(-index)
  # rownames(tmp) <- GetBatchCellNames(x, batch)

  tmp <- tmp %>% bind_rows() %>% column_to_rownames("cells")
  tmp
}
