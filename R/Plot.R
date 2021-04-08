#' plotRunningTime
#'
#' @param x CanekDebug object.
#'
#' @export
#'
plotRunningTime <- function(x) {
  d <- GetDebugInfo(x)
  tt <- GetTotalTime(x)
  ggplot(d, aes(batch, time)) +
    geom_col() +
    labs(y = paste0("time (", d$units[1], ")"), title = paste0("reference: ", d$reference[1], ", total time: ", round(tt, digits = 2), " (", units(tt), ")"))
}

#' plotTopGenes
#'
#' @param x CanekDebug object.
#' @param batch batch correction number/name.
#'
#' @export
#'
plotTopGenes <- function(x, batch = NULL) {
  d <- SummarizeTopGenes(x, batch = batch) %>% mutate(index = seq_along(correction))
  ggplot(d, aes(index, correction)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = gene), data = d %>% top_n(10, correction), min.segment.length = 0, color = "red")
}

#' plotCorrectionHist
#'
#' @param x CanekDebug object.
#' @param batch batch correction number/name.
#' @param by.membership whether facet plot by membership.
#'
#' @export
#'
plotCorrectionHist <- function(x, batch = 1, by.membership = TRUE) {
  tmp <- GetCorrectionMatrix(x, batch) %>% as_tibble()
  p <- tmp %>% gather(membership, correction) %>%
    ggplot(aes(correction)) +
    geom_histogram(binwidth = .1)
  if (by.membership)
    p <- p + facet_wrap(~membership)
  p + labs(title = paste0("batch: ", batch))
}

#' plotCorrectionHeatmap
#'
#' @param x CanekDebug object.
#' @param batch batch correction number/name.
#'
#' @export
#'
plotCorrectionHeatmap <- function(x, batch = 1) {
  tmp <- GetCorrectionMatrix(x, batch)
  Heatmap(tmp[rowSums(tmp) != 0, ], cluster_columns = FALSE, name = "correction", column_title = paste0("batch: ", batch))
}

#' plotNumberOfCorrectedGenes
#'
#' @param x CanekDebug object.
#'
#' @export
#'
plotNumberOfCorrectedGenes <- function(x) {
  batches <- GetBatchNames(x)
  d <- lapply(batches, function(batch) {
    tmp <- GetCorrectionMatrix(x, batch)
    data.frame(batch = basename(batch), genes = sum(rowSums(tmp) != 0))
  }) %>% bind_rows()
  d$batch <- factor(d$batch, levels = basename(batches))
  ggplot(d, aes(batch, genes)) +
    geom_col()
}

#' plotBatchMembership
#'
#' @param x CanekDebug object.
#' @param batch batch correction number/name.
#' @param add.ref whether to plot reference cells.
#' @param plot.membership whether to color by memberships.
#' @param plot.mnn whether to add MNN pairs.
#' @param ... arguments passed down to methods.
#'
#' @export
#'
plotBatchMembership <- function(x, batch = 1, reduction = "umapraw", add.ref = TRUE, plot.membership = TRUE, plot.mnn = FALSE, ...) {
  tmp <- GetDebugData(x)
  ref <- GetRefBatchName(x)
  if (is.numeric(batch))
    batch <- GetBatchNames(x)[batch]
  batchcells <- colnames(tmp[[batch]][["Query Batch (B2)"]])
  cells <- batchcells
  if (add.ref) {
    refcells <- colnames(tmp[[batch]][["Reference Batch (B1)"]])
    cells <- c(refcells, batchcells)
  }
  x <- x[, cells]

  #d <- data.frame(], batchcells[mnn[["queBatch-Cells-Index"]]])
  #g <- igraph::graph_from_data_frame(d)

  if (plot.membership) {
    membership <- GetMembershipIndex(x, batch)
    x <- AddMetaData(x, membership)
  } else {
    x$membership <- "Reference"
    x$membership[batchcells] <- batch
    x$membership <- factor(x$membership, levels = c("Reference", batch))
  }

  p <- DimPlot(x, group.by = "membership", reduction = "umapraw", ...)

  if (plot.mnn) {
    mnn <- GetMNNPairsData(x, batch)
    emb <- Embeddings(x, "umapraw")
    # ref_mnn <- refcells[mnn[, "refBatch-Cells-Index"]]
    # query_mnn <- batchcells[mnn[, "queBatch-Cells-Index"]]

    ref_mnn <- mnn[, "ref"]
    query_mnn <- mnn[, "query"]

    d <- as.data.frame(cbind(emb[ref_mnn, ], emb[query_mnn, ]))
    colnames(d) <- c("x", "y", "xend", "yend")

    if (plot.membership) {
      d$membership <- factor(membership[query_mnn, "membership"])
      p <- p + geom_segment(aes(x, y, xend = xend, yend = yend, color = membership), data = d, size = .1, alpha = .4)
    } else {
      p <- p + geom_segment(aes(x, y, xend = xend, yend = yend), data = d, size = .1, color = "grey", alpha = .4)
    }
  }
  p + labs(title = batch)
}
