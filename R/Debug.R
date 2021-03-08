#' SummarizeHVF
#'
#' @param x Seurat object.
#'
#' @export
#'
SummarizeHVF <- function(x) {
  if (inherits(x, "Seurat")) x <- list(Seurat = x)
  total_datasets <- length(x)
  x <- lapply(names(x), function(n) {
    tibble(dataset = n, feature = VariableFeatures(x[[n]]))
  }) %>% bind_rows()
  total_features <- length(unique(x$feature))
  x %>% gather(dataset, feature) %>% count(feature, sort = TRUE, name = "n_datasets") %>% count(n_datasets, name = "n_features") %>% mutate(total_datasets = total_datasets, percentage_datasets = n_datasets / total_datasets, total_features = total_features, percentage_features = 100 * n_features / total_features, cumsum_features = rev(cumsum(rev(n_features))), cumsum_percentage_features = rev(cumsum(rev(percentage_features)))) %>% select(n_datasets, total_datasets, percentage_datasets, n_features, cumsum_features, total_features, percentage_features, cumsum_percentage_features)
}


#' TabulateHVF
#'
#' @param x Seurat object.
#'
#' @export
#'
TabulateHVF <- function(x) {
  if (inherits(x, "Seurat")) x <- list(Seurat = x)
  lapply(names(x), function(n) {
    tibble(dataset = n, feature = VariableFeatures(x[[n]]), n = 1)
  }) %>% bind_rows() %>%
    group_by(feature) %>% mutate(total = sum(n)) %>%
    spread(dataset, n, fill = 0)
}

#' IntersectFeatures
#'
#' @param x Seurat object.
#' @param assay assay name.
#' @param slot slot name.
#' @param use.hvf whether to use HVF.
#' @param cutoff cutoff for gene expression.
#'
#' @export
#'
IntersectFeatures <- function(x, assay = "RNA", slot = "data", use.hvf = TRUE, cutoff = 0) {
  features <- lapply(names(x), function(n) {
    if (use.hvf) {
      return(VariableFeatures(x[[n]]))
    } else {
      tmp <- GetAssayData(x[[n]], assay = assay, slot = slot)
      return(rownames(tmp)[rowSums(tmp) > cutoff] )
    }
  })
  Reduce(intersect, features)
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
  tmp <- tmp %>% select(integration, reference, batch, everything())
  tmp
}

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


#' GetCorrectionMatrix
#'
#' @param x CanekDebug object.
#' @param batch batch correction number/name.
#'
#' @export
#'
GetCorrectionMatrix <- function(x, batch = 1) {
  tmp <- GetDebugData(x)[[batch]][["Correction Data"]][["Correction Matrix"]]
  rownames(tmp) <- VariableFeatures(x)
  colnames(tmp) <- seq_len(ncol(tmp))
  tmp
}

#' GetNumberOfBatches
#'
#' @param x CanekDebug object.
#'
#' @export
#'
GetNumberOfBatches <- function(x) {
  length(GetDebugData(x)) - 1
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
  }) %>% bind_rows() %>% group_by(gene) %>% summarize(correction = sum(correction)) %>% arrange(desc(correction))
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


#' plotBatchMembeship
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
plotBatchMembeship <- function(x, batch = 1, add.ref = TRUE, plot.membership = TRUE, plot.mnn = FALSE, ...) {
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
