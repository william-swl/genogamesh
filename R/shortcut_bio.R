#' add SingleR celltype annotation for Seurat object
#'
#' @param x Seurat object
#' @param ref SummarizedExperiment object, the annotation reference of SingleR
#'
#' @return Seurat object
#' @export
#'
SingleR_SE <- function(x, ref) { # nolint
  if (!is(x, "Seurat")) {
    stop("x should be a Seurat object")
  }

  if (!is(ref, "SummarizedExperiment")) {
    stop("ref should be a SummarizedExperiment object")
  }

  TBref_celltype <- SummarizedExperiment::colData(ref) %>% # nolint
    as_tibble() %>%
    dplyr::distinct(.data[["label.fine"]], .keep_all = TRUE)
  fine2main <- TBref_celltype %>%
    dplyr::pull(.data[["label.main"]], .data[["label.fine"]])

  # call SingleR
  TBannot <- SingleR::SingleR( # nolint
    test = Seurat::as.SingleCellExperiment(x),
    ref = ref, labels = ref$label.fine
  ) %>%
    parse_SingleR()
  # add main celltype
  TBannot <- TBannot %>% dplyr::mutate( # nolint
    celltype_main = fine2main[.data[["celltype"]]],
    .after = .data[["celltype"]]
  )
  # add meta.data
  res <- Seurat::AddMetaData(x, TBannot %>% baizer::c2r("cell"))

  return(res)
}



#' reduction from raw Seurat object created by read count matrix, including
#' normalization, variable features calling, scaling, PCA and UMAP
#'
#' @param x Seurat object
#' @param use_dim dims to use in UMAP reduction
#'
#' @return Seurat object
#' @export
#'
reduction_SE <- function(x, use_dim = 30) { # nolint
  # Normalizing the data
  x <- Seurat::NormalizeData(x,
    normalization.method = "LogNormalize", scale.factor = 10000
  )
  # Identification of highly variable features (feature selection)
  x <- Seurat::FindVariableFeatures(x,
    selection.method = "vst", nfeatures = 2000
  )
  # Scaling the data
  x <- Seurat::ScaleData(x)
  # Dimension reduction
  x <- Seurat::RunPCA(x,
    features = Seurat::VariableFeatures(x), verbose = FALSE
  )
  x <- Seurat::RunUMAP(x, dims = 1:use_dim)

  res <- Seurat::AddMetaData(
    x,
    as.data.frame(x@reductions$umap@cell.embeddings)
  )

  return(res)
}
