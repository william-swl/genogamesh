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




#' translate nucleotides into amino acids from the first character
#'
#' @param x nucleotides characters
#'
#' @return amino acids characters
#' @export
#'
#' @examples nt2aa(c("ATGAAA", "TTGCCC", "CTGTTT"))
nt2aa <- function(x) {
  ds <- Biostrings::DNAStringSet(x)
  res <- Biostrings::translate(ds,
    no.init.codon = TRUE, if.fuzzy.codon = "solve"
  ) %>%
    as.character()
  return(res)
}


#' build antigen map from sera titer data
#'
#' @param data titer tibble, colnames should be `id, antigen1, antigen2, ...`
#' @param sera_meta metadata of sera, colnames should be `id, arg1, arg2, ...`
#' @param ag_meta metadata of antigen, colnames should be `id, arg1, arg2, ...`
#' @param value_lim titer value limitation, default as `c(0, Inf)`
#' @param n_optim optimization runs to perform
#' @param seed random seed
#'
#' @return tibble
#' @export
#'
antigen_map <- function(data, sera_meta = NULL, ag_meta = NULL,
                        value_lim = c(0, Inf),
                        n_optim = 500, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }


  # prepare data
  data <- data %>%
    baizer::c2r("id") %>%
    mutate_all(~ ifelse(.x < value_lim[1], value_lim[1], .x)) %>%
    mutate_all(~ ifelse(.x > value_lim[2], value_lim[2], .x)) %>%
    mutate_all(as.character) %>%
    mutate_all(~ ifelse(is.na(.x), "*", .x)) %>%
    mutate_all(~ ifelse(.x == as.character(value_lim[1]),
      str_c("<", value_lim[1]), .x
    )) %>%
    mutate_all(~ ifelse(.x == as.character(value_lim[2]),
      str_c(">", value_lim[2]), .x
    )) %>%
    t()

  # construct map
  map <- Racmas::acmap(titer_table = data) %>%
    Racmas::optimizeMap(
      map = .,
      number_of_dimensions = 2,
      number_of_optimizations = n_optim,
      minimum_column_basis = "none",
      options = list(ignore_disconnected = TRUE)
    )

  # extract coords
  ag_data <- Racmas::agCoords(map, 1) %>%
    as_tibble(rownames = "id") %>%
    dplyr::mutate(datatype = "ag")
  sera_data <- Racmas::srCoords(map, 1) %>%
    as_tibble(rownames = "id") %>%
    dplyr::mutate(datatype = "sera")
  if (!is.null(sera_meta)) {
    sera_data <- sera_data %>% left_join(sera_meta, by = "id")
  }
  if (!is.null(ag_meta)) {
    ag_data <- ag_data %>% left_join(ag_meta, by = "id")
  }

  return(bind_rows(sera_data, ag_data))
}
