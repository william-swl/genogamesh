#' parse the output of SingleR
#'
#' @param x DFrame object from SingleR::SingleR()
#'
#' @return tibble
#' @export
#'

parse_SingleR <- function(x) { # nolint
  if (!is(x, "DFrame")) {
    stop("x should be the output of SingleR::SingleR()")
  }

  res <- x %>%
    as_tibble(rownames = "cell") %>%
    pivot_longer(-c("cell", "labels", "delta.next", "pruned.labels"),
      names_to = c(".value", "celltype"), names_sep = "cores."
    )

  # celltype in colnames and labels have some differences in symbol
  res <- res %>% filter(
    baizer::reg_join(.data[["labels"]], "[\\da-zA-Z]") ==
      baizer::reg_join(.data[["celltype"]], "[\\da-zA-Z]")
  )
  res <- res %>%
    dplyr::select(-"celltype") %>%
    dplyr::select(dplyr::all_of(c("cell",
      "celltype" = "labels", "celltype_score" = "s",
      "delta_next" = "delta.next",
      "celltype_pruned" = "pruned.labels"
    )))
  return(res)
}
