#' parse substitutions and deletions of amino acids
#'
#' @param mut_vector mutation vector
#'
#' @return list of mutation table
#' @export
#'
#' @examples
#' mut <- "G1:T10I,G2:D20N,G3:Q30E,G1:A40T,G3:P50L,G2:G60R"
#' parse_aa_mut(mut)
#'
parse_aa_mut <- function(mut_vector) {
  pattern <- "(\\w+):(\\w)(\\d+)([\\w-])"

  mut_vector %>%
    stringr::str_match_all(pattern) %>%
    purrr::map(
      ~ tibble::as_tibble(.x) %>%
        rlang::set_names(c("full", "gene", "ref", "pos", "alt"))
    )
}
