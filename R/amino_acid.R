#' return the index in the correspond codon (1, 2 or 3)
#'
#' @param x 1-based nt index
#' @param start ORF start
#' @param end ORF end
#'
#' @return 1, 2 or 3
#' @export
#'
#' @examples codon_index(5, 1, 12)
codon_index <- function(x, start, end) {
  if (any((end - start) %% 3 != 2)) {
    stop("pleas input 1-based close range for start and end")
  }

  (x - start) %% 3 + 1
}

#' return the codon range of the correspond codon
#'
#' @param x 1-based nt index
#' @param start ORF start
#' @param end ORF end
#'
#' @return codon range string
#' @export
#'
#' @examples codon_range(5, 1, 12)
codon_range <- function(x, start, end) {
  idx <- codon_index(x, start, end)
  str_c(x - idx + 1, x + 3 - idx, sep = "-")
}

#' return the amino acid index in this ORF
#'
#' @param x 1-based nt index
#' @param start ORF start
#' @param end ORF end
#'
#' @return amino acid index
#' @export
#'
#' @examples aa_index(5, 1, 12)
aa_index <- function(x, start, end) {
  floor((x - start) / 3) + 1
}
