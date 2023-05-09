#' write a sequence vector into fasta file
#'
#' @param x sequence vector
#' @param file fasta file
#'
#' @return NULL
#' @export
#'
write_fasta <- function(x, file) {
  if (is.null(names(x))) {
    names(x) <- seq_along(x)
  }

  res <- map2_chr(names(x), x, ~ str_c(">", .x, "\n", .y))
  out <- file(file, "w")
  writeLines(res, out)
  close(out)
}
