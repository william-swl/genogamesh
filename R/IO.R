#' read fasta file
#'
#' @param x file path
#' @param type one of `DNA, RNA, AA`
#'
#' @return XXXStringSet object
#' @export
#'
read_fasta <- function(x, type = "DNA") {
  if (type == "DNA") {
    Biostrings::readDNAStringSet(x, format = "fasta")
  } else if (type == "RNA") {
    Biostrings::readRNAStringSet(x, format = "fasta")
  } else if (type == "AA") {
    Biostrings::readAAStringSet(x, format = "fasta")
  }
}

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


#' read vcf file
#'
#' @param x file path
#'
#' @return vcfR object
#' @export
#'
read_vcf <- function(x) {
  vcfR::read.vcfR(x, verbose = FALSE)
}

#' read gff file
#'
#' @param x file path
#'
#' @return tibble
#' @export
#'
read_gff <- function(x) {
  readr::read_tsv(x,
    col_names = c(
      "seqid", "source", "type", "start", "end",
      "score", "strand", "phase", "attribute"
    )
  )
}
