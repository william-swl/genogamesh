% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_bio.R
\name{parse_ANARCI_aaseq}
\alias{parse_ANARCI_aaseq}
\title{parse sequences from ANARCI vdj output}
\usage{
parse_ANARCI_aaseq(
  x,
  chain,
  remove_gap = TRUE,
  scheme = "imgt",
  number_table = NULL,
  keep_number = FALSE
)
}
\arguments{
\item{x}{tibble from anarci_H.csv or anarci_KL.csv}

\item{chain}{one of 'H, L', for heavy chain or light chain}

\item{remove_gap}{remove the gap caused by numbering, default as TRUE}

\item{scheme}{antibody numbering scheme, one of 'imgt, chothia'}

\item{number_table}{custom antibody numbering system, please input a tibble
with three columns: region, start, end.}
}
\value{
tibble
}
\description{
parse sequences from ANARCI vdj output
}
