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


#' parse somatic hypermutation from igblast output
#'
#' @param igblast_fmt7 the filename of igblast result
#'
#' @return tibble
#' @export
#'
parse_IgBlast_shm <- function(igblast_fmt7) {
  alltxt <- readLines(igblast_fmt7)
  Ltxt <- baizer::split_vector(alltxt,
    grep("# IGBLASTN", alltxt),
    bounds = "[)"
  )
  Ltxt <- Ltxt[-1]

  parse_mut <- function(x) {
    contig_id <- str_subset(x, "# Query:") %>%
      str_replace("# Query: ", "")

    # VDJ, query id, subject id, % identity, alignment length,
    # mismatches, gap opens, gaps...
    vdj <- c(
      x %>% str_subset("^V\\t") %>% .[1],
      x %>% str_subset("^D\\t") %>% .[1],
      x %>% str_subset("^J\\t") %>% .[1]
    )
    vdj_name <- c("V", "D", "J")
    vdj_mut <- vdj %>%
      str_split("\t") %>%
      map_chr(~ if (all(is.na(.x))) NA_character_ else .x[6])

    # region, from, to, length, matches, mismatches, gaps, percent identity
    cdrfr <- x %>% str_subset("^[CDRF123]+-IMGT+")
    cdrfr_name <- cdrfr %>%
      str_split("\t") %>%
      map_chr(1) %>%
      str_replace("-IMGT", "") %>%
      str_replace(" \\(germline\\)", "")
    if (any(duplicated(cdrfr_name)) == TRUE) {
      stop(str_c("duplicated region names  ", contig_id))
    }
    cdrfr_mut <- cdrfr %>%
      str_split("\t") %>%
      map_chr(6)


    res <- tibble( # nolint
      sequence_id = contig_id,
      name = c(vdj_name, cdrfr_name),
      mut = c(vdj_mut, cdrfr_mut)
    ) %>%
      # fill NA as 0
      dplyr::mutate(mut = as.double(mut), mut = ifelse(is.na(mut), 0, mut)) %>%
      pivot_wider(names_from = "name", values_from = "mut")
  }

  res <- Ltxt %>% map_dfr(parse_mut)
  return(res)
}

#' parse sequences from CellRanger vdj output
#'
#' @param x tibble from airr_rearrangement.tsv
#' @param file output file, default is NULL and will return a tibble
#' @param fa_content if the extension of file is 'fa|fasta', select a column to
#' write in the file. Can be one of `sequence, seq_orf_nt, seq_orf_aa`
#'
#' @return tibble or NULL
#' @export
#'
parse_CellRanger_vdjseq <- function(x, file = NULL, fa_content = "sequence") {
  if (any(c(
    "cell_id", "sequence_id", "sequence",
    "v_sequence_start", "sequence_aa"
  ) %nin% colnames(x))) {
    stop("please read from airr_rearrangement.tsv!")
  }

  x <- x %>%
    dplyr::mutate(
      seq_orf_nt = str_sub(
        .data[["sequence"]],
        .data[["v_sequence_start"]]
      )
    ) %>%
    dplyr::mutate(
      seq_orf_aa =
        nt2aa(.data[["seq_orf_nt"]])
    )

  if (any(x[["seq_orf_aa"]] != x[["sequence_aa"]])) {
    stop("please check the sequences!")
  }

  if (any(!str_detect(x[["seq_orf_aa"]], "^M"))) {
    stop("please check the sequences!")
  }

  res <- x[, c(
    "cell_id", "sequence_id",
    "sequence", "seq_orf_nt", "seq_orf_aa"
  )]

  if (!is.null(file)) {
    out_format <- xfun::file_ext(file)
    if (out_format == "csv") {
      readr::write_excel_csv(res, file)
    } else if (out_format == "tsv") {
      readr::write_tsv(res, file)
    } else if (out_format %in% c("fa", "fasta")) {
      vector <- res %>% dplyr::pull(fa_content, "sequence_id")
      write_fasta(vector, file)
    } else {
      stop("file extension should be one of fa, csv and tsv")
    }
  } else {
    return(res)
  }
}
