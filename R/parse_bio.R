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

  # celltype in colnames have some duplication
  # e.g. "scores.NK.cells..NK.H..MCMV1." and
  # "scores.NK.cells..NK.H.MCMV1." in mouse reference
  res <- res %>%
    dplyr::arrange(.data[["cell"]], dplyr::desc(.data[["s"]])) %>%
    dplyr::distinct(.data[["cell"]], .keep_all = TRUE)

  # select columns
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

  # rename
  res <- res %>% rename_at(
    str_subset(colnames(.), "^FR\\d"),
    ~ str_replace(.x, "FR", "FWR")
  )
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



#' parse sequences from ANARCI vdj output
#'
#' @param x tibble from anarci_H.csv or anarci_KL.csv
#' @param remove_gap remove the gap caused by numbering, default as TRUE
#' @param chain one of 'H, L', for heavy chain or light chain
#' @param scheme antibody numbering scheme, one of 'imgt, chothia'
#' @param number_table custom antibody numbering system, please input a tibble
#' with three columns: region, start, end.
#'
#' @return tibble
#' @export
#'
parse_ANARCI_aaseq <- function(x, chain, remove_gap = TRUE,
                               scheme = "imgt", number_table = NULL,
                               keep_number = FALSE) {
  aa_cols <- colnames(x) %>% str_subset("^\\d+")
  res <- x %>%
    dplyr::select(sequence_id = "Id", dplyr::all_of(aa_cols)) %>%
    # trans TRUE, FALSE to the real character T, F
    dplyr::mutate(
      dplyr::across(
        dplyr::where(is.logical),
        ~ case_when(
          .x == TRUE ~ "T", .x == FALSE ~ "F",
          is.na(.x) ~ "-", TRUE ~ as.character(.x)
        )
      )
    ) %>%
    # trans NA to -, to avoid NA in sequence
    dplyr::mutate(
      dplyr::across(dplyr::all_of(aa_cols), ~ ifelse(is.na(.x), "-", .x))
    )

  # region seq aa
  region <- res %>%
    pivot_longer(-"sequence_id", names_to = "numbering", values_to = "aa") %>%
    dplyr::mutate(numbering = factor(.data[["numbering"]], aa_cols)) %>%
    dplyr::mutate(site = str_extract(.data[["numbering"]], "\\d+") %>%
      as.double()) %>%
    dplyr::arrange(.data[["sequence_id"]], .data[["numbering"]])

  if (is.null(number_table)) {
    number_table <- ab_numbering %>%
      dplyr::filter(scheme == .env[["scheme"]], chain == .env[["chain"]])
  } else if ("chain" %in% colnames(number_table)) {
    number_table <- number_table %>% dplyr::filter(chain == .env[["chain"]])
  }

  region <- region %>%
    left_join(number_table,
      by = dplyr::join_by(between(site, start, end, bounds = "[]"))
    ) %>%
    dplyr::filter(!is.na(region)) %>%
    dplyr::arrange(region)

  region <- region %>%
    group_split(.data[["sequence_id"]]) %>%
    map_dfr(~ dplyr::summarise(.x,
      aa = str_c(aa, collapse = ""),
      .by = c("sequence_id", "region")
    )) %>%
    pivot_wider(names_from = "region", values_from = "aa") %>%
    rename_at(
      dplyr::pull(number_table, "region") %>% as.character(),
      ~ str_c(.x, "_aa")
    )

  # V-domain seq aa
  res <- res %>% tidyr::unite("seq_align_aa", dplyr::all_of(aa_cols),
    sep = "", remove = !keep_number
  )

  res <- left_join(region, res, by = "sequence_id")

  if (remove_gap == TRUE) {
    res <- res %>% dplyr::mutate(dplyr::across(
      ends_with("_aa"),
      ~ str_replace_all(.x, "-", "")
    ))
  }

  return(res)
}



#' parse vcf file with the help of reference genome and annotations.
#' It is still under development, which can not process more than 3 nt
#' substitutions in a single record row of vcf file, and can not process indels
#'
#' @param vcf `vcfR` object from `read_vcf`
#' @param fa `DNAStringSet` object from `read_fasta`
#' @param gff `gff` tibble from `read_gff`
#'
#' @return tibble
#' @export
#'
parse_vcf <- function(vcf, fa, gff) {
  # get gff info
  vcf_fmt <- vcf@gt[, 2] %>%
    str_split(":") %>%
    baizer::list2df(
      rownames = FALSE,
      colnames = vcf@gt[1, 1] %>% str_split(":") %>% unlist()
    )

  res <- as_tibble(vcf@fix) %>%
    dplyr::mutate(POS = as.numeric(.data[["POS"]])) %>%
    bind_cols(vcf_fmt) %>%
    left_join(gff,
      by = join_by(
        "CHROM" == "seqid",
        between("POS", "start", "end", bounds = "[]")
      )
    )
  # fetch codon range
  res <- res %>%
    dplyr::mutate(
      codon_range =
        codon_range(
          .data[["POS"]],
          .data[["start"]],
          .data[["end"]]
        )
    ) %>%
    tidyr::separate(codon_range, into = c("codon_start", "codon_end"))

  # fetch codon ref
  res <- res %>% dplyr::mutate(
    codon_ref = as.character(fa[.data[["CHROM"]]]) %>%
      str_sub(.data[["codon_start"]], .data[["codon_end"]])
  )

  # alt codon index
  res <- res %>%
    dplyr::mutate(
      codon_idx_start = codon_index(
        .data[["POS"]], .data[["start"]], .data[["end"]]
      ),
      codon_idx_end =
        codon_index(
          .data[["POS"]], .data[["start"]],
          .data[["end"]]
        ) + nchar(.data[["ALT"]]) - 1
    ) %>%
    dplyr::mutate(
      codon_overflow =
        ifelse(.data[["codon_idx_end"]] > 3, TRUE, FALSE)
    )

  # alt codon
  res <- res %>% dplyr::mutate(
    codon_alt =
      baizer::str_replace_loc(
        .data[["codon_ref"]], .data[["codon_idx_start"]],
        .data[["codon_idx_end"]], .data[["ALT"]]
      )
  )

  # aa index
  res <- res %>% dplyr::mutate(
    aa_index = aa_index(.data[["POS"]], .data[["start"]], .data[["end"]])
  )

  # aa mut
  res <- res %>% dplyr::mutate(
    nonsynonymous = nt2aa(.data[["codon_ref"]]) != nt2aa(.data[["codon_alt"]]),
    aa_mut = str_c(
      nt2aa(.data[["codon_ref"]]),
      .data[["aa_index"]], nt2aa(.data[["codon_alt"]])
    )
  )

  # select and relocate
  res <- res %>%
    dplyr::select(-c(
      "ID", "source", "start", "end", "score",
      "strand", "phase", "codon_start", "codon_end",
      "codon_idx_start", "codon_idx_end", "aa_index"
    )) %>%
    dplyr::relocate("codon_overflow", .after = -1)

  n_overflow <- sum(res$codon_overflow)
  if (n_overflow > 0) {
    cat(str_glue("{n_overflow} records codon overflow,
                 please check by parse_vcf(...) %>% filter() == TRUE !"))
  } else {
    res <- res %>% dplyr::select(-"codon_overflow")
  }

  return(res)
}
