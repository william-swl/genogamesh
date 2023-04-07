#' a S4 class to operate mutation strings
#'
#' @slot names character. name of each item
#' @slot mstr character. raw mutation strings
#' @slot sep character. separation between two mutations
#' @slot mut list. splited mutation list
#'
#' @export
#'
setClass("mutstr", representation(
  names = "character",
  mstr = "character",
  sep = "character",
  mut = "list"
))

#' create a mutstr object
#'
#' @param mstr character vector of mutations
#' @param names item names
#' @param sep separation of each mutation
#'
#' @return mutstr
#' @export
#'
#' @examples
#'
#' raw_mut_string <- c(
#'   variant1 = "T10I,D20N,Q30E,A40T,P50L,G60R",
#'   variant2 = "T10I,D20-,Q30E,A40T,P50L,G60R,S80R",
#'   variant3 = "T10A,D20G,Q30E,A40T,P50L,G60R"
#' )
#'
#' m <- mutstr(raw_mut_string, sep = ",")
#'
#' m
#'
#' names(m)
#'
#' mstr(m)
#'
#' mut(m)
#'
#' m[1:2]
#'
#' m[[2]]
#'
#' intersect(m, m[1])
#'
#' setdiff(m, m[1])
#'
#' union(m, m[1])
#'
mutstr <- function(mstr = "", names = NULL, sep = ",") {
  if (!is.null(names)) {
    names(mstr) <- names
  } else if (is.null(names(mstr))) {
    names(mstr) <- stringr::str_c("mut-", seq_along(mstr))
  }
  mut <- stringr::str_split(mstr, sep)
  names(mut) <- names(mstr)
  new("mutstr", names = names(mstr), mstr = mstr, sep = sep, mut = mut)
}

# validation
setValidity("mutstr", function(object) {
  if (
    (length(object@mstr) != length(object@sep)) &&
      (length(object@sep) != 1)) {
    "the length of @sep must be 1 or same as @mstr"
  } else if (
    (length(object@mstr) != length(object@names)) &&
      (!is.null(object@names))) {
    "the length of @names must be 0 or same as @mstr"
  } else if (!is(object@mut, "list")) {
    "@mut must be a list"
  } else {
    TRUE
  }
})

# set @names
setMethod("names<-",
  signature = signature(x = "mutstr", value = "character"),
  function(x, value) {
    x@names <- value
    names(x@mstr) <- value
    names(x@mut) <- value
    validObject(x)
    x
  }
)


# get @mstr
setGeneric("mstr", function(x) standardGeneric("mstr"))
#' @export
setMethod("mstr", "mutstr", function(x) x@mstr)

# set @mstr
setGeneric("mstr<-", function(x, value) standardGeneric("mstr<-"))
#' @export
setMethod("mstr<-",
  signature = signature(x = "mutstr", value = "character"),
  function(x, value) {
    x@mstr <- value
    x@mut <- strsplit(x@mstr, x@sep)
    if (is.null(names(value))) {
      # if the new value have no names
      # use the raw names
      names(x) <- x@names[seq_along(value)]
    } else {
      # if the new value have names
      names(x) <- names(value)
    }
    validObject(x)
    x
  }
)

# get @sep
setGeneric("sep", function(x) standardGeneric("sep"))
#' @export
setMethod("sep", signature = signature(x = "mutstr"), function(x) x@sep)

# set @sep
setGeneric("sep<-", function(x, value) standardGeneric("sep<-"))
#' @export
setMethod("sep<-",
  signature = signature(x = "mutstr", value = "character"),
  function(x, value) {
    x@sep <- value
    x@mut <- strsplit(x@mstr, x@sep)
    validObject(x)
    x
  }
)

# get @mut
setGeneric("mut", function(x) standardGeneric("mut"))
#' @export
setMethod("mut", "mutstr", function(x) x@mut)

# set @mut
setGeneric("mut<-", function(x, value) standardGeneric("mut<-"))
#' @export
setMethod("mut<-",
  signature = signature(x = "mutstr", value = "list"),
  function(x, value) {
    x@mut <- value
    x@mstr <- purrr::map2_chr(x@mut, x@sep, ~ stringr::str_c(.x, collapse = .y))
    if (is.null(names(value))) {
      # if the new value have no names
      # use the raw names
      names(x) <- x@names[seq_along(value)]
    } else {
      # if the new value have names
      names(x) <- names(value)
    }
    validObject(x)
    x
  }
)

# set print format
setMethod("show", "mutstr", function(object) {
  last <- length(object@mstr)
  cat(is(object), last, sep = " ")
  # names
  cat("\n  @ names: ")
  cat(object@names, sep = " ")
  # sep
  cat("\n  @ sep: ", object@sep, sep = "")
  # mstr
  cat("\n  @ mstr:")
  cat("\n    [1]", object@mstr[[1]], sep = " ")
  if (last > 3) {
    cat("\n    ...")
  } else if (last %in% 2:3) {
    cat("\n    [2]", object@mstr[[2]], sep = " ")
  }
  if (last >= 3) {
    cat("\n    [", last, "] ", sep = "")
    cat(object@mstr[[last]], sep = " ")
  }
  # mut
  cat("\n  @ mut:")
  cat("\n    [1]", object@mut[[1]], sep = " ")
  if (last > 3) {
    cat("\n    ...")
  } else if (last %in% 2:3) {
    cat("\n    [2]", object@mut[[2]], sep = " ")
  }
  if (last >= 3) {
    cat("\n    [", last, "] ", sep = "")
    cat(object@mut[[last]], sep = " ")
  }

  cat("\n")
})


# subsettable
#' @export
setMethod("[",
  signature = signature(x = "mutstr", i = "numeric", j = "missing"),
  function(x, i, j = "missing") mutstr(mstr(x)[i])
)

#' @export
setMethod("[<-",
  signature = signature(
    x = "mutstr", i = "numeric",
    j = "missing", value = "list"
  ),
  function(x, i, j = "missing", value) {
    mut(x)[i] <- value
    x
  }
)

#' @export
setMethod("[[",
  signature =
    signature(x = "mutstr", i = "numeric", j = "missing"),
  function(x, i, j = "missing") mut(x)[[i]]
)

#' @export
setMethod("[[<-",
  signature = signature(
    x = "mutstr", i = "numeric",
    j = "missing", value = "character"
  ),
  function(x, i, j = "missing", value) {
    mut(x)[[i]] <- value
    x
  }
)

# intersect
setGeneric("intersect", function(x, y) standardGeneric("intersect"))
#' @export
setMethod("intersect",
  signature = signature(x = "mutstr", y = "mutstr"),
  function(x, y) {
    mut <- purrr::map2(mut(x), mut(y), ~ intersect(.x, .y))
    res <- mutstr()
    mut(res) <- mut
    names(res) <- names(mut)
    return(res)
  }
)

# setdiff
setGeneric("setdiff", function(x, y) standardGeneric("setdiff"))
#' @export
setMethod("setdiff",
  signature = signature(x = "mutstr", y = "mutstr"),
  function(x, y) {
    mut <- purrr::map2(mut(x), mut(y), ~ setdiff(.x, .y))
    res <- mutstr()
    mut(res) <- mut
    names(res) <- names(mut)
    return(res)
  }
)

# union
setGeneric("union", function(x, y) standardGeneric("union"))
#' @export
setMethod("union",
  signature = signature(x = "mutstr", y = "mutstr"),
  function(x, y) {
    mut <- purrr::map2(mut(x), mut(y), ~ union(.x, .y))
    res <- mutstr()
    mut(res) <- mut
    names(res) <- names(mut)
    return(res)
  }
)
