# Canonical module identifiers -------------------------------------------------

.network_aliases <- c(
  string = "String",
  physicalppin = "physicalPPIN",
  chengf = "chengF"
)

.method_aliases <- c(
  louvain = "Louvain",
  wf = "WF"
)

.normalize_token <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", trimws(as.character(x))))
}

#' Normalize a PPI network name
#'
#' @param x Character vector of network names.
#' @param output Return controlled display names or lowercase slugs.
#' @return A character vector.
#' @export
normalize_network_name <- function(x, output = c("display", "slug")) {
  output <- match.arg(output)
  key <- .normalize_token(x)
  unknown <- is.na(match(key, names(.network_aliases)))
  if (any(unknown)) {
    stop("Unknown network name(s): ", paste(unique(x[unknown]), collapse = ", "), call. = FALSE)
  }
  if (identical(output, "slug")) key else unname(.network_aliases[key])
}

#' Normalize a module-detection method name
#'
#' @param x Character vector of method names.
#' @param output Return controlled display names or lowercase slugs.
#' @return A character vector.
#' @export
normalize_method_name <- function(x, output = c("display", "slug")) {
  output <- match.arg(output)
  key <- .normalize_token(x)
  unknown <- is.na(match(key, names(.method_aliases)))
  if (any(unknown)) {
    stop("Unknown module method name(s): ", paste(unique(x[unknown]), collapse = ", "), call. = FALSE)
  }
  if (identical(output, "slug")) key else unname(.method_aliases[key])
}

.normalize_subtype <- function(x) {
  value <- toupper(trimws(as.character(x)))
  if (any(!grepl("^C[1-4]$", value))) {
    stop("Subtype must be one of C1, C2, C3 or C4.", call. = FALSE)
  }
  value
}

.normalize_module <- function(x) {
  value <- toupper(trimws(as.character(x)))
  if (any(!grepl("^M[0-9]+$", value))) {
    stop("Module must use the form M<number>.", call. = FALSE)
  }
  value
}

#' Build the canonical module UID
#'
#' @param network PPI network name.
#' @param method Module-detection method.
#' @param subtype Molecular subtype.
#' @param module Module label such as M10.
#' @return A character vector such as `physicalppin__louvain__C3__M10`.
#' @export
make_module_uid <- function(network, method, subtype, module) {
  lengths <- c(length(network), length(method), length(subtype), length(module))
  target_length <- max(lengths)
  if (any(target_length %% lengths != 0L)) {
    stop("Identifier fields must have compatible vector lengths.", call. = FALSE)
  }
  paste(
    normalize_network_name(network, "slug"),
    normalize_method_name(method, "slug"),
    .normalize_subtype(subtype),
    .normalize_module(module),
    sep = "__"
  )
}

#' Parse canonical module UIDs
#'
#' @param module_uid Character vector returned by `make_module_uid()`.
#' @return A data frame with canonical identity fields and a legacy identifier.
#' @export
parse_module_uid <- function(module_uid) {
  parsed <- lapply(as.character(module_uid), function(uid) {
    fields <- strsplit(uid, "__", fixed = TRUE)[[1]]
    if (length(fields) != 4L) {
      stop("Invalid module_uid: ", uid, call. = FALSE)
    }
    network <- normalize_network_name(fields[1], "display")
    method <- normalize_method_name(fields[2], "display")
    subtype <- .normalize_subtype(fields[3])
    module <- .normalize_module(fields[4])
    data.frame(
      module_uid = make_module_uid(network, method, subtype, module),
      legacy_module_id = paste(network, method, subtype, module, sep = "|"),
      network = network,
      method = method,
      subtype = subtype,
      module = module,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parsed)
}
