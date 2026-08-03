# Shared input/output helpers for steps 1-3 ---------------------------------

.upstream_read_tabular <- function(path, sheet = NULL) {
  if (length(path) != 1L || !file.exists(path)) {
    stop("Input table does not exist: ", path, call. = FALSE)
  }
  extension <- tolower(tools::file_ext(path))
  switch(
    extension,
    csv = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    tsv = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    txt = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    xlsx = {
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        stop("Reading .xlsx files requires the openxlsx package.", call. = FALSE)
      }
      arguments <- list(xlsxFile = path, check.names = FALSE)
      if (!is.null(sheet)) arguments$sheet <- sheet
      do.call(openxlsx::read.xlsx, arguments)
    },
    stop("Unsupported table extension: .", extension, call. = FALSE)
  )
}

.upstream_output_target <- function(path, label) {
  if (length(path) != 1L || is.na(path) || !nzchar(trimws(path))) {
    stop(label, " output directory must be one non-empty path.", call. = FALSE)
  }
  path <- path.expand(path)
  parent <- dirname(path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output parent directory: ", parent, call. = FALSE)
  }
  target <- file.path(normalizePath(parent, winslash = "/", mustWork = TRUE), basename(path))
  if (file.exists(target)) {
    stop(label, " output directory already exists: ", target, call. = FALSE)
  }
  target
}

.upstream_column <- function(data, requested, aliases, label) {
  lowered <- tolower(names(data))
  choices <- unique(tolower(c(requested, aliases)))
  hit <- match(choices, lowered, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (!length(hit)) {
    stop(label, " column was not found. Accepted names: ", paste(choices, collapse = ", "), call. = FALSE)
  }
  names(data)[hit[[1L]]]
}

.upstream_relative <- function(path) gsub("\\\\", "/", path)

