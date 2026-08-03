# Configuration ---------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

.deep_merge_lists <- function(base, override) {
  result <- base
  for (name in names(override)) {
    if (is.list(result[[name]]) && is.list(override[[name]])) {
      result[[name]] <- .deep_merge_lists(result[[name]], override[[name]])
    } else {
      result[[name]] <- override[[name]]
    }
  }
  result
}

#' Validate an ML-SnpDR configuration
#'
#' @param config Parsed configuration list.
#' @return The validated and normalized configuration.
#' @export
validate_mlsnpdr_config <- function(config) {
  required <- c(
    "project", "module_prefilter", "ml", "survival",
    "candidate_selection", "drug_response", "post_selection"
  )
  missing <- setdiff(required, names(config))
  if (length(missing)) {
    stop("Missing configuration section(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  mode <- config$project$mode %||% "paper_reproduction"
  if (!mode %in% c("paper_reproduction", "fast", "sensitivity")) {
    stop("project.mode must be paper_reproduction, fast or sensitivity.", call. = FALSE)
  }
  config$project$mode <- mode
  config$project$target_subtypes <- .normalize_subtype(config$project$target_subtypes)

  cutoff <- config$module_prefilter$min_size_exclusive
  if (length(cutoff) != 1L || !is.numeric(cutoff) || is.na(cutoff) || cutoff < 0) {
    stop("module_prefilter.min_size_exclusive must be one non-negative number.", call. = FALSE)
  }

  top_k <- config$ml$top_k_per_subtype
  if (length(top_k) != 1L || !is.numeric(top_k) || is.na(top_k) || top_k < 1) {
    stop("ml.top_k_per_subtype must be a positive integer.", call. = FALSE)
  }
  config$ml$top_k_per_subtype <- as.integer(top_k)

  if (!isTRUE(config$post_selection$selected_only)) {
    stop("post_selection.selected_only must be true to protect the downstream boundary.", call. = FALSE)
  }

  config
}

#' Read an ML-SnpDR YAML configuration
#'
#' @param path YAML file to read.
#' @param defaults Optional defaults YAML file merged before `path`.
#' @return A validated configuration list.
#' @export
read_mlsnpdr_config <- function(path, defaults = NULL) {
  if (!file.exists(path)) stop("Configuration file does not exist: ", path, call. = FALSE)
  config <- yaml::read_yaml(path)
  if (!is.null(defaults)) {
    if (!file.exists(defaults)) stop("Defaults file does not exist: ", defaults, call. = FALSE)
    config <- .deep_merge_lists(yaml::read_yaml(defaults), config)
  }
  validate_mlsnpdr_config(config)
}

