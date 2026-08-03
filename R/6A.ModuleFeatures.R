# Module features ------------------------------------------------------------

.feature_identity <- function(data, path) {
  lowered <- tolower(names(data))
  uid_column <- names(data)[match("module_uid", lowered, nomatch = 0L)]
  uid_column <- uid_column[nzchar(uid_column)]
  if (length(uid_column)) {
    identity <- parse_module_uid(trimws(as.character(data[[uid_column[[1L]]]])))
    return(list(identity = identity, source_identity_columns = uid_column[[1L]]))
  }

  module_id_column <- names(data)[match("module_id", lowered, nomatch = 0L)]
  module_id_column <- module_id_column[nzchar(module_id_column)]
  if (length(module_id_column)) {
    raw_ids <- trimws(as.character(data[[module_id_column[[1L]]]]))
    fields <- strsplit(raw_ids, "|", fixed = TRUE)
    invalid <- lengths(fields) != 4L
    if (any(invalid)) {
      stop("module_id must use Network|Method|Subtype|Module in: ", path, call. = FALSE)
    }
    matrix <- do.call(rbind, fields)
    identity <- data.frame(
      module_uid = make_module_uid(matrix[, 1L], matrix[, 2L], matrix[, 3L], matrix[, 4L]),
      legacy_module_id = paste(
        normalize_network_name(matrix[, 1L]),
        normalize_method_name(matrix[, 2L]),
        .normalize_subtype(matrix[, 3L]),
        .normalize_module(matrix[, 4L]),
        sep = "|"
      ),
      network = normalize_network_name(matrix[, 1L]),
      method = normalize_method_name(matrix[, 2L]),
      subtype = .normalize_subtype(matrix[, 3L]),
      module = .normalize_module(matrix[, 4L]),
      stringsAsFactors = FALSE
    )
    return(list(identity = identity, source_identity_columns = module_id_column[[1L]]))
  }

  aliases <- list(
    network = c("network"),
    method = c("method", "module_method"),
    subtype = c("subtype"),
    module = c("module", "module_name")
  )
  selected <- vapply(aliases, function(options) {
    hit <- match(options, lowered, nomatch = 0L)
    hit <- hit[hit > 0L]
    if (length(hit)) names(data)[hit[[1L]]] else NA_character_
  }, character(1))
  if (anyNA(selected)) {
    stop(
      "Feature input must contain module_uid, module_id, or Network/Method/Subtype/Module columns: ",
      path,
      call. = FALSE
    )
  }
  network <- normalize_network_name(data[[selected[["network"]]]])
  method <- normalize_method_name(data[[selected[["method"]]]])
  subtype <- .normalize_subtype(data[[selected[["subtype"]]]]
  )
  module <- .normalize_module(data[[selected[["module"]]]]
  )
  list(
    identity = data.frame(
      module_uid = make_module_uid(network, method, subtype, module),
      legacy_module_id = paste(network, method, subtype, module, sep = "|"),
      network = network,
      method = method,
      subtype = subtype,
      module = module,
      stringsAsFactors = FALSE
    ),
    source_identity_columns = unname(selected)
  )
}

.feature_default_columns <- function(data, identity_columns) {
  metadata <- c(
    identity_columns,
    names(data)[tolower(names(data)) %in% c(
      "module_uid", "module_id", "module_label", "legacy_module_id",
      "network", "method", "module_method", "subtype", "module", "module_name"
    )]
  )
  setdiff(names(data), unique(metadata))
}

#' Prepare the all-module feature matrix used by machine learning
#'
#' Step 6A joins an existing Core34 (or custom) feature table to the exact
#' modules in `module_manifest.tsv`. It validates one-to-one coverage, numeric
#' values, feature order and module sizes, then writes the only feature matrix
#' accepted by step 6B.
#'
#' @param module_manifest_file Step-4 `module_manifest.tsv`.
#' @param feature_file `.tsv`, `.csv` or `.xlsx` table containing either
#'   `module_uid`, `module_id` (`Network|Method|Subtype|Module`), or separate
#'   identity columns.
#' @param output_dir New step-6A output directory.
#' @param feature_columns Optional ordered feature column names. By default all
#'   non-identity/non-label columns are used in source order.
#' @param feature_set Name written to the feature schema.
#' @param expected_feature_count Required number of features; use `NULL` for a
#'   custom unrestricted feature set.
#' @param strict Require exact manifest coverage, no extra rows, finite feature
#'   values, verified step-4 hashes and matching `module_size`.
#' @return The canonical `module_features.tsv` table with output paths in
#'   attributes.
#' @export
prepare_module_features <- function(
  module_manifest_file,
  feature_file,
  output_dir,
  feature_columns = NULL,
  feature_set = "core34_v1",
  expected_feature_count = 34L,
  strict = TRUE
) {
  if (length(feature_file) != 1L || !file.exists(feature_file)) {
    stop("feature_file must be one existing tabular file.", call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("strict must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(expected_feature_count) &&
      (length(expected_feature_count) != 1L || !is.numeric(expected_feature_count) ||
       is.na(expected_feature_count) || expected_feature_count < 1 ||
       expected_feature_count != floor(expected_feature_count))) {
    stop("expected_feature_count must be NULL or one positive integer.", call. = FALSE)
  }
  if (length(feature_set) != 1L || is.na(feature_set) || !nzchar(trimws(feature_set))) {
    stop("feature_set must be one non-empty name.", call. = FALSE)
  }

  manifest <- read_module_manifest(
    module_manifest_file,
    check_files = TRUE,
    verify_hashes = strict
  )
  raw <- .annotation_read_tabular(feature_file)
  if (!nrow(raw)) stop("Feature input contains no modules: ", feature_file, call. = FALSE)
  identity_pack <- .feature_identity(raw, feature_file)
  identity <- identity_pack$identity
  if (anyDuplicated(identity$module_uid)) {
    stop("Feature input contains duplicate module identities.", call. = FALSE)
  }

  if (is.null(feature_columns)) {
    feature_columns <- .feature_default_columns(raw, identity_pack$source_identity_columns)
  }
  feature_columns <- as.character(feature_columns)
  if (!length(feature_columns) || any(!nzchar(feature_columns)) || anyDuplicated(feature_columns)) {
    stop("feature_columns must contain unique non-empty column names.", call. = FALSE)
  }
  missing_features <- setdiff(feature_columns, names(raw))
  if (length(missing_features)) {
    stop("Feature input is missing requested column(s): ", paste(missing_features, collapse = ", "), call. = FALSE)
  }
  if (!is.null(expected_feature_count) && length(feature_columns) != expected_feature_count) {
    stop(
      "Expected ", expected_feature_count, " feature columns but found ", length(feature_columns), ".",
      call. = FALSE
    )
  }

  missing_modules <- setdiff(manifest$module_uid, identity$module_uid)
  extra_modules <- setdiff(identity$module_uid, manifest$module_uid)
  if (length(missing_modules) || (isTRUE(strict) && length(extra_modules))) {
    stop(
      "Feature/manifest coverage mismatch: missing=", length(missing_modules),
      ", extra=", length(extra_modules),
      if (length(missing_modules)) paste0("; first missing=", missing_modules[[1L]]) else "",
      call. = FALSE
    )
  }
  match_index <- match(manifest$module_uid, identity$module_uid)
  raw <- raw[match_index, , drop = FALSE]
  identity <- identity[match_index, , drop = FALSE]
  values <- lapply(raw[, feature_columns, drop = FALSE], function(column) {
    suppressWarnings(as.numeric(as.character(column)))
  })
  values <- as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
  names(values) <- feature_columns
  missing_count <- rowSums(!is.finite(as.matrix(values)))
  size_match <- rep(TRUE, nrow(manifest))
  if ("module_size" %in% feature_columns) {
    size_match <- is.finite(values$module_size) & values$module_size == manifest$module_size
  }
  qc_pass <- missing_count == 0L & size_match
  qc_reason <- ifelse(
    qc_pass,
    "pass",
    paste0(
      ifelse(missing_count > 0L, paste0("non_finite_features=", missing_count), ""),
      ifelse(missing_count > 0L & !size_match, ";", ""),
      ifelse(!size_match, "module_size_mismatch", "")
    )
  )
  if (isTRUE(strict) && any(!qc_pass)) {
    failed <- manifest$module_uid[!qc_pass]
    stop(
      "Feature QC failed for ", length(failed), " module(s); first: ", failed[[1L]],
      " [", qc_reason[!qc_pass][[1L]], "]",
      call. = FALSE
    )
  }

  result <- cbind(
    manifest[, c("module_uid", "network", "method", "subtype", "module"), drop = FALSE],
    values,
    data.frame(
      feature_missing_count = as.integer(missing_count),
      feature_qc_pass = qc_pass,
      feature_qc_reason = qc_reason,
      stringsAsFactors = FALSE
    )
  )
  output_target <- .annotation_output_target(output_dir)
  staging_dir <- tempfile(
    pattern = paste0(".", basename(output_target), "-"),
    tmpdir = dirname(output_target)
  )
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create module-feature staging directory.", call. = FALSE)
  }
  on.exit({
    if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  feature_file_out <- file.path(staging_dir, "module_features.tsv")
  schema_file <- file.path(staging_dir, "feature_schema.json")
  qc_file <- file.path(staging_dir, "module_features_qc.tsv")
  .manifest_write_tsv(result, feature_file_out)
  .manifest_write_tsv(
    data.frame(
      module_uid = manifest$module_uid,
      source_row = match_index,
      source_identity = identity$legacy_module_id,
      feature_missing_count = missing_count,
      module_size_match = size_match,
      feature_qc_pass = qc_pass,
      feature_qc_reason = qc_reason,
      stringsAsFactors = FALSE
    ),
    qc_file
  )
  jsonlite::write_json(
    list(
      schema_version = "1.0.0",
      feature_set = feature_set,
      feature_count = length(feature_columns),
      feature_columns = feature_columns,
      identity_columns = c("module_uid", "network", "method", "subtype", "module"),
      qc_columns = c("feature_missing_count", "feature_qc_pass", "feature_qc_reason"),
      source_file = normalizePath(feature_file, winslash = "/", mustWork = TRUE),
      source_sha256 = digest::digest(feature_file, algo = "sha256", file = TRUE),
      manifest_file = normalizePath(module_manifest_file, winslash = "/", mustWork = TRUE),
      manifest_sha256 = digest::digest(module_manifest_file, algo = "sha256", file = TRUE)
    ),
    schema_file,
    auto_unbox = TRUE,
    pretty = TRUE
  )

  if (!file.rename(staging_dir, output_target)) {
    stop("Could not finalize module-feature output directory: ", output_target, call. = FALSE)
  }
  attr(result, "output_dir") <- output_target
  attr(result, "feature_file") <- file.path(output_target, basename(feature_file_out))
  attr(result, "schema_file") <- file.path(output_target, basename(schema_file))
  attr(result, "qc_file") <- file.path(output_target, basename(qc_file))
  result
}
