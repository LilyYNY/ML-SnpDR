# Machine-learning scores ----------------------------------------------------

.ml_probability_columns <- paste0("prob_C", 1:4)

.ml_score_type <- function(values) {
  values <- as.character(values)
  result <- ifelse(
    grepl("nested|oof", values, ignore.case = TRUE),
    "nested_oof_probability",
    "fitted_model_probability"
  )
  result[is.na(values) | !nzchar(trimws(values))] <- "fitted_model_probability"
  result
}

.ml_model_version <- function(raw, default) {
  lowered <- tolower(names(raw))
  candidate <- names(raw)[match(c("model_version", "selected_model", "model"), lowered, nomatch = 0L)]
  candidate <- candidate[nzchar(candidate)]
  if (length(candidate)) as.character(raw[[candidate[[1L]]]]) else rep(default, nrow(raw))
}

.ml_validate_features <- function(path) {
  features <- .annotation_read_tabular(path)
  required <- c(
    "module_uid", "network", "method", "subtype", "module",
    "feature_missing_count", "feature_qc_pass", "feature_qc_reason"
  )
  missing <- setdiff(required, names(features))
  if (length(missing)) {
    stop("module_features.tsv is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(features$module_uid)) stop("module_features.tsv contains duplicate module_uid values.", call. = FALSE)
  qc <- as.logical(features$feature_qc_pass)
  if (anyNA(qc) || any(!qc)) stop("All modules must pass feature QC before ML scoring.", call. = FALSE)
  identity <- parse_module_uid(features$module_uid)
  if (any(identity$subtype != features$subtype)) {
    stop("module_features.tsv subtype does not match module_uid.", call. = FALSE)
  }
  features
}

#' Prepare all-module ML probabilities and subtype-specific Top-K modules
#'
#' Step 6B validates an all-module score table against `module_features.tsv`,
#' calculates each module's probability for its true subtype, probability
#' margin and within-subtype rank, and writes a controlled Top-K subset. This
#' adapter accepts the nested OOF output of the paper model as well as custom
#' probability tables.
#'
#' @param module_features_file Step-6A `module_features.tsv`.
#' @param score_file Tabular all-module probabilities. It must contain
#'   `prob_C1` through `prob_C4` plus a module identity accepted by
#'   [prepare_module_features()].
#' @param output_dir New step-6B output directory.
#' @param top_k Number of eligible modules retained per true subtype.
#' @param require_predicted_subtype_match Require the predicted class to equal
#'   the module's true subtype before Top-K selection.
#' @param probability_min Optional minimum true-subtype probability.
#' @param margin_min Optional minimum true-subtype probability margin over the
#'   largest competing-class probability.
#' @param model_version Default model version when absent from `score_file`.
#' @param strict Require exact all-module coverage, probability rows summing to
#'   one and at least `top_k` eligible modules per subtype.
#' @return The canonical `ml_scores.tsv` with output paths in attributes.
#' @export
prepare_ml_scores <- function(
  module_features_file,
  score_file,
  output_dir,
  top_k = 10L,
  require_predicted_subtype_match = TRUE,
  probability_min = NULL,
  margin_min = NULL,
  model_version = "custom_model_v1",
  strict = TRUE
) {
  if (length(score_file) != 1L || !file.exists(score_file)) {
    stop("score_file must be one existing tabular file.", call. = FALSE)
  }
  if (length(top_k) != 1L || !is.numeric(top_k) || is.na(top_k) ||
      top_k < 1 || top_k != floor(top_k)) {
    stop("top_k must be a positive integer.", call. = FALSE)
  }
  if (!is.null(probability_min) &&
      (length(probability_min) != 1L || !is.numeric(probability_min) ||
       is.na(probability_min) || probability_min < 0 || probability_min > 1)) {
    stop("probability_min must be NULL or one value between 0 and 1.", call. = FALSE)
  }
  if (!is.null(margin_min) &&
      (length(margin_min) != 1L || !is.numeric(margin_min) ||
       is.na(margin_min) || margin_min < -1 || margin_min > 1)) {
    stop("margin_min must be NULL or one value between -1 and 1.", call. = FALSE)
  }
  if (!is.logical(require_predicted_subtype_match) ||
      length(require_predicted_subtype_match) != 1L ||
      is.na(require_predicted_subtype_match) ||
      !is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("require_predicted_subtype_match and strict must be TRUE or FALSE.", call. = FALSE)
  }

  features <- .ml_validate_features(module_features_file)
  raw <- .annotation_read_tabular(score_file)
  if (!nrow(raw)) stop("ML score input contains no modules.", call. = FALSE)
  identity_pack <- .feature_identity(raw, score_file)
  identity <- identity_pack$identity
  if (anyDuplicated(identity$module_uid)) stop("ML score input contains duplicate module identities.", call. = FALSE)
  missing_probabilities <- setdiff(.ml_probability_columns, names(raw))
  if (length(missing_probabilities)) {
    stop("ML score input is missing: ", paste(missing_probabilities, collapse = ", "), call. = FALSE)
  }
  missing_modules <- setdiff(features$module_uid, identity$module_uid)
  extra_modules <- setdiff(identity$module_uid, features$module_uid)
  if (length(missing_modules) || (isTRUE(strict) && length(extra_modules))) {
    stop(
      "ML score/feature coverage mismatch: missing=", length(missing_modules),
      ", extra=", length(extra_modules),
      call. = FALSE
    )
  }

  match_index <- match(features$module_uid, identity$module_uid)
  raw <- raw[match_index, , drop = FALSE]
  identity <- identity[match_index, , drop = FALSE]
  probabilities <- as.data.frame(
    lapply(raw[, .ml_probability_columns, drop = FALSE], function(x) suppressWarnings(as.numeric(as.character(x)))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(probabilities) <- .ml_probability_columns
  matrix <- as.matrix(probabilities)
  if (any(!is.finite(matrix)) || any(matrix < 0 | matrix > 1)) {
    stop("ML probabilities must be finite values between 0 and 1.", call. = FALSE)
  }
  row_sums <- rowSums(matrix)
  if (isTRUE(strict) && any(abs(row_sums - 1) > 1e-6)) {
    stop("Every strict-mode ML probability row must sum to 1.", call. = FALSE)
  }
  predicted <- paste0("C", max.col(matrix, ties.method = "first"))
  true_subtype <- as.character(features$subtype)
  true_column <- match(paste0("prob_", true_subtype), names(probabilities))
  target_probability <- matrix[cbind(seq_len(nrow(matrix)), true_column)]
  competing <- vapply(seq_len(nrow(matrix)), function(index) {
    max(matrix[index, -true_column[[index]], drop = TRUE])
  }, numeric(1))
  probability_margin <- target_probability - competing

  lowered <- tolower(names(raw))
  score_type_column <- names(raw)[match("score_type", lowered, nomatch = 0L)]
  score_type_column <- score_type_column[nzchar(score_type_column)]
  score_type <- if (length(score_type_column)) {
    .ml_score_type(raw[[score_type_column[[1L]]]])
  } else {
    rep("fitted_model_probability", nrow(raw))
  }
  result <- cbind(
    features[, c("module_uid", "network", "method", "subtype", "module"), drop = FALSE],
    data.frame(
      true_subtype = true_subtype,
      predicted_subtype = predicted,
      probabilities,
      target_subtype_probability = target_probability,
      probability_margin = probability_margin,
      rank_in_subtype = NA_integer_,
      score_type = score_type,
      model_version = .ml_model_version(raw, model_version),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
  for (subtype in sort(unique(result$true_subtype))) {
    indices <- which(result$true_subtype == subtype)
    ordering <- order(
      -result$target_subtype_probability[indices],
      -result$probability_margin[indices],
      result$module_uid[indices]
    )
    result$rank_in_subtype[indices[ordering]] <- seq_along(indices)
  }

  gate <- rep(TRUE, nrow(result))
  reasons <- rep("eligible", nrow(result))
  if (isTRUE(require_predicted_subtype_match)) {
    failed <- result$predicted_subtype != result$true_subtype
    gate[failed] <- FALSE
    reasons[failed] <- "predicted_subtype_mismatch"
  }
  if (!is.null(probability_min)) {
    failed <- result$target_subtype_probability < probability_min
    gate[failed] <- FALSE
    reasons[failed] <- ifelse(reasons[failed] == "eligible", "probability_below_min", paste(reasons[failed], "probability_below_min", sep = ";"))
  }
  if (!is.null(margin_min)) {
    failed <- result$probability_margin < margin_min
    gate[failed] <- FALSE
    reasons[failed] <- ifelse(reasons[failed] == "eligible", "margin_below_min", paste(reasons[failed], "margin_below_min", sep = ";"))
  }

  top_rows <- list()
  for (subtype in sort(unique(result$true_subtype))) {
    eligible <- which(result$true_subtype == subtype & gate)
    eligible <- eligible[order(
      result$rank_in_subtype[eligible],
      result$module_uid[eligible]
    )]
    if (isTRUE(strict) && length(eligible) < top_k) {
      stop(
        "Subtype ", subtype, " has only ", length(eligible),
        " eligible modules; top_k=", top_k, ".",
        call. = FALSE
      )
    }
    selected <- utils::head(eligible, as.integer(top_k))
    if (length(selected)) {
      top_rows[[subtype]] <- cbind(
        result[selected, , drop = FALSE],
        data.frame(
          target_subtype = subtype,
          top_k = as.integer(top_k),
          ml_gate = TRUE,
          ml_gate_reason = "passed_configured_ml_gates",
          stringsAsFactors = FALSE
        )
      )
    }
  }
  top <- if (length(top_rows)) do.call(rbind, top_rows) else result[FALSE, , drop = FALSE]
  rownames(top) <- NULL
  audit <- data.frame(
    module_uid = result$module_uid,
    ml_gate = gate,
    ml_gate_reason = reasons,
    selected_top_k = result$module_uid %in% top$module_uid,
    stringsAsFactors = FALSE
  )

  output_target <- .annotation_output_target(output_dir)
  staging_dir <- tempfile(
    pattern = paste0(".", basename(output_target), "-"),
    tmpdir = dirname(output_target)
  )
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create ML scoring staging directory.", call. = FALSE)
  }
  on.exit({
    if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  scores_file <- file.path(staging_dir, "ml_scores.tsv")
  top_file <- file.path(staging_dir, paste0("ml_top", as.integer(top_k), ".tsv"))
  audit_file <- file.path(staging_dir, "ml_gate_audit.tsv")
  .manifest_write_tsv(result, scores_file)
  .manifest_write_tsv(top, top_file)
  .manifest_write_tsv(audit, audit_file)
  jsonlite::write_json(
    list(
      schema_version = "1.0.0",
      source_score_file = normalizePath(score_file, winslash = "/", mustWork = TRUE),
      source_score_sha256 = digest::digest(score_file, algo = "sha256", file = TRUE),
      module_features_file = normalizePath(module_features_file, winslash = "/", mustWork = TRUE),
      top_k = as.integer(top_k),
      require_predicted_subtype_match = require_predicted_subtype_match,
      probability_min = probability_min,
      margin_min = margin_min,
      selected_module_number = nrow(top)
    ),
    file.path(staging_dir, "ml_scoring_parameters.json"),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  if (!file.rename(staging_dir, output_target)) {
    stop("Could not finalize ML scoring output directory: ", output_target, call. = FALSE)
  }
  attr(result, "output_dir") <- output_target
  attr(result, "scores_file") <- file.path(output_target, basename(scores_file))
  attr(result, "top_file") <- file.path(output_target, basename(top_file))
  attr(result, "audit_file") <- file.path(output_target, basename(audit_file))
  result
}

#' Run the bundled Python nested-OOF machine-learning scorer
#'
#' @param module_features_file Step-6A feature matrix.
#' @param feature_schema_file Step-6A feature schema.
#' @param output_dir New step-6B output directory.
#' @param python_executable Python executable with the `pyproject.toml`
#'   dependencies installed.
#' @param mode `fast` uses the paper-selected GradientBoosting parameters;
#'   `paper` repeats the full inner parameter grid.
#' @param seed,outer_splits,outer_repeats,inner_splits,top_k,n_jobs Nested-CV
#'   settings.
#' @param require_predicted_subtype_match Apply the ML eligibility gate.
#' @return The generated `ml_scores.tsv` table with output paths in attributes.
#' @export
run_nested_ml_scoring <- function(
  module_features_file,
  feature_schema_file,
  output_dir,
  python_executable = "python",
  mode = c("fast", "paper"),
  seed = 20260513L,
  outer_splits = 5L,
  outer_repeats = 5L,
  inner_splits = 3L,
  top_k = 10L,
  n_jobs = 1L,
  require_predicted_subtype_match = TRUE
) {
  mode <- match.arg(mode)
  script <- system.file("python", "mlsnpdr", "scoring.py", package = "MLSnpDR")
  if (!nzchar(script) || !file.exists(script)) {
    stop("Could not locate the bundled Python scoring script.", call. = FALSE)
  }
  for (path in c(module_features_file, feature_schema_file)) {
    if (!file.exists(path)) stop("ML scoring input does not exist: ", path, call. = FALSE)
  }
  arguments <- c(
    shQuote(script),
    "--features", shQuote(normalizePath(module_features_file, winslash = "/", mustWork = TRUE)),
    "--schema", shQuote(normalizePath(feature_schema_file, winslash = "/", mustWork = TRUE)),
    "--output", shQuote(normalizePath(output_dir, winslash = "/", mustWork = FALSE)),
    "--mode", mode,
    "--seed", as.integer(seed),
    "--outer-splits", as.integer(outer_splits),
    "--outer-repeats", as.integer(outer_repeats),
    "--inner-splits", as.integer(inner_splits),
    "--top-k", as.integer(top_k),
    "--n-jobs", as.integer(n_jobs)
  )
  if (!isTRUE(require_predicted_subtype_match)) {
    arguments <- c(arguments, "--allow-predicted-subtype-mismatch")
  }
  command_output <- system2(
    command = python_executable,
    args = arguments,
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(command_output, "status") %||% 0L
  if (status != 0L) {
    stop(
      "Python nested-OOF scoring failed (status ", status, "):\n",
      paste(command_output, collapse = "\n"),
      call. = FALSE
    )
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  scores_file <- file.path(output_dir, "ml_scores.tsv")
  top_file <- file.path(output_dir, paste0("ml_top", as.integer(top_k), ".tsv"))
  scores <- .annotation_read_tabular(scores_file)
  attr(scores, "output_dir") <- output_dir
  attr(scores, "scores_file") <- scores_file
  attr(scores, "top_file") <- top_file
  attr(scores, "metadata_file") <- file.path(output_dir, "ml_scoring_metadata.json")
  attr(scores, "model_file") <- file.path(output_dir, "fitted_model.joblib")
  scores
}
