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
    "project", "paths", "pipeline", "differential_expression", "network_construction",
    "module_division", "module_prefilter", "annotation", "ml", "survival",
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
  config$pipeline$network_methods <- normalize_network_name(
    config$pipeline$network_methods,
    output = "display"
  )
  config$pipeline$module_methods <- normalize_method_name(
    config$pipeline$module_methods,
    output = "display"
  )
  if (!isTRUE(config$pipeline$strict_inputs)) {
    stop("pipeline.strict_inputs must be true for reproducible pipeline execution.", call. = FALSE)
  }
  if (!is.logical(config$pipeline$write_module_plots) ||
      length(config$pipeline$write_module_plots) != 1L ||
      is.na(config$pipeline$write_module_plots)) {
    stop("pipeline.write_module_plots must be true or false.", call. = FALSE)
  }

  dep_probabilities <- c(
    detection_threshold = config$differential_expression$detection_threshold,
    p_adjust_threshold = config$differential_expression$p_adjust_threshold
  )
  if (any(!is.finite(dep_probabilities)) || any(dep_probabilities < 0) ||
      any(dep_probabilities > 1)) {
    stop("Differential-expression probability thresholds must be between 0 and 1.", call. = FALSE)
  }
  dep_fc <- config$differential_expression$fc_threshold
  if (length(dep_fc) != 1L || !is.finite(dep_fc) || dep_fc <= 1) {
    stop("differential_expression.fc_threshold must be greater than 1.", call. = FALSE)
  }
  if (!config$differential_expression$p_adjust_method %in% stats::p.adjust.methods) {
    stop("differential_expression.p_adjust_method is unsupported.", call. = FALSE)
  }
  dep_pseudocount <- config$differential_expression$pseudocount
  if (length(dep_pseudocount) != 1L || !is.finite(dep_pseudocount) || dep_pseudocount < 0) {
    stop("differential_expression.pseudocount must be non-negative.", call. = FALSE)
  }
  if (!is.logical(config$differential_expression$write_legacy_excel) ||
      length(config$differential_expression$write_legacy_excel) != 1L ||
      is.na(config$differential_expression$write_legacy_excel)) {
    stop("differential_expression.write_legacy_excel must be true or false.", call. = FALSE)
  }
  labels <- tolower(gsub("-", "_", trimws(as.character(config$network_construction$include_labels))))
  if (!length(labels) || any(!nzchar(labels))) {
    stop("network_construction.include_labels must contain at least one label.", call. = FALSE)
  }
  config$network_construction$include_labels <- unique(labels)
  ppi_score <- config$network_construction$ppi_score_threshold
  if (length(ppi_score) != 1L || !is.finite(ppi_score)) {
    stop("network_construction.ppi_score_threshold must be finite.", call. = FALSE)
  }
  division_seed <- config$module_division$seed
  if (length(division_seed) != 1L || !is.finite(division_seed) ||
      division_seed != floor(division_seed)) {
    stop("module_division.seed must be an integer.", call. = FALSE)
  }
  config$module_division$seed <- as.integer(division_seed)
  wf_cutoff <- config$module_division$wf_pvalue_cutoff
  if (length(wf_cutoff) != 1L || !is.finite(wf_cutoff) || wf_cutoff <= 0 || wf_cutoff > 1) {
    stop("module_division.wf_pvalue_cutoff must be in (0, 1].", call. = FALSE)
  }

  cutoff <- config$module_prefilter$min_size_exclusive
  if (length(cutoff) != 1L || !is.numeric(cutoff) || is.na(cutoff) || cutoff < 0) {
    stop("module_prefilter.min_size_exclusive must be one non-negative number.", call. = FALSE)
  }

  annotation_probabilities <- c(
    pvalue_cutoff = config$annotation$pvalue_cutoff,
    p_adjust_cutoff = config$annotation$p_adjust_cutoff,
    qvalue_cutoff = config$annotation$qvalue_cutoff
  )
  if (any(!is.finite(annotation_probabilities)) ||
      any(annotation_probabilities < 0) || any(annotation_probabilities > 1)) {
    stop("annotation probability cutoffs must be between 0 and 1.", call. = FALSE)
  }
  if (!config$annotation$p_adjust_method %in% stats::p.adjust.methods) {
    stop("annotation.p_adjust_method is not supported by stats::p.adjust().", call. = FALSE)
  }
  annotation_integers <- c(
    min_gene_set_size = config$annotation$min_gene_set_size,
    max_gene_set_size = config$annotation$max_gene_set_size,
    top_n = config$annotation$top_n
  )
  if (any(!is.finite(annotation_integers)) || any(annotation_integers < 1) ||
      any(annotation_integers != floor(annotation_integers)) ||
      annotation_integers[["min_gene_set_size"]] > annotation_integers[["max_gene_set_size"]]) {
    stop("annotation gene-set sizes and top_n must be valid positive integers.", call. = FALSE)
  }
  config$annotation$min_gene_set_size <- as.integer(config$annotation$min_gene_set_size)
  config$annotation$max_gene_set_size <- as.integer(config$annotation$max_gene_set_size)
  config$annotation$top_n <- as.integer(config$annotation$top_n)
  config$annotation$databases <- toupper(trimws(as.character(config$annotation$databases)))
  if (!length(config$annotation$databases) || any(!nzchar(config$annotation$databases))) {
    stop("annotation.databases must contain at least one database name.", call. = FALSE)
  }
  if (!isTRUE(config$annotation$write_module_files)) {
    stop("annotation.write_module_files must be true for auditable output.", call. = FALSE)
  }

  top_k <- config$ml$top_k_per_subtype
  if (length(top_k) != 1L || !is.numeric(top_k) || is.na(top_k) || top_k < 1) {
    stop("ml.top_k_per_subtype must be a positive integer.", call. = FALSE)
  }
  config$ml$top_k_per_subtype <- as.integer(top_k)
  expected_features <- config$ml$expected_feature_count
  if (length(expected_features) != 1L || !is.numeric(expected_features) ||
      is.na(expected_features) || expected_features < 1 ||
      expected_features != floor(expected_features)) {
    stop("ml.expected_feature_count must be a positive integer.", call. = FALSE)
  }
  config$ml$expected_feature_count <- as.integer(expected_features)
  cv_integers <- c(
    outer_splits = config$ml$outer_splits,
    outer_repeats = config$ml$outer_repeats,
    inner_splits = config$ml$inner_splits,
    n_jobs = config$ml$n_jobs
  )
  if (any(!is.finite(cv_integers)) || any(cv_integers < 1) ||
      any(cv_integers != floor(cv_integers))) {
    stop("ML cross-validation settings must be positive integers.", call. = FALSE)
  }
  config$ml$outer_splits <- as.integer(config$ml$outer_splits)
  config$ml$outer_repeats <- as.integer(config$ml$outer_repeats)
  config$ml$inner_splits <- as.integer(config$ml$inner_splits)
  config$ml$n_jobs <- as.integer(config$ml$n_jobs)

  if (!identical(config$drug_response$analysis_scope, "all_manifest_modules")) {
    stop(
      "drug_response.analysis_scope must be all_manifest_modules so step 6 covers every step-4 module.",
      call. = FALSE
    )
  }
  config$drug_response$panels <- toupper(trimws(as.character(config$drug_response$panels)))
  if (!length(config$drug_response$panels) || any(!nzchar(config$drug_response$panels)) ||
      anyDuplicated(config$drug_response$panels)) {
    stop("drug_response.panels must contain unique non-empty panel names.", call. = FALSE)
  }
  drug_cutoff <- config$drug_response$significance_cutoff
  if (length(drug_cutoff) != 1L || !is.numeric(drug_cutoff) || is.na(drug_cutoff) ||
      drug_cutoff < 0 || drug_cutoff > 1) {
    stop("drug_response.significance_cutoff must be between 0 and 1.", call. = FALSE)
  }
  if (!config$drug_response$p_adjust_method %in% c("none", stats::p.adjust.methods)) {
    stop("drug_response.p_adjust_method is unsupported.", call. = FALSE)
  }

  select_n <- config$candidate_selection$select_n_per_subtype
  if (length(select_n) != 1L || !is.numeric(select_n) || is.na(select_n) || select_n != 1) {
    stop("candidate_selection.select_n_per_subtype must be 1.", call. = FALSE)
  }
  config$candidate_selection$select_n_per_subtype <- 1L

  if (!isTRUE(config$post_selection$selected_only)) {
    stop("post_selection.selected_only must be true to protect the downstream boundary.", call. = FALSE)
  }
  if (!config$post_selection$binding_score_direction %in% c("lower_better", "higher_better")) {
    stop("post_selection.binding_score_direction must be lower_better or higher_better.", call. = FALSE)
  }
  final_top_n <- config$post_selection$final_top_n_per_subtype
  if (length(final_top_n) != 1L || !is.numeric(final_top_n) || is.na(final_top_n) ||
      final_top_n < 1 || final_top_n != floor(final_top_n)) {
    stop("post_selection.final_top_n_per_subtype must be a positive integer.", call. = FALSE)
  }
  config$post_selection$final_top_n_per_subtype <- as.integer(final_top_n)

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
