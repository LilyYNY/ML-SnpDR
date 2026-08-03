# Module triage --------------------------------------------------------------

.triage_find_column <- function(data, aliases, label, required = TRUE) {
  lowered <- tolower(names(data))
  hit <- match(tolower(aliases), lowered, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (!length(hit)) {
    if (required) stop("Could not find ", label, " column in survival input.", call. = FALSE)
    return(NA_character_)
  }
  names(data)[hit[[1L]]]
}

.triage_numeric <- function(data, column) {
  if (is.na(column)) return(rep(NA_real_, nrow(data)))
  suppressWarnings(as.numeric(as.character(data[[column]])))
}

.triage_standardize_survival <- function(raw, manifest, path, strict) {
  identity_pack <- .feature_identity(raw, path)
  identity <- identity_pack$identity
  if (anyDuplicated(identity$module_uid)) stop("Survival input contains duplicate module identities.", call. = FALSE)
  missing <- setdiff(manifest$module_uid, identity$module_uid)
  extra <- setdiff(identity$module_uid, manifest$module_uid)
  if (length(missing) || (isTRUE(strict) && length(extra))) {
    stop(
      "Survival/manifest coverage mismatch: missing=", length(missing), ", extra=", length(extra),
      call. = FALSE
    )
  }
  raw <- raw[match(manifest$module_uid, identity$module_uid), , drop = FALSE]
  columns <- list(
    n_samples = .triage_find_column(raw, c("n_samples"), "n_samples"),
    n_events = .triage_find_column(raw, c("n_events"), "n_events"),
    matched_gene_fraction = .triage_find_column(raw, c("matched_gene_fraction"), "matched_gene_fraction"),
    cox_hr_per_sd = .triage_find_column(raw, c("cox_hr_per_sd", "cox_HR_per_1SD"), "Cox HR per SD"),
    cox_p = .triage_find_column(raw, c("cox_p", "cox_P"), "Cox p-value"),
    cox_fdr = .triage_find_column(raw, c("cox_fdr", "cox_FDR"), "Cox FDR"),
    cutpoint = .triage_find_column(raw, c("cutpoint", "optimal_cutpoint"), "cutpoint"),
    logrank_p = .triage_find_column(raw, c("logrank_p", "optimal_logrank_P"), "log-rank p-value"),
    logrank_fdr = .triage_find_column(raw, c("logrank_fdr", "optimal_logrank_FDR"), "log-rank FDR"),
    hr_high_vs_low = .triage_find_column(raw, c("hr_high_vs_low", "optimal_HR_high_vs_low"), "high-vs-low HR"),
    survival_direction = .triage_find_column(raw, c("survival_direction", "optimal_direction"), "survival direction")
  )
  result <- data.frame(
    module_uid = manifest$module_uid,
    n_samples = as.integer(.triage_numeric(raw, columns$n_samples)),
    n_events = as.integer(.triage_numeric(raw, columns$n_events)),
    matched_gene_fraction = .triage_numeric(raw, columns$matched_gene_fraction),
    cox_hr_per_sd = .triage_numeric(raw, columns$cox_hr_per_sd),
    cox_p = .triage_numeric(raw, columns$cox_p),
    cox_fdr = .triage_numeric(raw, columns$cox_fdr),
    cutpoint = .triage_numeric(raw, columns$cutpoint),
    logrank_p = .triage_numeric(raw, columns$logrank_p),
    logrank_fdr = .triage_numeric(raw, columns$logrank_fdr),
    hr_high_vs_low = .triage_numeric(raw, columns$hr_high_vs_low),
    survival_direction = trimws(as.character(raw[[columns$survival_direction]])),
    stringsAsFactors = FALSE
  )
  numeric_columns <- setdiff(names(result), c("module_uid", "survival_direction"))
  if (isTRUE(strict) && any(!is.finite(as.matrix(result[, numeric_columns, drop = FALSE])))) {
    stop("Survival input contains non-finite required values.", call. = FALSE)
  }
  result
}

.triage_resolve_output_path <- function(path, table_file) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  .drug_absolute_path(path, dirname(normalizePath(table_file, winslash = "/", mustWork = TRUE)))
}

.triage_sort_candidates <- function(data, sort_by) {
  supported <- list(
    drug_number_desc = function(x) -x$drug_number,
    drug_response_density_desc = function(x) -x$drug_response_density,
    target_probability_desc = function(x) -x$target_subtype_probability,
    probability_margin_desc = function(x) -x$probability_margin,
    logrank_p_asc = function(x) x$logrank_p,
    ml_rank_asc = function(x) x$rank_in_subtype
  )
  unknown <- setdiff(sort_by, names(supported))
  if (length(unknown)) {
    stop("Unsupported candidate sort key(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  keys <- lapply(sort_by, function(name) supported[[name]](data))
  ordering <- do.call(order, c(keys, list(data$module_uid, na.last = TRUE)))
  data[ordering, , drop = FALSE]
}

.triage_stepwise <- function(evidence) {
  gates <- c("ml_gate", "prognosis_gate", "module_size_gate", "drug_response_gate")
  reasons <- paste0(gates, "_reason")
  rows <- vector("list", nrow(evidence) * length(gates))
  index <- 0L
  for (row in seq_len(nrow(evidence))) {
    cumulative <- TRUE
    for (gate_index in seq_along(gates)) {
      index <- index + 1L
      passed <- isTRUE(evidence[[gates[[gate_index]]]][[row]])
      cumulative <- cumulative && passed
      rows[[index]] <- data.frame(
        module_uid = evidence$module_uid[[row]],
        subtype = evidence$subtype[[row]],
        gate_order = gate_index,
        gate_name = gates[[gate_index]],
        gate_pass = passed,
        cumulative_pass = cumulative,
        gate_reason = evidence[[reasons[[gate_index]]]][[row]],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#' Integrate ML Top-K, survival and drug evidence and select one module/subtype
#'
#' Step 6C starts only from `ml_top*.tsv`, applies prognosis and module-size
#' gates, joins the configured primary drug panel, and ranks surviving modules
#' by explicit sort keys. It writes a complete evidence table, gate-by-gate
#' audit trail and exactly one self-contained selected module per subtype.
#'
#' @param module_manifest_file Step-4 manifest.
#' @param ml_top_file Step-6B `ml_top*.tsv`.
#' @param survival_file All-module survival table. Both canonical columns and
#'   the current LUAD `optimal_*`/`cox_*` names are accepted.
#' @param drug_response_summary_file Step-6 summary table.
#' @param output_dir New step-6C output directory.
#' @param primary_drug_panel Panel used for candidate ranking.
#' @param min_module_size Minimum module size, inclusive.
#' @param required_direction Required survival direction.
#' @param logrank_p_max Maximum nominal log-rank p-value.
#' @param logrank_fdr_max Optional maximum log-rank FDR.
#' @param sort_by Ordered candidate sort keys.
#' @param select_n_per_subtype Must be one for the selected-only boundary.
#' @param strict Require complete inputs and at least one passing candidate for
#'   every subtype represented in `ml_top_file`.
#' @return The selected-module table with output paths in attributes.
#' @export
triage_modules <- function(
  module_manifest_file,
  ml_top_file,
  survival_file,
  drug_response_summary_file,
  output_dir,
  primary_drug_panel = "PRISM",
  min_module_size = 30L,
  required_direction = "High_score_worse",
  logrank_p_max = 0.05,
  logrank_fdr_max = NULL,
  sort_by = c("drug_number_desc", "drug_response_density_desc", "target_probability_desc"),
  select_n_per_subtype = 1L,
  strict = TRUE
) {
  if (select_n_per_subtype != 1L) {
    stop("select_n_per_subtype must be 1 to protect the selected-only downstream boundary.", call. = FALSE)
  }
  if (length(min_module_size) != 1L || !is.numeric(min_module_size) ||
      is.na(min_module_size) || min_module_size < 1 || min_module_size != floor(min_module_size)) {
    stop("min_module_size must be a positive integer.", call. = FALSE)
  }
  for (item in list(logrank_p_max = logrank_p_max, logrank_fdr_max = logrank_fdr_max)) {
    if (!is.null(item) && (length(item) != 1L || !is.numeric(item) || is.na(item) || item < 0 || item > 1)) {
      stop("Log-rank thresholds must be NULL or values between 0 and 1.", call. = FALSE)
    }
  }
  manifest <- read_module_manifest(module_manifest_file, check_files = TRUE, verify_hashes = strict)
  ml_top <- .annotation_read_tabular(ml_top_file)
  ml_required <- c(
    "module_uid", "target_subtype", "target_subtype_probability", "probability_margin",
    "rank_in_subtype", "ml_gate", "ml_gate_reason"
  )
  missing_ml <- setdiff(ml_required, names(ml_top))
  if (length(missing_ml)) stop("ML Top-K input is missing: ", paste(missing_ml, collapse = ", "), call. = FALSE)
  if (anyDuplicated(ml_top$module_uid)) stop("ML Top-K input contains duplicate module_uid values.", call. = FALSE)
  if (any(!ml_top$module_uid %in% manifest$module_uid)) stop("ML Top-K contains modules absent from the manifest.", call. = FALSE)
  parsed_top <- parse_module_uid(ml_top$module_uid)
  if (any(parsed_top$subtype != ml_top$target_subtype)) {
    stop("ML target_subtype must match module_uid subtype.", call. = FALSE)
  }

  survival_raw <- .annotation_read_tabular(survival_file)
  survival <- .triage_standardize_survival(survival_raw, manifest, survival_file, strict)
  drug <- .annotation_read_tabular(drug_response_summary_file)
  drug_required <- c(
    "module_uid", "drug_panel", "drug_number", "tested_drug_number",
    "drug_response_density", "drn_edge_number", "drn_file", "drn_info_file",
    "analysis_status", "status_reason"
  )
  missing_drug <- setdiff(drug_required, names(drug))
  if (length(missing_drug)) stop("Drug-response summary is missing: ", paste(missing_drug, collapse = ", "), call. = FALSE)
  drug$drug_panel <- toupper(trimws(as.character(drug$drug_panel)))
  primary_drug_panel <- toupper(trimws(primary_drug_panel))
  drug <- drug[drug$drug_panel == primary_drug_panel, , drop = FALSE]
  if (anyDuplicated(drug$module_uid)) stop("Primary drug panel has duplicate module_uid rows.", call. = FALSE)
  if (isTRUE(strict) && !setequal(drug$module_uid, manifest$module_uid)) {
    stop("Primary drug panel does not cover every manifest module.", call. = FALSE)
  }

  manifest_subset <- manifest[match(ml_top$module_uid, manifest$module_uid), , drop = FALSE]
  survival_subset <- survival[match(ml_top$module_uid, survival$module_uid), , drop = FALSE]
  drug_subset <- drug[match(ml_top$module_uid, drug$module_uid), , drop = FALSE]
  if (isTRUE(strict) && (anyNA(survival_subset$module_uid) || anyNA(drug_subset$module_uid))) {
    stop("ML Top-K modules are missing survival or primary drug evidence.", call. = FALSE)
  }
  evidence <- cbind(
    manifest_subset[, c("module_uid", "network", "method", "subtype", "module", "module_size"), drop = FALSE],
    ml_top[, setdiff(names(ml_top), c("module_uid", "network", "method", "subtype", "module")), drop = FALSE],
    survival_subset[, setdiff(names(survival_subset), "module_uid"), drop = FALSE],
    data.frame(
      primary_drug_panel = drug_subset$drug_panel,
      drug_number = as.integer(drug_subset$drug_number),
      tested_drug_number = as.integer(drug_subset$tested_drug_number),
      drug_response_density = as.numeric(drug_subset$drug_response_density),
      drn_edge_number = as.integer(drug_subset$drn_edge_number),
      source_drn_file = vapply(
        drug_subset$drn_file,
        .triage_resolve_output_path,
        character(1),
        table_file = drug_response_summary_file
      ),
      source_drn_info_file = vapply(
        drug_subset$drn_info_file,
        .triage_resolve_output_path,
        character(1),
        table_file = drug_response_summary_file
      ),
      drug_analysis_status = drug_subset$analysis_status,
      drug_status_reason = drug_subset$status_reason,
      stringsAsFactors = FALSE
    )
  )
  evidence$ml_gate <- as.logical(evidence$ml_gate)
  evidence$ml_gate_reason <- ifelse(evidence$ml_gate, "ML Top-K eligible", evidence$ml_gate_reason)
  evidence$prognosis_gate <- evidence$survival_direction == required_direction &
    evidence$logrank_p <= logrank_p_max
  if (!is.null(logrank_fdr_max)) {
    evidence$prognosis_gate <- evidence$prognosis_gate & evidence$logrank_fdr <= logrank_fdr_max
  }
  evidence$prognosis_gate_reason <- ifelse(
    evidence$prognosis_gate,
    "required_direction_and_logrank_threshold_passed",
    paste0(
      "direction=", evidence$survival_direction,
      ";logrank_p=", signif(evidence$logrank_p, 5),
      if (!is.null(logrank_fdr_max)) paste0(";logrank_fdr=", signif(evidence$logrank_fdr, 5)) else ""
    )
  )
  evidence$module_size_gate <- evidence$module_size >= min_module_size
  evidence$module_size_gate_reason <- ifelse(
    evidence$module_size_gate,
    paste0("module_size>=", min_module_size),
    paste0("module_size<", min_module_size)
  )
  evidence$drug_response_gate <- evidence$drug_analysis_status == "completed" &
    is.finite(evidence$drug_number) & is.finite(evidence$drug_response_density)
  evidence$drug_response_gate_reason <- ifelse(
    evidence$drug_response_gate,
    paste0("primary_panel=", primary_drug_panel, ";drug_evidence_available"),
    paste0("primary_panel=", primary_drug_panel, ";", evidence$drug_analysis_status)
  )
  evidence$all_gates_pass <- evidence$ml_gate & evidence$prognosis_gate &
    evidence$module_size_gate & evidence$drug_response_gate

  selected_rows <- list()
  for (subtype in sort(unique(evidence$subtype))) {
    candidates <- evidence[evidence$subtype == subtype & evidence$all_gates_pass, , drop = FALSE]
    if (!nrow(candidates)) {
      if (isTRUE(strict)) stop("No module passed all triage gates for subtype ", subtype, ".", call. = FALSE)
      next
    }
    candidates <- .triage_sort_candidates(candidates, sort_by)
    selected_rows[[subtype]] <- candidates[1L, , drop = FALSE]
  }
  if (!length(selected_rows)) stop("No subtype has a selected module.", call. = FALSE)
  selected_evidence <- do.call(rbind, selected_rows)
  rownames(selected_evidence) <- NULL
  stepwise <- .triage_stepwise(evidence)

  output_target <- .annotation_output_target(output_dir)
  staging_dir <- tempfile(
    pattern = paste0(".", basename(output_target), "-"),
    tmpdir = dirname(output_target)
  )
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create module-triage staging directory.", call. = FALSE)
  }
  on.exit({
    if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  selected_rows_out <- vector("list", nrow(selected_evidence))
  for (index in seq_len(nrow(selected_evidence))) {
    row <- selected_evidence[index, , drop = FALSE]
    source_manifest <- manifest[manifest$module_uid == row$module_uid, , drop = FALSE]
    relative_dir <- file.path("selected", row$subtype, row$module_uid)
    target_dir <- file.path(staging_dir, relative_dir)
    if (!dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create selected-module directory: ", target_dir, call. = FALSE)
    }
    sources <- c(
      node_file = source_manifest$node_file_abs,
      edge_file = source_manifest$edge_file_abs,
      drn_file = row$source_drn_file,
      drn_info_file = row$source_drn_info_file
    )
    if (any(!file.exists(sources))) {
      stop("Selected module references missing node/edge/DRN files: ", row$module_uid, call. = FALSE)
    }
    targets <- c(
      node_file = file.path(target_dir, "nodes.tsv"),
      edge_file = file.path(target_dir, "edges.tsv"),
      drn_file = file.path(target_dir, "drn.tsv"),
      drn_info_file = file.path(target_dir, "drn_info.tsv")
    )
    copied <- mapply(file.copy, sources, targets, MoreArgs = list(overwrite = FALSE), USE.NAMES = FALSE)
    if (!all(copied)) stop("Could not copy selected module inputs: ", row$module_uid, call. = FALSE)
    relative_paths <- gsub("\\\\", "/", file.path(relative_dir, basename(targets)))
    names(relative_paths) <- names(targets)
    selected_rows_out[[index]] <- data.frame(
      module_uid = row$module_uid,
      network = row$network,
      method = row$method,
      subtype = row$subtype,
      module = row$module,
      module_size = row$module_size,
      primary_drug_panel = primary_drug_panel,
      node_file = relative_paths[["node_file"]],
      edge_file = relative_paths[["edge_file"]],
      drn_file = relative_paths[["drn_file"]],
      drn_info_file = relative_paths[["drn_info_file"]],
      target_subtype_probability = row$target_subtype_probability,
      probability_margin = row$probability_margin,
      ml_rank_in_subtype = row$rank_in_subtype,
      survival_direction = row$survival_direction,
      logrank_p = row$logrank_p,
      logrank_fdr = row$logrank_fdr,
      hr_high_vs_low = row$hr_high_vs_low,
      drug_number = row$drug_number,
      tested_drug_number = row$tested_drug_number,
      drug_response_density = row$drug_response_density,
      drn_edge_number = row$drn_edge_number,
      selection_rank = 1L,
      selection_reason = paste0(
        "ML Top-K; prognosis ", required_direction, " with logrank_p<=", logrank_p_max,
        "; module_size>=", min_module_size, "; ranked by ", paste(sort_by, collapse = ">")
      ),
      stringsAsFactors = FALSE
    )
  }
  selected <- do.call(rbind, selected_rows_out)
  rownames(selected) <- NULL

  survival_file_out <- file.path(staging_dir, "module_survival.tsv")
  evidence_file <- file.path(staging_dir, "module_evidence.tsv")
  stepwise_file <- file.path(staging_dir, "module_filtering_stepwise.tsv")
  selected_file <- file.path(staging_dir, "selected_modules.tsv")
  .manifest_write_tsv(survival, survival_file_out)
  .manifest_write_tsv(evidence, evidence_file)
  .manifest_write_tsv(stepwise, stepwise_file)
  .manifest_write_tsv(selected, selected_file)

  if (!file.rename(staging_dir, output_target)) {
    stop("Could not finalize module-triage output directory: ", output_target, call. = FALSE)
  }
  attr(selected, "output_dir") <- output_target
  attr(selected, "selected_file") <- file.path(output_target, basename(selected_file))
  attr(selected, "survival_file") <- file.path(output_target, basename(survival_file_out))
  attr(selected, "evidence_file") <- file.path(output_target, basename(evidence_file))
  attr(selected, "stepwise_file") <- file.path(output_target, basename(stepwise_file))
  selected
}
