# Perturbation scores --------------------------------------------------------

.perturbation_sensitivity_input <- function(path) {
  data <- .annotation_read_tabular(path)
  lowered <- tolower(names(data))
  aliases <- list(
    module_uid = c("module_uid"),
    target_name = c("target_name", "target.name", "target name", "protein", "orf_name"),
    sensitivity = c("sensitivity", "sens", "prs_sensitivity")
  )
  selected <- vapply(aliases, function(options) {
    hit <- match(options, lowered, nomatch = 0L)
    hit <- hit[hit > 0L]
    if (length(hit)) names(data)[hit[[1L]]] else NA_character_
  }, character(1))
  if (anyNA(selected)) {
    stop("Sensitivity input must contain module_uid, target and sensitivity.", call. = FALSE)
  }
  result <- data.frame(
    module_uid = trimws(as.character(data[[selected[["module_uid"]]]])),
    target_name = trimws(as.character(data[[selected[["target_name"]]]])),
    sensitivity = suppressWarnings(as.numeric(as.character(data[[selected[["sensitivity"]]]]))),
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(result$sensitivity)) || any(!nzchar(result$target_name))) {
    stop("Sensitivity input contains invalid values.", call. = FALSE)
  }
  if (anyDuplicated(paste(result$module_uid, result$target_name, sep = "\r"))) {
    stop("Sensitivity input contains duplicate module/target rows.", call. = FALSE)
  }
  result
}

.perturbation_read_edges <- function(path) {
  edges <- .annotation_read_tabular(path)
  if (!all(c("node1", "node2") %in% names(edges))) {
    if (all(c("gene1", "gene2") %in% names(edges))) {
      names(edges)[match(c("gene1", "gene2"), names(edges))] <- c("node1", "node2")
    } else {
      stop("Selected module edge file must contain node1/node2 or gene1/gene2.", call. = FALSE)
    }
  }
  unique(edges[, c("node1", "node2"), drop = FALSE])
}

#' Combine binding and perturbation sensitivity scores for selected modules
#'
#' Step 9 reads only selected-module edge files and step-8 binding scores.
#' Target sensitivity can be imported from an ENM/PRS result table or returned
#' by a user function receiving `(edges, selected_module_row)`. The upstream
#' perturbation score definition, binding score multiplied by sensitivity, is
#' retained.
#'
#' @param selected_modules_file Step-6C selected modules.
#' @param binding_scores_file Step-8 unified binding scores.
#' @param output_base New step-9 output directory.
#' @param sensitivity_file Optional table with `module_uid`, target/protein and
#'   sensitivity columns.
#' @param sensitivity_function Optional function called per selected module
#'   with its edge table and selected-module row; it must return target and
#'   sensitivity columns.
#' @param top_n Number of final drug-target rows retained per subtype.
#' @param strict Require sensitivity for every target in the binding table.
#' @return Unified perturbation-score table with final-candidate paths in
#'   attributes.
#' @export
process_prs_dti <- function(
  selected_modules_file,
  binding_scores_file,
  output_base,
  sensitivity_file = NULL,
  sensitivity_function = NULL,
  top_n = 10L,
  strict = TRUE
) {
  if (!xor(is.null(sensitivity_file), is.null(sensitivity_function))) {
    stop("Supply exactly one of sensitivity_file or sensitivity_function.", call. = FALSE)
  }
  if (!is.null(sensitivity_function) && !is.function(sensitivity_function)) {
    stop("sensitivity_function must be a function.", call. = FALSE)
  }
  if (length(top_n) != 1L || !is.numeric(top_n) || is.na(top_n) || top_n < 1 || top_n != floor(top_n)) {
    stop("top_n must be a positive integer.", call. = FALSE)
  }
  selected <- .selected_read(selected_modules_file, verify_files = TRUE)
  binding <- .annotation_read_tabular(binding_scores_file)
  binding_required <- c(
    "module_uid", "target_name", "drug_name", "binding_score",
    "binding_rank", "score_direction"
  )
  missing <- setdiff(binding_required, names(binding))
  if (length(missing)) stop("Binding scores are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (any(!binding$module_uid %in% selected$module_uid)) {
    stop("Binding scores contain modules outside selected_modules.tsv.", call. = FALSE)
  }
  if (isTRUE(strict) && !setequal(unique(binding$module_uid), selected$module_uid)) {
    stop("Binding scores do not cover every selected module.", call. = FALSE)
  }
  binding$binding_score <- suppressWarnings(as.numeric(binding$binding_score))
  if (any(!is.finite(binding$binding_score))) stop("Binding scores contain non-finite values.", call. = FALSE)

  sensitivity <- if (!is.null(sensitivity_file)) {
    .perturbation_sensitivity_input(sensitivity_file)
  } else {
    rows <- vector("list", nrow(selected))
    for (index in seq_len(nrow(selected))) {
      module <- selected[index, , drop = FALSE]
      edges <- .perturbation_read_edges(module$edge_file_abs)
      value <- sensitivity_function(edges, module)
      if (!is.data.frame(value)) stop("sensitivity_function must return a data frame.", call. = FALSE)
      lowered <- tolower(names(value))
      target_hit <- match(c("target_name", "protein", "node", "orf_name"), lowered, nomatch = 0L)
      target_hit <- target_hit[target_hit > 0L]
      sensitivity_hit <- match(c("sensitivity", "sens", "prs_sensitivity"), lowered, nomatch = 0L)
      sensitivity_hit <- sensitivity_hit[sensitivity_hit > 0L]
      if (!length(target_hit) || !length(sensitivity_hit)) {
        stop("sensitivity_function output must contain target and sensitivity columns.", call. = FALSE)
      }
      rows[[index]] <- data.frame(
        module_uid = module$module_uid,
        target_name = trimws(as.character(value[[target_hit[[1L]]]])),
        sensitivity = suppressWarnings(as.numeric(value[[sensitivity_hit[[1L]]]])),
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }
  if (any(!sensitivity$module_uid %in% selected$module_uid)) {
    stop("Sensitivity input contains modules outside selected_modules.tsv.", call. = FALSE)
  }
  if (anyDuplicated(paste(sensitivity$module_uid, sensitivity$target_name, sep = "\r"))) {
    stop("Sensitivity input contains duplicate module/target rows.", call. = FALSE)
  }

  binding_key <- paste(binding$module_uid, binding$target_name, sep = "\r")
  sensitivity_key <- paste(sensitivity$module_uid, sensitivity$target_name, sep = "\r")
  sensitivity_index <- match(binding_key, sensitivity_key)
  if (isTRUE(strict) && anyNA(sensitivity_index)) {
    first <- which(is.na(sensitivity_index))[[1L]]
    stop(
      "Sensitivity is missing for ", binding$module_uid[[first]], "/", binding$target_name[[first]], ".",
      call. = FALSE
    )
  }
  result <- binding[!is.na(sensitivity_index), , drop = FALSE]
  sensitivity_index <- sensitivity_index[!is.na(sensitivity_index)]
  result$sensitivity <- sensitivity$sensitivity[sensitivity_index]
  if (any(!is.finite(result$sensitivity))) stop("Sensitivity contains non-finite values.", call. = FALSE)
  result$perturbation_score <- result$binding_score * result$sensitivity
  result$perturbation_rank_in_module <- NA_integer_
  for (uid in unique(result$module_uid)) {
    indices <- which(result$module_uid == uid)
    ordering <- order(-result$perturbation_score[indices], result$binding_rank[indices], result$drug_name[indices])
    result$perturbation_rank_in_module[indices[ordering]] <- seq_along(indices)
  }
  selected_identity <- selected[, c("module_uid", "subtype", "network", "method", "module"), drop = FALSE]
  result <- merge(selected_identity, result, by = "module_uid", all.y = TRUE, sort = FALSE)
  result <- result[order(result$subtype, result$perturbation_rank_in_module, result$module_uid), , drop = FALSE]
  rownames(result) <- NULL

  final_rows <- lapply(split(result, result$subtype), function(group) {
    group <- group[order(-group$perturbation_score, group$binding_rank, group$module_uid), , drop = FALSE]
    utils::head(group, as.integer(top_n))
  })
  final <- do.call(rbind, final_rows)
  rownames(final) <- NULL
  final$final_rank_in_subtype <- stats::ave(
    -final$perturbation_score,
    final$subtype,
    FUN = function(x) rank(x, ties.method = "first")
  )

  output_target <- .annotation_output_target(output_base)
  staging_dir <- tempfile(pattern = paste0(".", basename(output_target), "-"), tmpdir = dirname(output_target))
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create perturbation-score staging directory.", call. = FALSE)
  }
  on.exit({
    if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  scores_file <- file.path(staging_dir, "perturbation_scores.tsv")
  final_file <- file.path(staging_dir, "final_candidates.tsv")
  sensitivity_file_out <- file.path(staging_dir, "target_sensitivity.tsv")
  .manifest_write_tsv(result, scores_file)
  .manifest_write_tsv(final, final_file)
  .manifest_write_tsv(sensitivity, sensitivity_file_out)
  if (!file.rename(staging_dir, output_target)) {
    stop("Could not finalize perturbation-score output directory: ", output_target, call. = FALSE)
  }
  attr(result, "output_dir") <- output_target
  attr(result, "scores_file") <- file.path(output_target, basename(scores_file))
  attr(result, "final_file") <- file.path(output_target, basename(final_file))
  attr(result, "sensitivity_file") <- file.path(output_target, basename(sensitivity_file_out))
  result
}
