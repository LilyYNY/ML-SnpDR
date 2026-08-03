# Binding scores -------------------------------------------------------------

.binding_score_input <- function(path) {
  data <- .annotation_read_tabular(path)
  lowered <- tolower(names(data))
  aliases <- list(
    module_uid = c("module_uid"),
    target_name = c("target_name", "target name", "protein", "target"),
    drug_name = c("drug_name", "drug name", "drug", "compound"),
    binding_score = c("binding_score", "binding score", "score", "prediction")
  )
  selected <- vapply(aliases, function(options) {
    hit <- match(options, lowered, nomatch = 0L)
    hit <- hit[hit > 0L]
    if (length(hit)) names(data)[hit[[1L]]] else NA_character_
  }, character(1))
  if (anyNA(selected)) {
    stop("Binding-score input must contain module_uid, target, drug and binding score.", call. = FALSE)
  }
  result <- data.frame(
    module_uid = trimws(as.character(data[[selected[["module_uid"]]]])),
    target_name = trimws(as.character(data[[selected[["target_name"]]]])),
    drug_name = trimws(as.character(data[[selected[["drug_name"]]]])),
    binding_score = suppressWarnings(as.numeric(as.character(data[[selected[["binding_score"]]]]))),
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(result$binding_score)) || any(!nzchar(result$target_name)) || any(!nzchar(result$drug_name))) {
    stop("Binding-score input contains invalid values.", call. = FALSE)
  }
  if (anyDuplicated(paste(result$module_uid, result$target_name, result$drug_name, sep = "\r"))) {
    stop("Binding-score input contains duplicate module/target/drug rows.", call. = FALSE)
  }
  result
}

.binding_read_seq_manifest <- function(path, selected) {
  manifest <- .annotation_read_tabular(path)
  required <- c("module_uid", "sequence_file", "smiles_file", "analysis_status")
  missing <- setdiff(required, names(manifest))
  if (length(missing)) stop("seq_smiles_manifest.tsv is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyDuplicated(manifest$module_uid)) stop("Sequence/SMILES manifest contains duplicate module_uid values.", call. = FALSE)
  if (!setequal(manifest$module_uid, selected$module_uid)) {
    stop("Sequence/SMILES manifest must cover exactly the selected modules.", call. = FALSE)
  }
  root <- dirname(normalizePath(path, winslash = "/", mustWork = TRUE))
  manifest$sequence_file_abs <- vapply(manifest$sequence_file, .drug_absolute_path, character(1), root = root)
  manifest$smiles_file_abs <- vapply(manifest$smiles_file, .drug_absolute_path, character(1), root = root)
  if (any(!file.exists(c(manifest$sequence_file_abs, manifest$smiles_file_abs)))) {
    stop("Sequence/SMILES manifest references missing files.", call. = FALSE)
  }
  manifest[match(selected$module_uid, manifest$module_uid), , drop = FALSE]
}

#' Predict or import binding scores for selected modules only
#'
#' Step 8 constructs drug-target pairs from each selected DRN and joins the
#' step-7 sequences/SMILES. Scores may be supplied in a tabular file (for
#' example DeepPurpose output) or calculated by a user function receiving the
#' module's DPI table. This keeps the selected-module boundary independent of
#' the Python environment.
#'
#' @param selected_modules_file Step-6C `selected_modules.tsv`.
#' @param seq_smiles_manifest_file Step-7 manifest.
#' @param output_base New step-8 output directory.
#' @param binding_score_file Optional all-selected-module score table with
#'   `module_uid`, target/protein, drug and score columns.
#' @param predictor Optional function called once per module with a data frame
#'   containing `module_uid`, `target_name`, `drug_name`, `target_seq` and
#'   `drug_smiles`. It must return one numeric score per row or a data frame
#'   containing `binding_score`.
#' @param score_direction Whether lower or higher binding scores rank first.
#' @param strict Require one score for every selected DRN pair.
#' @return Unified `binding_scores.tsv` with output paths in attributes.
#' @export
predict_BA <- function(
  selected_modules_file,
  seq_smiles_manifest_file,
  output_base,
  binding_score_file = NULL,
  predictor = NULL,
  score_direction = c("lower_better", "higher_better"),
  strict = TRUE
) {
  score_direction <- match.arg(score_direction)
  if (!xor(is.null(binding_score_file), is.null(predictor))) {
    stop("Supply exactly one of binding_score_file or predictor.", call. = FALSE)
  }
  if (!is.null(predictor) && !is.function(predictor)) stop("predictor must be a function.", call. = FALSE)
  selected <- .selected_read(selected_modules_file, verify_files = TRUE)
  seq_manifest <- .binding_read_seq_manifest(seq_smiles_manifest_file, selected)
  supplied <- if (!is.null(binding_score_file)) .binding_score_input(binding_score_file) else NULL
  if (!is.null(supplied) && any(!supplied$module_uid %in% selected$module_uid)) {
    stop("Binding-score input contains modules outside selected_modules.tsv.", call. = FALSE)
  }

  output_target <- .annotation_output_target(output_base)
  staging_dir <- tempfile(pattern = paste0(".", basename(output_target), "-"), tmpdir = dirname(output_target))
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create binding-score staging directory.", call. = FALSE)
  }
  on.exit({
    if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  score_rows <- vector("list", nrow(selected))
  manifest_rows <- vector("list", nrow(selected))
  for (index in seq_len(nrow(selected))) {
    module <- selected[index, , drop = FALSE]
    drn <- .drug_read_drn(module$drn_file_abs, "none", 1)
    pairs <- unique(drn[, c("protein", "drug"), drop = FALSE])
    names(pairs) <- c("target_name", "drug_name")
    sequences <- .annotation_read_tabular(seq_manifest$sequence_file_abs[[index]])
    smiles <- .annotation_read_tabular(seq_manifest$smiles_file_abs[[index]])
    if (!all(c("node", "sequence") %in% names(sequences)) ||
        !all(c("node", "SMILES") %in% names(smiles))) {
      stop("Step-7 module files have invalid columns for ", module$module_uid, ".", call. = FALSE)
    }
    dpi <- data.frame(
      module_uid = module$module_uid,
      target_name = pairs$target_name,
      drug_name = pairs$drug_name,
      target_seq = sequences$sequence[match(pairs$target_name, sequences$node)],
      drug_smiles = smiles$SMILES[match(pairs$drug_name, smiles$node)],
      stringsAsFactors = FALSE
    )
    complete <- stats::complete.cases(dpi[, c("target_seq", "drug_smiles")]) &
      nzchar(dpi$target_seq) & nzchar(dpi$drug_smiles)
    if (isTRUE(strict) && any(!complete)) {
      stop("DPI sequence/SMILES coverage is incomplete for ", module$module_uid, ".", call. = FALSE)
    }
    dpi <- dpi[complete, , drop = FALSE]
    if (!nrow(dpi)) stop("No scoreable DPI pairs for ", module$module_uid, ".", call. = FALSE)

    if (!is.null(supplied)) {
      module_scores <- supplied[supplied$module_uid == module$module_uid, , drop = FALSE]
      pair_key <- paste(dpi$target_name, dpi$drug_name, sep = "\r")
      score_key <- paste(module_scores$target_name, module_scores$drug_name, sep = "\r")
      score_index <- match(pair_key, score_key)
      if (isTRUE(strict) && anyNA(score_index)) {
        stop("Binding-score file does not cover every selected DRN pair for ", module$module_uid, ".", call. = FALSE)
      }
      dpi$binding_score <- module_scores$binding_score[score_index]
      dpi <- dpi[is.finite(dpi$binding_score), , drop = FALSE]
    } else {
      prediction <- predictor(dpi)
      values <- if (is.data.frame(prediction) && "binding_score" %in% names(prediction)) {
        prediction$binding_score
      } else {
        prediction
      }
      values <- suppressWarnings(as.numeric(values))
      if (length(values) != nrow(dpi) || any(!is.finite(values))) {
        stop("predictor must return one finite binding score per DPI row.", call. = FALSE)
      }
      dpi$binding_score <- values
    }
    ordering <- if (score_direction == "lower_better") {
      order(dpi$binding_score, dpi$target_name, dpi$drug_name)
    } else {
      order(-dpi$binding_score, dpi$target_name, dpi$drug_name)
    }
    dpi <- dpi[ordering, , drop = FALSE]
    dpi$binding_rank <- seq_len(nrow(dpi))
    dpi$score_direction <- score_direction
    score_rows[[index]] <- dpi

    relative_dir <- file.path("selected", module$subtype, module$module_uid)
    module_dir <- file.path(staging_dir, relative_dir)
    if (!dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create binding module directory.", call. = FALSE)
    }
    dpi_relative <- gsub("\\\\", "/", file.path(relative_dir, "dpi_info.tsv"))
    score_relative <- gsub("\\\\", "/", file.path(relative_dir, "binding_scores.tsv"))
    .manifest_write_tsv(dpi[, setdiff(names(dpi), c("binding_score", "binding_rank", "score_direction")), drop = FALSE], file.path(staging_dir, dpi_relative))
    .manifest_write_tsv(dpi, file.path(staging_dir, score_relative))
    manifest_rows[[index]] <- data.frame(
      module_uid = module$module_uid,
      subtype = module$subtype,
      dpi_pair_number = nrow(dpi),
      dpi_info_file = dpi_relative,
      binding_score_file = score_relative,
      analysis_status = "completed",
      status_reason = "",
      stringsAsFactors = FALSE
    )
  }
  scores <- do.call(rbind, score_rows)
  rownames(scores) <- NULL
  binding_manifest <- do.call(rbind, manifest_rows)
  scores_file <- file.path(staging_dir, "binding_scores.tsv")
  manifest_file <- file.path(staging_dir, "binding_score_manifest.tsv")
  .manifest_write_tsv(scores, scores_file)
  .manifest_write_tsv(binding_manifest, manifest_file)
  if (!file.rename(staging_dir, output_target)) {
    stop("Could not finalize binding-score output directory: ", output_target, call. = FALSE)
  }
  attr(scores, "output_dir") <- output_target
  attr(scores, "scores_file") <- file.path(output_target, basename(scores_file))
  attr(scores, "manifest_file") <- file.path(output_target, basename(manifest_file))
  scores
}
