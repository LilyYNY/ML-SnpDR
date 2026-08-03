# Drug response --------------------------------------------------------------

.drug_absolute_path <- function(path, root) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  absolute_pattern <- "^(?:[A-Za-z]:[/\\\\]|/|\\\\\\\\)"
  resolved <- if (grepl(absolute_pattern, path)) path else file.path(root, path)
  normalizePath(resolved, winslash = "/", mustWork = FALSE)
}

.drug_standardize_index <- function(index, root = getwd()) {
  required <- c("module_uid", "drug_panel", "drn_file", "drn_info_file")
  missing <- setdiff(required, names(index))
  if (length(missing)) {
    stop(
      "Drug-response input index is missing column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  optional <- c("drug_level_file", "prediction_file")
  for (column in optional) if (!column %in% names(index)) index[[column]] <- NA_character_
  index <- index[, c(required, optional), drop = FALSE]
  index[] <- lapply(index, as.character)
  index$module_uid <- trimws(index$module_uid)
  index$drug_panel <- toupper(trimws(index$drug_panel))
  for (column in c("drn_file", "drn_info_file", optional)) {
    index[[column]] <- vapply(index[[column]], .drug_absolute_path, character(1), root = root)
  }
  if (any(!nzchar(index$module_uid)) || any(!nzchar(index$drug_panel))) {
    stop("module_uid and drug_panel cannot be blank in the drug-response index.", call. = FALSE)
  }
  key <- paste(index$module_uid, index$drug_panel, sep = "\r")
  if (anyDuplicated(key)) {
    stop("Drug-response input index contains duplicate module_uid/drug_panel rows.", call. = FALSE)
  }
  index
}

.drug_index_from_roots <- function(manifest, drug_response_roots) {
  if (is.list(drug_response_roots)) {
    roots <- unlist(drug_response_roots, use.names = TRUE)
  } else {
    roots <- drug_response_roots
  }
  if (!is.character(roots) || !length(roots) || is.null(names(roots)) ||
      any(!nzchar(names(roots))) || anyDuplicated(toupper(names(roots)))) {
    stop(
      "drug_response_roots must be a named character vector/list, for example ",
      "c(PRISM = 'path/to/DrugResponse_PRISM').",
      call. = FALSE
    )
  }
  panels <- toupper(names(roots))
  roots <- vapply(roots, function(path) {
    if (length(path) != 1L || is.na(path) || !dir.exists(path)) {
      stop("Drug-response root does not exist: ", path, call. = FALSE)
    }
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }, character(1))
  names(roots) <- panels

  rows <- vector("list", nrow(manifest) * length(roots))
  index <- 0L
  for (manifest_index in seq_len(nrow(manifest))) {
    module <- manifest[manifest_index, , drop = FALSE]
    method_file <- tolower(module$method)
    network_file <- tolower(module$network)
    prefix <- paste(method_file, network_file, module$subtype, module$module, sep = "_")
    for (panel in names(roots)) {
      index <- index + 1L
      module_dir <- file.path(
        roots[[panel]],
        module$subtype,
        module$method,
        network_file,
        module$module
      )
      rows[[index]] <- data.frame(
        module_uid = module$module_uid,
        drug_panel = panel,
        drn_file = file.path(module_dir, paste0("DRN_", prefix, ".txt")),
        drn_info_file = file.path(module_dir, paste0("DRN_info_", prefix, ".txt")),
        drug_level_file = file.path(module_dir, paste0("drug_response_level_", prefix, ".csv")),
        prediction_file = file.path(module_dir, "calcPhenotype_Output", "DrugPredictions.csv"),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

.drug_read_index <- function(path) {
  raw <- .annotation_read_tabular(path)
  .drug_standardize_index(
    raw,
    root = dirname(normalizePath(path, winslash = "/", mustWork = TRUE))
  )
}

.drug_missing_files <- function(index) {
  required <- c("drn_file", "drn_info_file", "drug_level_file", "prediction_file")
  missing <- lapply(seq_len(nrow(index)), function(i) {
    fields <- required[!vapply(index[i, required, drop = FALSE], function(path) {
      length(path) == 1L && !is.na(path) && file.exists(path) && !dir.exists(path)
    }, logical(1))]
    if (!length(fields)) return(NULL)
    data.frame(
      module_uid = index$module_uid[[i]],
      drug_panel = index$drug_panel[[i]],
      missing_fields = paste(fields, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  missing <- missing[!vapply(missing, is.null, logical(1))]
  if (length(missing)) do.call(rbind, missing) else data.frame()
}

.drug_read_drn <- function(path, p_adjust_method, significance_cutoff) {
  raw <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  aliases <- list(
    protein = c("protein", "target", "gene"),
    drug = c("drug", "compound"),
    p_value = c("p_value", "pvalue", "p.value")
  )
  lowered <- tolower(names(raw))
  selected <- vapply(aliases, function(options) {
    hit <- match(options, lowered, nomatch = 0L)
    hit <- hit[hit > 0L]
    if (length(hit)) names(raw)[hit[[1L]]] else NA_character_
  }, character(1))
  if (anyNA(selected)) {
    stop("DRN file must contain protein, drug and pvalue/p_value columns: ", path, call. = FALSE)
  }
  drn <- data.frame(
    protein = trimws(as.character(raw[[selected[["protein"]]]])),
    drug = trimws(as.character(raw[[selected[["drug"]]]])),
    p_value = suppressWarnings(as.numeric(raw[[selected[["p_value"]]]])),
    stringsAsFactors = FALSE
  )
  if (nrow(drn) && (any(!is.finite(drn$p_value)) || any(drn$p_value < 0 | drn$p_value > 1) ||
                    any(!nzchar(drn$protein)) || any(!nzchar(drn$drug)))) {
    stop("DRN contains invalid protein, drug or p-value values: ", path, call. = FALSE)
  }
  drn$p_adjust <- if (identical(p_adjust_method, "none")) {
    drn$p_value
  } else {
    stats::p.adjust(drn$p_value, method = p_adjust_method)
  }
  drn$significant <- drn$p_adjust <= significance_cutoff
  drn <- unique(drn)
  drn[order(drn$p_adjust, drn$p_value, drn$protein, drn$drug), , drop = FALSE]
}

.drug_read_info <- function(path) {
  info <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  if (!all(c("node", "type") %in% names(info))) {
    stop("DRN info file must contain node and type columns: ", path, call. = FALSE)
  }
  info <- unique(data.frame(
    node = trimws(as.character(info$node)),
    type = tolower(trimws(as.character(info$type))),
    stringsAsFactors = FALSE
  ))
  if (any(!nzchar(info$node)) || any(!info$type %in% c("protein", "drug"))) {
    stop("DRN info contains invalid node or type values: ", path, call. = FALSE)
  }
  info[order(match(info$type, c("protein", "drug")), info$node), , drop = FALSE]
}

.drug_level_counts <- function(path) {
  level <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("label", "count") %in% names(level))) {
    stop("Drug response level file must contain label and count columns: ", path, call. = FALSE)
  }
  labels <- tolower(trimws(as.character(level$label)))
  counts <- suppressWarnings(as.numeric(level$count))
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) {
    stop("Drug response level counts must be non-negative integers: ", path, call. = FALSE)
  }
  significant <- sum(counts[labels %in% c("sig", "significant")])
  list(significant = as.integer(significant), tested = as.integer(sum(counts)))
}

.drug_prediction_count <- function(path) {
  header <- utils::read.csv(
    path,
    nrows = 1L,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  max(0L, ncol(header) - 1L)
}

.drug_empty_hits <- function() {
  data.frame(
    module_uid = character(),
    drug_panel = character(),
    drug = character(),
    evidence_type = character(),
    effect_size = numeric(),
    effect_direction = character(),
    p_value = numeric(),
    p_adjust = numeric(),
    significant = logical(),
    tested_samples = integer(),
    drn_edge_number = integer(),
    stringsAsFactors = FALSE
  )
}

.drug_hits_from_drn <- function(drn, module_uid, panel) {
  if (!nrow(drn)) return(.drug_empty_hits())
  rows <- lapply(split(drn, drn$drug), function(group) {
    best <- group[order(group$p_adjust, group$p_value), , drop = FALSE][1L, ]
    data.frame(
      module_uid = module_uid,
      drug_panel = panel,
      drug = best$drug,
      evidence_type = "protein_drug_drn",
      effect_size = NA_real_,
      effect_direction = "not_available_in_legacy_output",
      p_value = best$p_value,
      p_adjust = best$p_adjust,
      significant = any(group$significant),
      tested_samples = NA_integer_,
      drn_edge_number = nrow(group),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result[order(result$p_adjust, result$p_value, result$drug), , drop = FALSE]
}

#' Standardize all-module drug-response results and drug-response networks
#'
#' Step 6 reads exactly the modules in `module_manifest.tsv`. Existing
#' subnetDR-style results can be supplied as named panel roots, or every input
#' file can be declared in an index table. The function writes canonical DRNs,
#' a drug-level evidence table, one summary row per module/panel and coverage
#' QC. No module is selected or ranked at this stage.
#'
#' @param module_manifest_file Explicit step-4 manifest path.
#' @param drug_response_path New step-6 output directory.
#' @param drug_response_roots Optional named paths to subnetDR-style output
#'   roots, for example `c(PRISM = "DrugResponse_PRISM")`.
#' @param drug_response_index_file Optional `.tsv`, `.csv` or `.xlsx` index
#'   with `module_uid`, `drug_panel`, `drn_file`, `drn_info_file`,
#'   `drug_level_file` and `prediction_file`. Paths may be relative to the
#'   index. Supply this or `drug_response_roots`, not both.
#' @param panels Optional panel filter.
#' @param network_methods Optional network filter.
#' @param module_methods Optional module-algorithm filter.
#' @param subtypes Optional subtype filter.
#' @param significance_cutoff Adjusted p-value cutoff for canonical DRN edges.
#' @param p_adjust_method `"none"` for legacy subnetDR DRNs or a method from
#'   [stats::p.adjust.methods].
#' @param strict Require complete module-by-panel input coverage and verify
#'   step-4 hashes. When `FALSE`, missing inputs produce explicit status rows.
#' @return The canonical `drug_response_summary.tsv` table with paths to all
#'   output contracts in attributes.
#' @export
drug_response_analysis <- function(
  module_manifest_file,
  drug_response_path,
  drug_response_roots = NULL,
  drug_response_index_file = NULL,
  panels = NULL,
  network_methods = NULL,
  module_methods = NULL,
  subtypes = NULL,
  significance_cutoff = 0.05,
  p_adjust_method = "none",
  strict = TRUE
) {
  if (!xor(is.null(drug_response_roots), is.null(drug_response_index_file))) {
    stop("Supply exactly one of drug_response_roots or drug_response_index_file.", call. = FALSE)
  }
  if (length(significance_cutoff) != 1L || !is.numeric(significance_cutoff) ||
      is.na(significance_cutoff) || significance_cutoff < 0 || significance_cutoff > 1) {
    stop("significance_cutoff must be between 0 and 1.", call. = FALSE)
  }
  allowed_adjustments <- c("none", stats::p.adjust.methods)
  if (length(p_adjust_method) != 1L || !p_adjust_method %in% allowed_adjustments) {
    stop("p_adjust_method must be 'none' or one of stats::p.adjust.methods.", call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("strict must be TRUE or FALSE.", call. = FALSE)
  }

  manifest <- read_module_manifest(
    module_manifest_file,
    check_files = TRUE,
    verify_hashes = strict
  )
  keep <- rep(TRUE, nrow(manifest))
  if (!is.null(network_methods)) {
    keep <- keep & manifest$network %in% normalize_network_name(network_methods, output = "display")
  }
  if (!is.null(module_methods)) {
    keep <- keep & manifest$method %in% normalize_method_name(module_methods, output = "display")
  }
  if (!is.null(subtypes)) keep <- keep & manifest$subtype %in% .normalize_subtype(subtypes)
  manifest <- manifest[keep, , drop = FALSE]
  if (!nrow(manifest)) stop("No manifest modules matched the requested filters.", call. = FALSE)

  input_index <- if (!is.null(drug_response_index_file)) {
    .drug_read_index(drug_response_index_file)
  } else {
    .drug_index_from_roots(manifest, drug_response_roots)
  }
  if (!is.null(panels)) {
    panels <- toupper(trimws(as.character(panels)))
    input_index <- input_index[input_index$drug_panel %in% panels, , drop = FALSE]
  } else {
    panels <- sort(unique(input_index$drug_panel))
  }
  if (!length(panels) || !nrow(input_index)) stop("No drug panels remained after filtering.", call. = FALSE)
  input_index <- input_index[input_index$module_uid %in% manifest$module_uid, , drop = FALSE]

  expected <- expand.grid(
    module_uid = manifest$module_uid,
    drug_panel = panels,
    stringsAsFactors = FALSE
  )
  expected$key <- paste(expected$module_uid, expected$drug_panel, sep = "\r")
  input_index$key <- paste(input_index$module_uid, input_index$drug_panel, sep = "\r")
  coverage_missing <- expected[!expected$key %in% input_index$key, c("module_uid", "drug_panel")]
  if (nrow(coverage_missing) && isTRUE(strict)) {
    stop(
      "Drug-response index does not cover every manifest module/panel; first missing: ",
      coverage_missing$module_uid[[1L]], "/", coverage_missing$drug_panel[[1L]],
      call. = FALSE
    )
  }
  if (nrow(coverage_missing)) {
    coverage_missing$drn_file <- NA_character_
    coverage_missing$drn_info_file <- NA_character_
    coverage_missing$drug_level_file <- NA_character_
    coverage_missing$prediction_file <- NA_character_
    coverage_missing$key <- paste(coverage_missing$module_uid, coverage_missing$drug_panel, sep = "\r")
    input_index <- rbind(input_index, coverage_missing[, names(input_index), drop = FALSE])
  }
  input_index <- input_index[match(expected$key, input_index$key), , drop = FALSE]
  input_index$key <- NULL
  missing_files <- .drug_missing_files(input_index)
  if (nrow(missing_files) && isTRUE(strict)) {
    stop(
      "Missing drug-response input file(s); first incomplete row: ",
      missing_files$module_uid[[1L]], "/", missing_files$drug_panel[[1L]],
      " [", missing_files$missing_fields[[1L]], "]",
      call. = FALSE
    )
  }

  output_dir <- .annotation_output_target(drug_response_path)
  staging_dir <- tempfile(
    pattern = paste0(".", basename(output_dir), "-"),
    tmpdir = dirname(output_dir)
  )
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create drug-response staging directory.", call. = FALSE)
  }
  on.exit({
    if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  summary_rows <- vector("list", nrow(input_index))
  hit_rows <- vector("list", nrow(input_index))
  source_rows <- input_index
  for (index in seq_len(nrow(input_index))) {
    input <- input_index[index, , drop = FALSE]
    module <- manifest[manifest$module_uid == input$module_uid, , drop = FALSE]
    missing <- c("drn_file", "drn_info_file", "drug_level_file", "prediction_file")[
      !vapply(input[, c("drn_file", "drn_info_file", "drug_level_file", "prediction_file"), drop = FALSE],
              function(path) !is.na(path) && file.exists(path), logical(1))
    ]
    relative_dir <- file.path("modules", input$module_uid, input$drug_panel)
    drn_relative <- file.path(relative_dir, "drn.tsv")
    info_relative <- file.path(relative_dir, "drn_info.tsv")
    if (length(missing)) {
      summary_rows[[index]] <- data.frame(
        module_uid = module$module_uid,
        network = module$network,
        method = module$method,
        subtype = module$subtype,
        module = module$module,
        module_size = module$module_size,
        drug_panel = input$drug_panel,
        drug_number = NA_integer_,
        tested_drug_number = NA_integer_,
        drug_response_density = NA_real_,
        drn_drug_number = NA_integer_,
        drn_edge_number = NA_integer_,
        p_adjust_method = p_adjust_method,
        significance_cutoff = significance_cutoff,
        drn_file = NA_character_,
        drn_info_file = NA_character_,
        analysis_status = "missing_input",
        status_reason = paste("Missing", paste(missing, collapse = ", ")),
        stringsAsFactors = FALSE
      )
      hit_rows[[index]] <- .drug_empty_hits()
      next
    }

    drn <- .drug_read_drn(input$drn_file, p_adjust_method, significance_cutoff)
    info <- .drug_read_info(input$drn_info_file)
    level <- .drug_level_counts(input$drug_level_file)
    prediction_count <- .drug_prediction_count(input$prediction_file)
    if (level$tested != prediction_count) {
      stop(
        "Tested-drug count differs between drug level and predictions for ",
        input$module_uid, "/", input$drug_panel, ": ", level$tested, " vs ", prediction_count,
        call. = FALSE
      )
    }
    module_output <- file.path(staging_dir, relative_dir)
    if (!dir.create(module_output, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create canonical DRN directory: ", module_output, call. = FALSE)
    }
    .manifest_write_tsv(drn, file.path(staging_dir, drn_relative))
    .manifest_write_tsv(info, file.path(staging_dir, info_relative))
    hit_rows[[index]] <- .drug_hits_from_drn(drn, input$module_uid, input$drug_panel)
    summary_rows[[index]] <- data.frame(
      module_uid = module$module_uid,
      network = module$network,
      method = module$method,
      subtype = module$subtype,
      module = module$module,
      module_size = module$module_size,
      drug_panel = input$drug_panel,
      drug_number = level$significant,
      tested_drug_number = level$tested,
      drug_response_density = level$significant / module$module_size,
      drn_drug_number = length(unique(drn$drug[drn$significant])),
      drn_edge_number = sum(drn$significant),
      p_adjust_method = p_adjust_method,
      significance_cutoff = significance_cutoff,
      drn_file = gsub("\\\\", "/", drn_relative),
      drn_info_file = gsub("\\\\", "/", info_relative),
      analysis_status = "completed",
      status_reason = "",
      stringsAsFactors = FALSE
    )
  }

  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL
  hits_nonempty <- hit_rows[vapply(hit_rows, nrow, integer(1)) > 0L]
  hits <- if (length(hits_nonempty)) do.call(rbind, hits_nonempty) else .drug_empty_hits()
  rownames(hits) <- NULL
  source_rows$key <- NULL
  source_rows$input_status <- ifelse(
    paste(source_rows$module_uid, source_rows$drug_panel, sep = "\r") %in%
      paste(missing_files$module_uid %||% character(), missing_files$drug_panel %||% character(), sep = "\r"),
    "missing_input",
    "available"
  )
  coverage <- stats::aggregate(
    summary$analysis_status == "completed",
    by = list(drug_panel = summary$drug_panel),
    FUN = sum
  )
  names(coverage)[[2L]] <- "completed_module_number"
  coverage$expected_module_number <- nrow(manifest)
  coverage$coverage_fraction <- coverage$completed_module_number / coverage$expected_module_number
  coverage$coverage_pass <- coverage$completed_module_number == coverage$expected_module_number

  summary_file <- file.path(staging_dir, "drug_response_summary.tsv")
  hits_file <- file.path(staging_dir, "drug_response_hits.tsv")
  source_file <- file.path(staging_dir, "drug_response_source_index.tsv")
  coverage_file <- file.path(staging_dir, "drug_response_coverage.tsv")
  .manifest_write_tsv(summary, summary_file)
  .manifest_write_tsv(hits, hits_file)
  .manifest_write_tsv(source_rows, source_file)
  .manifest_write_tsv(coverage, coverage_file)

  if (!file.rename(staging_dir, output_dir)) {
    stop("Could not finalize drug-response output directory: ", output_dir, call. = FALSE)
  }
  attr(summary, "output_dir") <- output_dir
  attr(summary, "summary_file") <- file.path(output_dir, basename(summary_file))
  attr(summary, "hits_file") <- file.path(output_dir, basename(hits_file))
  attr(summary, "source_index_file") <- file.path(output_dir, basename(source_file))
  attr(summary, "coverage_file") <- file.path(output_dir, basename(coverage_file))
  summary
}
