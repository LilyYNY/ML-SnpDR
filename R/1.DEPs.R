# Step 1: differential expression -------------------------------------------

#' Run subtype-versus-rest differential expression analysis
#'
#' This chaining-compatible implementation preserves subnetDR's Wilcoxon
#' rank-sum test, BH-adjusted P value, detection-rate filter and fold-change
#' labels. It additionally writes one canonical long TSV used directly by
#' step 2, a significant-only table, a summary and an optional subnetDR-style
#' Excel workbook.
#'
#' @param input_dir Optional compatibility root containing `expr_name` and
#'   `pheno_name`. Explicit files take precedence.
#' @param expr_name,pheno_name File names relative to `input_dir`.
#' @param expression_file,phenotype_file Explicit expression and phenotype
#'   files (`.tsv`, `.csv`, `.txt` or `.xlsx`).
#' @param output_dir New step-1 output directory.
#' @param gene_column,sample_column,subtype_column Input identity columns.
#' @param detection_threshold Minimum non-missing, non-zero fraction required
#'   in both the target subtype and all other samples.
#' @param fc_threshold Fold-change threshold; down-regulation uses its inverse.
#' @param p_threshold Adjusted-P threshold used for labels.
#' @param p_adjust_method Method accepted by [stats::p.adjust()].
#' @param pseudocount Pseudocount used only for the reported log2 fold change.
#' @param subtypes Optional subtype subset. Defaults to all phenotype subtypes.
#' @param write_legacy_excel Write `diff_expression_results_all.xlsx` with the
#'   original subnetDR sheet layout.
#' @param strict Require complete sample matching and at least two samples in
#'   each target/rest group.
#' @return Canonical differential-expression table with output paths stored in
#'   attributes.
#' @export
run_diff_expr_analysis <- function(
  input_dir = NULL,
  expr_name = "expression.xlsx",
  pheno_name = "subtype.xlsx",
  expression_file = NULL,
  phenotype_file = NULL,
  output_dir = NULL,
  gene_column = NULL,
  sample_column = "Sample",
  subtype_column = "Subtype",
  detection_threshold = 0.75,
  fc_threshold = 1.5,
  p_threshold = 0.05,
  p_adjust_method = "BH",
  pseudocount = 1,
  subtypes = NULL,
  write_legacy_excel = TRUE,
  strict = TRUE
) {
  if (is.null(expression_file)) {
    if (is.null(input_dir)) stop("Supply expression_file or input_dir.", call. = FALSE)
    expression_file <- file.path(input_dir, expr_name)
  }
  if (is.null(phenotype_file)) {
    if (is.null(input_dir)) stop("Supply phenotype_file or input_dir.", call. = FALSE)
    phenotype_file <- file.path(input_dir, pheno_name)
  }
  if (is.null(output_dir)) {
    if (is.null(input_dir)) stop("Supply output_dir when input_dir is not used.", call. = FALSE)
    output_dir <- file.path(input_dir, "DEPS_result")
  }
  probabilities <- c(detection_threshold, p_threshold)
  if (any(!is.finite(probabilities)) || any(probabilities < 0) || any(probabilities > 1)) {
    stop("detection_threshold and p_threshold must be between 0 and 1.", call. = FALSE)
  }
  if (length(fc_threshold) != 1L || !is.finite(fc_threshold) || fc_threshold <= 1) {
    stop("fc_threshold must be one finite value greater than 1.", call. = FALSE)
  }
  if (!p_adjust_method %in% stats::p.adjust.methods) {
    stop("Unsupported p_adjust_method.", call. = FALSE)
  }
  if (length(pseudocount) != 1L || !is.finite(pseudocount) || pseudocount < 0) {
    stop("pseudocount must be one non-negative finite value.", call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict) ||
      !is.logical(write_legacy_excel) || length(write_legacy_excel) != 1L || is.na(write_legacy_excel)) {
    stop("strict and write_legacy_excel must be TRUE or FALSE.", call. = FALSE)
  }

  expression <- .upstream_read_tabular(expression_file)
  phenotype <- .upstream_read_tabular(phenotype_file)
  if (!nrow(expression) || ncol(expression) < 3L) {
    stop("Expression input must contain one gene column and at least two sample columns.", call. = FALSE)
  }
  gene_column <- if (is.null(gene_column)) names(expression)[[1L]] else
    .upstream_column(expression, gene_column, c("gene", "protein", "symbol", "gene_symbol"), "Gene")
  if (!gene_column %in% names(expression)) stop("Expression gene column was not found: ", gene_column, call. = FALSE)
  sample_column <- .upstream_column(phenotype, sample_column, c("sample", "sample_id", "sampleid"), "Sample")
  subtype_column <- .upstream_column(phenotype, subtype_column, c("subtype", "class", "group"), "Subtype")

  genes <- trimws(as.character(expression[[gene_column]]))
  if (any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("Expression gene identifiers must be non-empty and unique.", call. = FALSE)
  }
  sample_data <- expression[, setdiff(names(expression), gene_column), drop = FALSE]
  numeric_data <- lapply(sample_data, function(x) suppressWarnings(as.numeric(as.character(x))))
  numeric_data <- as.data.frame(numeric_data, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(numeric_data) <- genes
  nonmissing_original <- !is.na(sample_data)
  conversion_failed <- nonmissing_original & is.na(numeric_data)
  if (any(conversion_failed)) stop("Expression matrix contains non-numeric values.", call. = FALSE)

  phenotype[[sample_column]] <- trimws(as.character(phenotype[[sample_column]]))
  phenotype[[subtype_column]] <- .normalize_subtype(phenotype[[subtype_column]])
  if (any(!nzchar(phenotype[[sample_column]])) || anyDuplicated(phenotype[[sample_column]])) {
    stop("Phenotype sample identifiers must be non-empty and unique.", call. = FALSE)
  }
  missing_samples <- setdiff(phenotype[[sample_column]], names(numeric_data))
  if (length(missing_samples) && isTRUE(strict)) {
    stop("Phenotype samples missing from expression matrix: ", paste(missing_samples, collapse = ", "), call. = FALSE)
  }
  phenotype <- phenotype[phenotype[[sample_column]] %in% names(numeric_data), , drop = FALSE]
  if (!nrow(phenotype)) stop("No phenotype samples matched the expression matrix.", call. = FALSE)
  if (is.null(subtypes)) subtypes <- unique(phenotype[[subtype_column]])
  subtypes <- sort(unique(.normalize_subtype(subtypes)))
  absent <- setdiff(subtypes, phenotype[[subtype_column]])
  if (length(absent)) stop("Requested subtype(s) absent from phenotype: ", paste(absent, collapse = ", "), call. = FALSE)

  result_rows <- vector("list", length(subtypes))
  for (index in seq_along(subtypes)) {
    subtype <- subtypes[[index]]
    target_samples <- phenotype[[sample_column]][phenotype[[subtype_column]] == subtype]
    other_samples <- phenotype[[sample_column]][phenotype[[subtype_column]] != subtype]
    if (isTRUE(strict) && (length(target_samples) < 2L || length(other_samples) < 2L)) {
      stop("Subtype ", subtype, " requires at least two target and two rest samples.", call. = FALSE)
    }
    target <- as.matrix(numeric_data[, target_samples, drop = FALSE])
    other <- as.matrix(numeric_data[, other_samples, drop = FALSE])
    detection_target <- rowMeans(!is.na(target) & target != 0)
    detection_other <- rowMeans(!is.na(other) & other != 0)
    testable <- detection_target >= detection_threshold & detection_other >= detection_threshold
    mean_target <- rowMeans(target, na.rm = TRUE)
    mean_other <- rowMeans(other, na.rm = TRUE)
    mean_target[!testable] <- NA_real_
    mean_other[!testable] <- NA_real_
    fold_change <- ifelse(testable & mean_other != 0, mean_target / mean_other, NA_real_)
    log_numerator <- mean_target + pseudocount
    log_denominator <- mean_other + pseudocount
    log2_fold_change <- ifelse(
      testable & log_numerator > 0 & log_denominator > 0,
      log2(log_numerator / log_denominator),
      NA_real_
    )
    p_value <- rep(NA_real_, length(genes))
    for (gene_index in which(testable)) {
      x <- target[gene_index, ]
      y <- other[gene_index, ]
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
      if (length(x) && length(y)) {
        p_value[[gene_index]] <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value)
      }
    }
    p_adjust <- stats::p.adjust(p_value, method = p_adjust_method)
    label <- rep("non_significant", length(genes))
    label[testable & is.finite(p_adjust) & p_adjust < p_threshold &
            is.finite(fold_change) & fold_change > fc_threshold] <- "up"
    label[testable & is.finite(p_adjust) & p_adjust < p_threshold &
            is.finite(fold_change) & fold_change < 1 / fc_threshold] <- "down"
    result_rows[[index]] <- data.frame(
      subtype = subtype,
      gene = genes,
      mean_target = mean_target,
      mean_other = mean_other,
      fold_change = fold_change,
      log2_fold_change = log2_fold_change,
      p_value = p_value,
      p_adjust = p_adjust,
      detection_target = detection_target,
      detection_other = detection_other,
      label = label,
      n_target = length(target_samples),
      n_other = length(other_samples),
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, result_rows)
  result <- result[order(result$subtype, result$p_adjust, -abs(result$log2_fold_change), result$gene, na.last = TRUE), , drop = FALSE]
  rownames(result) <- NULL
  significant <- result[result$label %in% c("up", "down"), , drop = FALSE]
  summary_rows <- lapply(split(result, result$subtype), function(group) {
    counts <- table(factor(group$label, levels = c("up", "down", "non_significant")))
    data.frame(
      subtype = unique(group$subtype),
      label = names(counts),
      gene_number = as.integer(counts),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL

  output_target <- .upstream_output_target(output_dir, "Differential-expression")
  staging_dir <- tempfile(pattern = paste0(".", basename(output_target), "-"), tmpdir = dirname(output_target))
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create differential-expression staging directory.", call. = FALSE)
  }
  on.exit(if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
  result_file <- file.path(staging_dir, "differential_expression.tsv")
  significant_file <- file.path(staging_dir, "differential_expression_significant.tsv")
  summary_file <- file.path(staging_dir, "differential_expression_summary.tsv")
  .manifest_write_tsv(result, result_file)
  .manifest_write_tsv(significant, significant_file)
  .manifest_write_tsv(summary, summary_file)

  legacy_file <- NULL
  if (isTRUE(write_legacy_excel)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("write_legacy_excel=TRUE requires the openxlsx package.", call. = FALSE)
    }
    legacy_file <- file.path(staging_dir, "diff_expression_results_all.xlsx")
    workbook <- openxlsx::createWorkbook()
    for (subtype in subtypes) {
      group <- result[result$subtype == subtype, , drop = FALSE]
      legacy <- data.frame(
        Gene = group$gene,
        mean_group1 = group$mean_target,
        mean_group2 = group$mean_other,
        fc = group$fold_change,
        log2FC = group$log2_fold_change,
        p_value = group$p_value,
        adj_p = group$p_adjust,
        detection_group1 = group$detection_target,
        detection_group2 = group$detection_other,
        label = ifelse(group$label == "non_significant", "non-significant", group$label),
        stringsAsFactors = FALSE
      )
      stats <- summary[summary$subtype == subtype, c("label", "gene_number"), drop = FALSE]
      names(stats) <- c("Label", "Count")
      openxlsx::addWorksheet(workbook, paste0(subtype, "_DiffResults"))
      openxlsx::writeData(workbook, paste0(subtype, "_DiffResults"), legacy)
      openxlsx::addWorksheet(workbook, paste0(subtype, "_Statistics"))
      openxlsx::writeData(workbook, paste0(subtype, "_Statistics"), stats)
    }
    openxlsx::saveWorkbook(workbook, legacy_file, overwrite = TRUE)
  }
  if (!file.rename(staging_dir, output_target)) {
    stop("Could not finalize differential-expression output directory.", call. = FALSE)
  }
  attr(result, "output_dir") <- output_target
  attr(result, "result_file") <- file.path(output_target, basename(result_file))
  attr(result, "significant_file") <- file.path(output_target, basename(significant_file))
  attr(result, "summary_file") <- file.path(output_target, basename(summary_file))
  attr(result, "legacy_excel_file") <- if (is.null(legacy_file)) NULL else file.path(output_target, basename(legacy_file))
  result
}

