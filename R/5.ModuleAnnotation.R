# Module annotation ----------------------------------------------------------

.annotation_manifest_path <- function(module_manifest_file, base_input_path) {
  if (!is.null(module_manifest_file)) {
    candidates <- module_manifest_file
  } else if (length(base_input_path) == 1L && file.exists(base_input_path) &&
             !dir.exists(base_input_path)) {
    candidates <- base_input_path
  } else {
    candidates <- c(
      file.path(base_input_path, "standardized", "module_manifest.tsv"),
      file.path(base_input_path, "module_manifest.tsv")
    )
  }
  candidates <- candidates[file.exists(candidates) & !dir.exists(candidates)]
  if (!length(candidates)) {
    stop(
      "Could not find module_manifest.tsv. Supply module_manifest_file or a ",
      "ModuleSelection directory containing standardized/module_manifest.tsv.",
      call. = FALSE
    )
  }
  normalizePath(candidates[[1L]], winslash = "/", mustWork = TRUE)
}

.annotation_output_target <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("base_output_path must be one non-empty path.", call. = FALSE)
  }
  path <- path.expand(path)
  parent <- dirname(path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output parent directory: ", parent, call. = FALSE)
  }
  parent <- normalizePath(parent, winslash = "/", mustWork = TRUE)
  target <- file.path(parent, basename(path))
  if (file.exists(target)) {
    stop(
      "base_output_path already exists; choose a new directory to avoid overwriting: ",
      target,
      call. = FALSE
    )
  }
  target
}

.annotation_read_tabular <- function(path) {
  if (length(path) != 1L || !file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }
  extension <- tolower(tools::file_ext(path))
  if (extension == "csv") {
    return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (extension %in% c("tsv", "txt")) {
    return(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (extension %in% c("xlsx", "xlsm")) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required to read Excel input files.", call. = FALSE)
    }
    return(openxlsx::read.xlsx(path))
  }
  stop("Supported tabular formats are .csv, .tsv, .txt, .xlsx and .xlsm.", call. = FALSE)
}

.annotation_read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(trimws(lines))]
  rows <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(fields) < 3L) {
      stop("Every GMT row must contain a term, description and at least one gene.", call. = FALSE)
    }
    data.frame(
      database = "CUSTOM",
      category = "CUSTOM",
      term_id = fields[[1L]],
      description = fields[[2L]],
      gene = fields[-c(1L, 2L)],
      stringsAsFactors = FALSE
    )
  })
  if (!length(rows)) stop("GMT file contains no gene sets: ", path, call. = FALSE)
  do.call(rbind, rows)
}

.annotation_standardize_gene_sets <- function(gene_sets, source = "gene_sets") {
  if (!is.data.frame(gene_sets)) {
    stop(source, " must be a data frame.", call. = FALSE)
  }
  aliases <- list(
    database = c("database", "db", "collection"),
    category = c("category", "subcategory", "subcollection"),
    term_id = c("term_id", "term", "id", "gs_id", "gs_name"),
    description = c("description", "term_name", "gs_name"),
    gene = c("gene", "gene_symbol", "symbol")
  )
  lowered <- tolower(names(gene_sets))
  selected <- vapply(aliases, function(options) {
    hit <- match(options, lowered, nomatch = 0L)
    hit <- hit[hit > 0L]
    if (length(hit)) names(gene_sets)[hit[[1L]]] else NA_character_
  }, character(1))
  if (is.na(selected[["term_id"]]) || is.na(selected[["gene"]])) {
    stop(
      source, " must contain term_id (or term/gs_name) and gene (or gene_symbol) columns.",
      call. = FALSE
    )
  }

  result <- data.frame(
    database = if (is.na(selected[["database"]])) "CUSTOM" else gene_sets[[selected[["database"]]]],
    category = if (is.na(selected[["category"]])) "CUSTOM" else gene_sets[[selected[["category"]]]],
    term_id = gene_sets[[selected[["term_id"]]]],
    description = if (is.na(selected[["description"]])) {
      gene_sets[[selected[["term_id"]]]]
    } else {
      gene_sets[[selected[["description"]]]]
    },
    gene = gene_sets[[selected[["gene"]]]],
    stringsAsFactors = FALSE
  )
  result[] <- lapply(result, function(x) trimws(as.character(x)))
  result$database <- toupper(result$database)
  result$category[!nzchar(result$category)] <- result$database[!nzchar(result$category)]
  result$description[!nzchar(result$description)] <- result$term_id[!nzchar(result$description)]
  keep <- stats::complete.cases(result[, c("database", "term_id", "gene")]) &
    nzchar(result$database) & nzchar(result$term_id) & nzchar(result$gene)
  result <- unique(result[keep, , drop = FALSE])
  if (!nrow(result)) stop(source, " contains no usable gene-set memberships.", call. = FALSE)
  rownames(result) <- NULL
  result
}

.annotation_read_gene_sets <- function(gene_set_file) {
  extension <- tolower(tools::file_ext(gene_set_file))
  raw <- if (extension == "gmt") {
    .annotation_read_gmt(gene_set_file)
  } else {
    .annotation_read_tabular(gene_set_file)
  }
  .annotation_standardize_gene_sets(raw, source = gene_set_file)
}

.annotation_msigdb <- function(databases, species) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop(
      "Package 'msigdbr' is required when gene_sets/gene_set_file is not supplied.",
      call. = FALSE
    )
  }
  supported <- c("GO_BP", "KEGG", "HALLMARK")
  unknown <- setdiff(databases, supported)
  if (length(unknown)) {
    stop(
      "Built-in MSigDB loading supports GO_BP, KEGG and HALLMARK. ",
      "Supply gene_sets for: ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  load_one <- function(database) {
    definition <- switch(
      database,
      GO_BP = list(collection = "C5", subcollection = "GO:BP"),
      KEGG = list(collection = "C2", subcollection = "CP:KEGG_LEGACY"),
      HALLMARK = list(collection = "H", subcollection = NULL)
    )
    formal_names <- names(formals(msigdbr::msigdbr))
    args <- list(species = species)
    if ("collection" %in% formal_names) {
      args$collection <- definition$collection
      if (!is.null(definition$subcollection)) args$subcollection <- definition$subcollection
    } else {
      args$category <- definition$collection
      if (!is.null(definition$subcollection)) args$subcategory <- definition$subcollection
    }
    raw <- do.call(msigdbr::msigdbr, args)
    data.frame(
      database = database,
      category = if ("gs_subcollection" %in% names(raw)) raw$gs_subcollection else database,
      term_id = if ("gs_id" %in% names(raw)) raw$gs_id else raw$gs_name,
      description = raw$gs_name,
      gene = raw$gene_symbol,
      stringsAsFactors = FALSE
    )
  }
  .annotation_standardize_gene_sets(
    do.call(rbind, lapply(databases, load_one)),
    source = "msigdbr"
  )
}

.annotation_read_genes <- function(path) {
  data <- .annotation_read_tabular(path)
  if (!ncol(data)) stop("Background gene file contains no columns: ", path, call. = FALSE)
  gene_column <- names(data)[tolower(names(data)) %in% c("gene", "gene_symbol", "symbol")]
  column <- if (length(gene_column)) gene_column[[1L]] else names(data)[[1L]]
  unique(trimws(as.character(data[[column]])))
}

.annotation_empty_result <- function() {
  data.frame(
    module_uid = character(),
    network = character(),
    method = character(),
    subtype = character(),
    module = character(),
    database = character(),
    category = character(),
    term_id = character(),
    description = character(),
    gene_ratio = character(),
    background_ratio = character(),
    rich_factor = numeric(),
    fold_enrichment = numeric(),
    p_value = numeric(),
    p_adjust = numeric(),
    q_value = numeric(),
    gene_count = integer(),
    gene_ids = character(),
    analysis_status = character(),
    status_reason = character(),
    stringsAsFactors = FALSE
  )
}

.annotation_status_row <- function(manifest_row, status, reason) {
  row <- .annotation_empty_result()[NA_integer_, , drop = FALSE]
  row[1L, ] <- NA
  row$module_uid <- manifest_row$module_uid
  row$network <- manifest_row$network
  row$method <- manifest_row$method
  row$subtype <- manifest_row$subtype
  row$module <- manifest_row$module
  row$database <- "NONE"
  row$category <- "NONE"
  row$term_id <- "NONE"
  row$description <- "No significant enrichment term"
  row$gene_count <- 0L
  row$gene_ids <- ""
  row$analysis_status <- status
  row$status_reason <- reason
  row
}

.annotation_enrich_module <- function(
  genes,
  manifest_row,
  gene_sets,
  universe,
  pvalue_cutoff,
  p_adjust_method,
  p_adjust_cutoff,
  qvalue_cutoff,
  min_set_size,
  max_set_size
) {
  genes <- intersect(unique(genes), universe)
  if (!length(genes)) {
    return(.annotation_status_row(
      manifest_row,
      "completed_no_testable_genes",
      "No module genes overlapped the annotation universe."
    ))
  }

  split_key <- paste(gene_sets$database, gene_sets$term_id, sep = "\r")
  terms <- split(gene_sets, split_key)
  terms <- terms[vapply(terms, function(x) {
    size <- length(unique(x$gene[x$gene %in% universe]))
    size >= min_set_size && size <= max_set_size
  }, logical(1))]
  if (!length(terms)) {
    return(.annotation_status_row(
      manifest_row,
      "completed_no_testable_terms",
      "No gene set passed the configured size limits."
    ))
  }

  total <- length(universe)
  query_size <- length(genes)
  rows <- lapply(terms, function(term) {
    term_genes <- unique(term$gene[term$gene %in% universe])
    overlap <- sort(intersect(genes, term_genes))
    overlap_size <- length(overlap)
    set_size <- length(term_genes)
    p_value <- if (overlap_size) {
      stats::phyper(overlap_size - 1L, set_size, total - set_size, query_size, lower.tail = FALSE)
    } else {
      1
    }
    rich_factor <- if (set_size) overlap_size / set_size else NA_real_
    fold_enrichment <- if (query_size && set_size) {
      (overlap_size / query_size) / (set_size / total)
    } else {
      NA_real_
    }
    data.frame(
      module_uid = manifest_row$module_uid,
      network = manifest_row$network,
      method = manifest_row$method,
      subtype = manifest_row$subtype,
      module = manifest_row$module,
      database = term$database[[1L]],
      category = term$category[[1L]],
      term_id = term$term_id[[1L]],
      description = term$description[[1L]],
      gene_ratio = paste0(overlap_size, "/", query_size),
      background_ratio = paste0(set_size, "/", total),
      rich_factor = rich_factor,
      fold_enrichment = fold_enrichment,
      p_value = p_value,
      p_adjust = NA_real_,
      q_value = NA_real_,
      gene_count = overlap_size,
      gene_ids = paste(overlap, collapse = "/"),
      analysis_status = "completed",
      status_reason = "",
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$p_adjust <- stats::ave(
    result$p_value,
    result$database,
    FUN = function(x) stats::p.adjust(x, method = p_adjust_method)
  )
  result$q_value <- result$p_adjust
  keep <- result$p_value <= pvalue_cutoff &
    result$p_adjust <= p_adjust_cutoff &
    result$q_value <= qvalue_cutoff &
    result$gene_count > 0L
  result <- result[keep, , drop = FALSE]
  if (!nrow(result)) {
    return(.annotation_status_row(
      manifest_row,
      "completed_no_significant_terms",
      "No term passed p-value, adjusted p-value and q-value thresholds."
    ))
  }
  result <- result[order(result$database, result$p_adjust, result$p_value, result$term_id), , drop = FALSE]
  rownames(result) <- NULL
  result
}

#' Annotate every module in a canonical module manifest
#'
#' This chaining-compatible replacement for subnetDR's
#' `functional_annotation()` reads the exact modules emitted by step 4 instead
#' of rediscovering aggregate node files. It writes one unified table, a top-N
#' table, module-level files and QC. A custom gene-set table/GMT can be supplied
#' for reproducible or offline runs; otherwise GO biological process, KEGG and
#' Hallmark sets are loaded with `msigdbr`.
#'
#' @param subtype_file Optional phenotype file used only to filter subtypes.
#' @param base_input_path A `module_manifest.tsv` file or ModuleSelection root.
#' @param base_output_path New annotation output directory.
#' @param network_method Networks retained from the manifest.
#' @param module_method Module algorithms retained from the manifest.
#' @param module_manifest_file Optional explicit manifest path; takes priority
#'   over `base_input_path`.
#' @param gene_set_file Optional `.tsv`, `.csv`, `.xlsx` or `.gmt` gene-set
#'   file. Tabular files require `term_id` and `gene`; `database`, `category`
#'   and `description` are optional.
#' @param gene_sets Optional in-memory gene-set data frame with the same
#'   columns as `gene_set_file`.
#' @param background_gene_file Optional tabular file containing background
#'   genes in `gene`, `gene_symbol`, `symbol` or the first column.
#' @param background_genes Optional character vector of background genes.
#' @param databases Optional database names to retain. When gene sets are not
#'   supplied, defaults to `GO_BP`, `KEGG` and `HALLMARK`.
#' @param species Species passed to `msigdbr` when built-in sets are used.
#' @param subtypes Optional explicit subtype filter.
#' @param subtype_column Subtype column in `subtype_file`.
#' @param pvalueCutoff Raw p-value threshold.
#' @param pAdjustMethod Method passed to [stats::p.adjust()].
#' @param pAdjustCutoff Adjusted p-value threshold.
#' @param qvalueCutoff Q-value threshold; BH adjusted p-values are reported as
#'   q-values in this over-representation implementation.
#' @param minGSSize,maxGSSize Inclusive gene-set size limits.
#' @param top_n Number of terms retained per module and database in the top
#'   table.
#' @param strict Whether to verify module file checksums.
#' @param write_module_files Whether to write one result file per module.
#' @return The unified annotation table with output file paths in attributes.
#' @export
functional_annotation <- function(
  subtype_file = NULL,
  base_input_path = "./ModuleSelection/",
  base_output_path = "./FunctionalAnnotation/",
  network_method = c("String", "physicalPPIN", "chengF"),
  module_method = c("Louvain", "WF"),
  module_manifest_file = NULL,
  gene_set_file = NULL,
  gene_sets = NULL,
  background_gene_file = NULL,
  background_genes = NULL,
  databases = NULL,
  species = "Homo sapiens",
  subtypes = NULL,
  subtype_column = "Subtype",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  pAdjustCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  top_n = 15,
  strict = TRUE,
  write_module_files = TRUE
) {
  exclusive_inputs <- sum(!vapply(list(gene_set_file, gene_sets), is.null, logical(1)))
  if (exclusive_inputs > 1L) stop("Supply gene_set_file or gene_sets, not both.", call. = FALSE)
  background_inputs <- sum(!vapply(list(background_gene_file, background_genes), is.null, logical(1)))
  if (background_inputs > 1L) {
    stop("Supply background_gene_file or background_genes, not both.", call. = FALSE)
  }
  numeric_parameters <- list(
    pvalueCutoff = pvalueCutoff,
    pAdjustCutoff = pAdjustCutoff,
    qvalueCutoff = qvalueCutoff
  )
  invalid_probability <- vapply(numeric_parameters, function(x) {
    length(x) != 1L || !is.numeric(x) || is.na(x) || x < 0 || x > 1
  }, logical(1))
  if (any(invalid_probability)) {
    stop(
      paste(names(numeric_parameters)[invalid_probability], collapse = ", "),
      " must be numeric values between 0 and 1.",
      call. = FALSE
    )
  }
  if (!pAdjustMethod %in% stats::p.adjust.methods) {
    stop("pAdjustMethod must be one of stats::p.adjust.methods.", call. = FALSE)
  }
  integer_parameters <- c(minGSSize = minGSSize, maxGSSize = maxGSSize, top_n = top_n)
  if (any(!is.finite(integer_parameters)) || any(integer_parameters < 1) ||
      any(integer_parameters != floor(integer_parameters))) {
    stop("minGSSize, maxGSSize and top_n must be positive integers.", call. = FALSE)
  }
  if (minGSSize > maxGSSize) stop("minGSSize cannot exceed maxGSSize.", call. = FALSE)
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict) ||
      !is.logical(write_module_files) || length(write_module_files) != 1L ||
      is.na(write_module_files)) {
    stop("strict and write_module_files must be TRUE or FALSE.", call. = FALSE)
  }

  manifest_file <- .annotation_manifest_path(module_manifest_file, base_input_path)
  manifest <- read_module_manifest(
    manifest_file,
    check_files = TRUE,
    verify_hashes = strict
  )
  networks <- unique(normalize_network_name(network_method, output = "display"))
  methods <- unique(normalize_method_name(module_method, output = "display"))
  if (!is.null(subtype_file) && !is.null(subtypes)) {
    stop("Supply subtype_file or subtypes, not both.", call. = FALSE)
  }
  if (!is.null(subtype_file)) {
    subtypes <- .selection_read_subtypes(subtype_file, subtype_column)
  } else if (!is.null(subtypes)) {
    subtypes <- unique(.normalize_subtype(subtypes))
  }
  keep <- manifest$network %in% networks & manifest$method %in% methods
  if (!is.null(subtypes)) keep <- keep & manifest$subtype %in% subtypes
  manifest <- manifest[keep, , drop = FALSE]
  if (!nrow(manifest)) stop("No manifest modules matched the requested filters.", call. = FALSE)

  if (!is.null(gene_set_file)) {
    gene_sets <- .annotation_read_gene_sets(gene_set_file)
  } else if (!is.null(gene_sets)) {
    gene_sets <- .annotation_standardize_gene_sets(gene_sets)
  } else {
    databases <- databases %||% c("GO_BP", "KEGG", "HALLMARK")
    databases <- toupper(databases)
    gene_sets <- .annotation_msigdb(databases, species)
  }
  if (!is.null(databases)) {
    databases <- toupper(trimws(as.character(databases)))
    gene_sets <- gene_sets[gene_sets$database %in% databases, , drop = FALSE]
    if (!nrow(gene_sets)) {
      stop("No gene sets matched databases: ", paste(databases, collapse = ", "), call. = FALSE)
    }
  }

  if (!is.null(background_gene_file)) {
    background_genes <- .annotation_read_genes(background_gene_file)
  }
  set_genes <- unique(gene_sets$gene)
  universe <- if (is.null(background_genes)) {
    set_genes
  } else {
    unique(trimws(as.character(background_genes)))
  }
  universe <- sort(intersect(universe[nzchar(universe)], set_genes))
  if (!length(universe)) stop("Annotation background has no overlap with gene sets.", call. = FALSE)

  output_dir <- .annotation_output_target(base_output_path)
  staging_dir <- tempfile(
    pattern = paste0(".", basename(output_dir), "-"),
    tmpdir = dirname(output_dir)
  )
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create annotation staging directory.", call. = FALSE)
  }
  on.exit({
    if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  result_rows <- vector("list", nrow(manifest))
  qc_rows <- vector("list", nrow(manifest))
  for (index in seq_len(nrow(manifest))) {
    manifest_row <- manifest[index, , drop = FALSE]
    nodes <- utils::read.delim(
      manifest_row$node_file_abs,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    )
    if (!all(c("node", "module") %in% names(nodes))) {
      stop("Module node file must contain node and module columns: ", manifest_row$node_file_abs, call. = FALSE)
    }
    module_genes <- unique(trimws(as.character(nodes$node)))
    module_result <- .annotation_enrich_module(
      genes = module_genes,
      manifest_row = manifest_row,
      gene_sets = gene_sets,
      universe = universe,
      pvalue_cutoff = pvalueCutoff,
      p_adjust_method = pAdjustMethod,
      p_adjust_cutoff = pAdjustCutoff,
      qvalue_cutoff = qvalueCutoff,
      min_set_size = as.integer(minGSSize),
      max_set_size = as.integer(maxGSSize)
    )
    result_rows[[index]] <- module_result
    significant <- module_result$analysis_status == "completed"
    qc_rows[[index]] <- data.frame(
      module_uid = manifest_row$module_uid,
      module_gene_number = length(module_genes),
      tested_gene_number = length(intersect(module_genes, universe)),
      tested_term_number = length(unique(paste(gene_sets$database, gene_sets$term_id))),
      significant_term_number = sum(significant),
      analysis_status = module_result$analysis_status[[1L]],
      status_reason = module_result$status_reason[[1L]],
      stringsAsFactors = FALSE
    )
    if (isTRUE(write_module_files)) {
      module_dir <- file.path(staging_dir, "modules", manifest_row$module_uid)
      if (!dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)) {
        stop("Could not create module annotation directory: ", module_dir, call. = FALSE)
      }
      .manifest_write_tsv(module_result, file.path(module_dir, "module_annotation.tsv"))
    }
  }

  result <- do.call(rbind, result_rows)
  rownames(result) <- NULL
  completed <- result[result$analysis_status == "completed", , drop = FALSE]
  if (nrow(completed)) {
    group_key <- paste(completed$module_uid, completed$database, sep = "\r")
    top <- do.call(rbind, lapply(split(completed, group_key), function(group) {
      group <- group[order(group$p_adjust, group$p_value, group$term_id), , drop = FALSE]
      utils::head(group, as.integer(top_n))
    }))
    rownames(top) <- NULL
  } else {
    top <- .annotation_empty_result()
  }
  qc <- do.call(rbind, qc_rows)
  rownames(qc) <- NULL

  annotation_file <- file.path(staging_dir, "module_annotation.tsv")
  top_file <- file.path(staging_dir, paste0("module_annotation_top", as.integer(top_n), ".tsv"))
  qc_file <- file.path(staging_dir, "module_annotation_qc.tsv")
  .manifest_write_tsv(result, annotation_file)
  .manifest_write_tsv(top, top_file)
  .manifest_write_tsv(qc, qc_file)
  .manifest_write_tsv(
    data.frame(
      parameter = c(
        "manifest_file", "gene_set_source", "databases", "species", "background_size",
        "pvalue_cutoff", "p_adjust_method", "p_adjust_cutoff", "qvalue_cutoff",
        "min_gene_set_size", "max_gene_set_size", "top_n"
      ),
      value = c(
        manifest_file,
        if (!is.null(gene_set_file)) normalizePath(gene_set_file, winslash = "/", mustWork = TRUE) else if (exclusive_inputs) "in_memory" else "msigdbr",
        paste(sort(unique(gene_sets$database)), collapse = ","),
        species,
        length(universe),
        pvalueCutoff,
        pAdjustMethod,
        pAdjustCutoff,
        qvalueCutoff,
        minGSSize,
        maxGSSize,
        top_n
      ),
      stringsAsFactors = FALSE
    ),
    file.path(staging_dir, "module_annotation_parameters.tsv")
  )

  if (!file.rename(staging_dir, output_dir)) {
    stop("Could not finalize annotation output directory: ", output_dir, call. = FALSE)
  }
  attr(result, "output_dir") <- output_dir
  attr(result, "annotation_file") <- file.path(output_dir, basename(annotation_file))
  attr(result, "top_file") <- file.path(output_dir, basename(top_file))
  attr(result, "qc_file") <- file.path(output_dir, basename(qc_file))
  result
}
