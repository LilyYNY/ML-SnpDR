# Step 2: subtype-specific PPI networks -------------------------------------

.network_read_diff <- function(path) {
  extension <- tolower(tools::file_ext(path))
  if (extension == "xlsx") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Reading legacy differential-expression Excel files requires openxlsx.", call. = FALSE)
    }
    sheets <- openxlsx::getSheetNames(path)
    sheets <- sheets[grepl("_DiffResults$", sheets)]
    if (!length(sheets)) stop("No *_DiffResults sheets were found in diff_file.", call. = FALSE)
    rows <- lapply(sheets, function(sheet) {
      data <- .upstream_read_tabular(path, sheet = sheet)
      gene <- .upstream_column(data, "Gene", c("gene", "protein", "symbol"), "Gene")
      label <- .upstream_column(data, "label", c("change", "direction"), "Label")
      data.frame(
        subtype = sub("_DiffResults$", "", sheet),
        gene = trimws(as.character(data[[gene]])),
        label = tolower(gsub("-", "_", trimws(as.character(data[[label]])))),
        stringsAsFactors = FALSE
      )
    })
    result <- do.call(rbind, rows)
  } else {
    data <- .upstream_read_tabular(path)
    subtype <- .upstream_column(data, "subtype", c("class", "group"), "Subtype")
    gene <- .upstream_column(data, "gene", c("Gene", "protein", "symbol"), "Gene")
    label <- .upstream_column(data, "label", c("change", "direction"), "Label")
    result <- data.frame(
      subtype = as.character(data[[subtype]]),
      gene = trimws(as.character(data[[gene]])),
      label = tolower(gsub("-", "_", trimws(as.character(data[[label]])))),
      stringsAsFactors = FALSE
    )
  }
  result$subtype <- .normalize_subtype(result$subtype)
  result$label[result$label %in% c("non_significant", "nonsignificant", "ns", "not")] <- "non_significant"
  if (any(!nzchar(result$gene))) stop("Differential-expression input contains empty gene identifiers.", call. = FALSE)
  unique(result)
}

.network_source_index <- function(ppi_index_file, ppi_sources) {
  if (!xor(is.null(ppi_index_file), is.null(ppi_sources))) {
    stop("Supply exactly one of ppi_index_file or ppi_sources.", call. = FALSE)
  }
  if (!is.null(ppi_sources)) {
    if (is.null(names(ppi_sources)) || any(!nzchar(names(ppi_sources)))) {
      stop("ppi_sources must be a named vector such as c(String='string.tsv').", call. = FALSE)
    }
    index <- data.frame(
      network = names(ppi_sources),
      edge_file = as.character(ppi_sources),
      stringsAsFactors = FALSE
    )
    root <- getwd()
  } else {
    index <- .upstream_read_tabular(ppi_index_file)
    required <- c("network", "edge_file")
    missing <- setdiff(required, names(index))
    if (length(missing)) stop("PPI index is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    root <- dirname(normalizePath(ppi_index_file, winslash = "/", mustWork = TRUE))
  }
  optional <- c(
    "node1_column", "node2_column", "score_column", "score_min", "delimiter",
    "mapping_file", "mapping_id_column", "mapping_symbol_column", "mapping_delimiter"
  )
  for (column in setdiff(optional, names(index))) index[[column]] <- NA
  index$network <- normalize_network_name(index$network, output = "display")
  if (anyDuplicated(index$network)) stop("PPI index must contain one row per network.", call. = FALSE)
  resolve <- function(path) {
    if (is.na(path) || !nzchar(trimws(path))) return(NA_character_)
    path <- path.expand(as.character(path))
    if (!grepl("^[A-Za-z]:[/\\\\]|^/", path)) path <- file.path(root, path)
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }
  index$edge_file <- vapply(index$edge_file, resolve, character(1))
  index$mapping_file <- vapply(index$mapping_file, resolve, character(1))
  index
}

.network_read_source_table <- function(path, delimiter = NA_character_) {
  extension <- tolower(tools::file_ext(path))
  if (extension %in% c("xlsx", "csv")) return(.upstream_read_tabular(path))
  if (is.na(delimiter) || !nzchar(delimiter) || delimiter %in% c("tab", "\\t")) {
    return(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = ""))
  }
  if (tolower(delimiter) %in% c("space", "whitespace")) delimiter <- ""
  utils::read.table(
    path,
    header = TRUE,
    sep = delimiter,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

.network_optional_column <- function(data, requested, aliases, required, label) {
  if (!is.na(requested) && nzchar(trimws(requested))) {
    if (!requested %in% names(data)) stop(label, " column does not exist: ", requested, call. = FALSE)
    return(requested)
  }
  hits <- match(tolower(aliases), tolower(names(data)), nomatch = 0L)
  hits <- hits[hits > 0L]
  if (length(hits)) return(names(data)[hits[[1L]]])
  if (isTRUE(required)) stop(label, " column could not be identified.", call. = FALSE)
  NA_character_
}

.network_standardize_source <- function(row, default_score_min) {
  edges <- .network_read_source_table(row$edge_file, row$delimiter)
  node1_column <- .network_optional_column(
    edges, row$node1_column,
    c("node1", "protein1", "protein.a", "protein_a", "gene1", "from"),
    TRUE, "PPI node1"
  )
  node2_column <- .network_optional_column(
    edges, row$node2_column,
    c("node2", "protein2", "protein.b", "protein_b", "gene2", "to"),
    TRUE, "PPI node2"
  )
  score_column <- .network_optional_column(
    edges, row$score_column,
    c("score", "combined_score", "combined.score", "confidence"),
    FALSE, "PPI score"
  )
  result <- data.frame(
    node1 = trimws(as.character(edges[[node1_column]])),
    node2 = trimws(as.character(edges[[node2_column]])),
    stringsAsFactors = FALSE
  )
  score_min <- suppressWarnings(as.numeric(as.character(row$score_min)))
  if (!is.finite(score_min)) score_min <- default_score_min
  if (!is.na(score_column)) {
    score <- suppressWarnings(as.numeric(as.character(edges[[score_column]])))
    if (any(is.na(score) & !is.na(edges[[score_column]]))) {
      stop("PPI score column contains non-numeric values for ", row$network, ".", call. = FALSE)
    }
    result <- result[is.finite(score) & score > score_min, , drop = FALSE]
  }

  if (!is.na(row$mapping_file)) {
    mapping <- .network_read_source_table(row$mapping_file, row$mapping_delimiter)
    id_column <- .network_optional_column(
      mapping, row$mapping_id_column,
      c("id", "protein_id", "protein.id", "geneid", "gene_id"),
      TRUE, "Mapping ID"
    )
    symbol_column <- .network_optional_column(
      mapping, row$mapping_symbol_column,
      c("symbol", "gene", "gene_symbol", "preferred_name"),
      TRUE, "Mapping symbol"
    )
    ids <- trimws(as.character(mapping[[id_column]]))
    symbols <- trimws(as.character(mapping[[symbol_column]]))
    if (anyDuplicated(ids)) {
      conflicts <- vapply(split(symbols, ids), function(x) length(unique(x[nzchar(x)])) > 1L, logical(1))
      if (any(conflicts)) stop("PPI mapping contains conflicting duplicate IDs.", call. = FALSE)
      keep <- !duplicated(ids)
      ids <- ids[keep]
      symbols <- symbols[keep]
    }
    result$node1 <- symbols[match(result$node1, ids)]
    result$node2 <- symbols[match(result$node2, ids)]
  }
  result <- result[stats::complete.cases(result) & nzchar(result$node1) & nzchar(result$node2), , drop = FALSE]
  result <- result[result$node1 != result$node2, , drop = FALSE]
  if (nrow(result)) {
    endpoints <- t(apply(result, 1L, sort))
    result$node1 <- endpoints[, 1L]
    result$node2 <- endpoints[, 2L]
    result <- unique(result[order(result$node1, result$node2), , drop = FALSE])
  }
  rownames(result) <- NULL
  result
}

#' Construct subtype-specific PPI networks from differential proteins
#'
#' Step 2 preserves subnetDR's rule that both endpoints must be subtype-up
#' proteins, while replacing fixed `./PPI` paths with a named source vector or
#' an explicit PPI index. Raw ID networks may declare a mapping file in the
#' index. The canonical network manifest is passed directly to step 3.
#'
#' @param diff_file Step-1 canonical TSV or subnetDR Excel workbook.
#' @param ppi_index_file Optional PPI index with `network` and `edge_file` plus
#'   optional column, score, delimiter and mapping metadata.
#' @param ppi_sources Optional named vector of already symbol-level edge files.
#' @param output_dir New step-2 output directory.
#' @param ppi_method Networks to construct.
#' @param ppiScore Default strict lower score threshold when a source declares
#'   a score column but no row-specific `score_min`.
#' @param include_labels Differential-expression labels used as network nodes.
#' @param subtypes Optional subtype subset.
#' @param strict Require every requested network source and subtype.
#' @return Network manifest with paths stored in attributes.
#' @export
run_network_construction <- function(
  diff_file,
  ppi_index_file = NULL,
  ppi_sources = NULL,
  output_dir = "./Netconstruct_result/",
  ppi_method = c("String", "physicalPPIN", "chengF"),
  ppiScore = 400,
  include_labels = "up",
  subtypes = NULL,
  strict = TRUE
) {
  if (length(ppiScore) != 1L || !is.finite(ppiScore)) stop("ppiScore must be one finite number.", call. = FALSE)
  if (!length(include_labels) || any(!nzchar(trimws(include_labels)))) {
    stop("include_labels must contain at least one label.", call. = FALSE)
  }
  include_labels <- tolower(gsub("-", "_", trimws(include_labels)))
  diff <- .network_read_diff(diff_file)
  sources <- .network_source_index(ppi_index_file, ppi_sources)
  networks <- unique(normalize_network_name(ppi_method, output = "display"))
  missing_networks <- setdiff(networks, sources$network)
  if (length(missing_networks) && isTRUE(strict)) {
    stop("PPI source(s) missing for: ", paste(missing_networks, collapse = ", "), call. = FALSE)
  }
  networks <- intersect(networks, sources$network)
  if (!length(networks)) stop("No requested PPI source is available.", call. = FALSE)
  if (is.null(subtypes)) subtypes <- unique(diff$subtype)
  subtypes <- sort(unique(.normalize_subtype(subtypes)))
  absent <- setdiff(subtypes, diff$subtype)
  if (length(absent)) stop("Differential-expression results lack subtype(s): ", paste(absent, collapse = ", "), call. = FALSE)

  standardized <- stats::setNames(vector("list", length(networks)), networks)
  for (network in networks) {
    standardized[[network]] <- .network_standardize_source(
      sources[sources$network == network, , drop = FALSE],
      ppiScore
    )
  }

  output_target <- .upstream_output_target(output_dir, "Network-construction")
  staging_dir <- tempfile(pattern = paste0(".", basename(output_target), "-"), tmpdir = dirname(output_target))
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create network-construction staging directory.", call. = FALSE)
  }
  on.exit(if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
  manifest_rows <- list()
  row_index <- 0L
  for (subtype in subtypes) {
    proteins <- sort(unique(diff$gene[diff$subtype == subtype & diff$label %in% include_labels]))
    if (!length(proteins) && isTRUE(strict)) stop("No selected differential proteins for ", subtype, ".", call. = FALSE)
    for (network in networks) {
      source_edges <- standardized[[network]]
      selected <- source_edges[
        source_edges$node1 %in% proteins & source_edges$node2 %in% proteins,
        , drop = FALSE
      ]
      relative_dir <- file.path(subtype, network)
      if (!dir.create(file.path(staging_dir, relative_dir), recursive = TRUE, showWarnings = FALSE)) {
        stop("Could not create subtype-network output directory.", call. = FALSE)
      }
      relative_file <- .upstream_relative(file.path(relative_dir, paste0("ppi_", subtype, ".txt")))
      absolute_file <- file.path(staging_dir, relative_file)
      .manifest_write_tsv(selected, absolute_file)
      row_index <- row_index + 1L
      source_row <- sources[sources$network == network, , drop = FALSE]
      manifest_rows[[row_index]] <- data.frame(
        network_uid = paste0(normalize_network_name(network, output = "slug"), "__", subtype),
        network = network,
        subtype = subtype,
        selected_labels = paste(include_labels, collapse = ","),
        selected_protein_number = length(proteins),
        source_edge_number = nrow(source_edges),
        node_count = length(unique(c(selected$node1, selected$node2))),
        edge_count = nrow(selected),
        ppi_file = relative_file,
        ppi_sha256 = digest::digest(absolute_file, algo = "sha256", file = TRUE),
        source_edge_file = source_row$edge_file,
        analysis_status = if (nrow(selected)) "completed" else "completed_empty_network",
        status_reason = if (nrow(selected)) "" else "No source PPI edge joined two selected differential proteins.",
        stringsAsFactors = FALSE
      )
    }
  }
  manifest <- do.call(rbind, manifest_rows)
  rownames(manifest) <- NULL
  manifest_file <- file.path(staging_dir, "network_manifest.tsv")
  .manifest_write_tsv(manifest, manifest_file)
  if (!file.rename(staging_dir, output_target)) {
    stop("Could not finalize network-construction output directory.", call. = FALSE)
  }
  attr(manifest, "output_dir") <- output_target
  attr(manifest, "manifest_file") <- file.path(output_target, basename(manifest_file))
  manifest
}
