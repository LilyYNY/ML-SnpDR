# Step 3: module division ----------------------------------------------------

.division_read_ppi <- function(path) {
  data <- .upstream_read_tabular(path)
  node1 <- .upstream_column(data, "node1", c("gene1", "protein1", "from"), "PPI node1")
  node2 <- .upstream_column(data, "node2", c("gene2", "protein2", "to"), "PPI node2")
  result <- data.frame(
    node1 = trimws(as.character(data[[node1]])),
    node2 = trimws(as.character(data[[node2]])),
    stringsAsFactors = FALSE
  )
  result <- result[nzchar(result$node1) & nzchar(result$node2) & result$node1 != result$node2, , drop = FALSE]
  if (nrow(result)) {
    endpoints <- t(apply(result, 1L, sort))
    result$node1 <- endpoints[, 1L]
    result$node2 <- endpoints[, 2L]
    result <- unique(result[order(result$node1, result$node2), , drop = FALSE])
  }
  rownames(result) <- NULL
  result
}

.division_louvain <- function(graph, seed) {
  set.seed(seed)
  clustering <- igraph::cluster_louvain(graph, weights = NULL, resolution = 1)
  data.frame(
    node = names(igraph::membership(clustering)),
    module = as.integer(igraph::membership(clustering)),
    stringsAsFactors = FALSE
  )
}

.division_wf <- function(graph, seed, pvalue_cutoff) {
  set.seed(seed)
  gn <- igraph::cluster_edge_betweenness(graph, weights = NULL)
  set.seed(seed)
  lp <- igraph::cluster_label_prop(graph, weights = NULL)
  gn_groups <- split(names(igraph::membership(gn)), igraph::membership(gn))
  lp_groups <- split(names(igraph::membership(lp)), igraph::membership(lp))
  gn_groups <- gn_groups[lengths(gn_groups) > 1L]
  lp_groups <- lp_groups[lengths(lp_groups) > 1L]
  if (!length(gn_groups) || !length(lp_groups)) {
    stop("WF could not form non-singleton GN and label-propagation groups.", call. = FALSE)
  }
  universe_size <- igraph::vcount(graph)
  candidates <- list()
  candidate_index <- 0L
  for (gn_index in seq_along(gn_groups)) {
    for (lp_index in seq_along(lp_groups)) {
      overlap <- intersect(gn_groups[[gn_index]], lp_groups[[lp_index]])
      if (length(overlap) < 2L) next
      p_value <- stats::phyper(
        length(overlap) - 1L,
        length(gn_groups[[gn_index]]),
        universe_size - length(gn_groups[[gn_index]]),
        length(lp_groups[[lp_index]]),
        lower.tail = FALSE
      )
      if (is.finite(p_value) && p_value < pvalue_cutoff) {
        candidate_index <- candidate_index + 1L
        candidates[[candidate_index]] <- data.frame(
          gn_index = gn_index,
          lp_index = lp_index,
          overlap_number = length(overlap),
          p_value = p_value,
          nodes = I(list(sort(overlap))),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(candidates)) stop("WF found no significant GN/LP overlap module.", call. = FALSE)
  candidate_table <- do.call(rbind, candidates)
  ordering <- order(candidate_table$p_value, -candidate_table$overlap_number,
                    candidate_table$gn_index, candidate_table$lp_index)
  candidate_table <- candidate_table[ordering, , drop = FALSE]
  assigned <- character()
  modules <- list()
  for (index in seq_len(nrow(candidate_table))) {
    nodes <- setdiff(candidate_table$nodes[[index]], assigned)
    if (length(nodes) < 2L) next
    modules[[length(modules) + 1L]] <- nodes
    assigned <- c(assigned, nodes)
  }
  if (!length(modules)) stop("WF significant overlaps were not disjoint after assignment.", call. = FALSE)
  do.call(rbind, lapply(seq_along(modules), function(index) {
    data.frame(node = modules[[index]], module = index, stringsAsFactors = FALSE)
  }))
}

.division_annotate_edges <- function(edges, node_module) {
  membership <- stats::setNames(node_module$module, node_module$node)
  module1 <- unname(membership[edges$node1])
  module2 <- unname(membership[edges$node2])
  module <- ifelse(!is.na(module1) & !is.na(module2) & module1 == module2, module1, 0L)
  data.frame(node1 = edges$node1, node2 = edges$node2, module = as.integer(module), stringsAsFactors = FALSE)
}

#' Divide one subtype PPI network into modules
#'
#' @param ppi_file_path Symbol-level edge table with `node1,node2`.
#' @param output_base_path Existing output directory for this group.
#' @param network_method Network name.
#' @param module_method `Louvain` or `WF`.
#' @param seed Reproducibility seed.
#' @param wf_pvalue_cutoff Hypergeometric overlap threshold for WF consensus.
#' @return One-row processing summary with absolute output paths.
#' @export
subtype_module <- function(
  ppi_file_path,
  output_base_path,
  network_method,
  module_method,
  seed = 123L,
  wf_pvalue_cutoff = 0.05
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Step 3 requires the igraph package.", call. = FALSE)
  }
  network <- normalize_network_name(network_method, output = "display")
  method <- normalize_method_name(module_method, output = "display")
  if (length(wf_pvalue_cutoff) != 1L || !is.finite(wf_pvalue_cutoff) ||
      wf_pvalue_cutoff <= 0 || wf_pvalue_cutoff > 1) {
    stop("wf_pvalue_cutoff must be in (0, 1].", call. = FALSE)
  }
  edges <- .division_read_ppi(ppi_file_path)
  if (!nrow(edges)) stop("PPI network contains no usable edges: ", ppi_file_path, call. = FALSE)
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE)
  graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)
  node_module <- if (method == "Louvain") {
    .division_louvain(graph, as.integer(seed))
  } else {
    .division_wf(graph, as.integer(seed), wf_pvalue_cutoff)
  }
  node_module <- node_module[order(node_module$module, node_module$node), , drop = FALSE]
  if (anyDuplicated(node_module$node)) stop("Module division assigned a node more than once.", call. = FALSE)
  edge_module <- .division_annotate_edges(edges, node_module)
  intramodule <- edge_module[edge_module$module != 0L, , drop = FALSE]
  if (!dir.exists(output_base_path) && !dir.create(output_base_path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create module-division group directory.", call. = FALSE)
  }
  node_file <- file.path(output_base_path, paste0("node_Module_", network, "_", method, ".txt"))
  edge_file <- file.path(output_base_path, paste0("edges_", network, "_", method, ".txt"))
  intramodule_file <- file.path(output_base_path, paste0("edge_Module_", network, "_", method, ".txt"))
  .manifest_write_tsv(node_module, node_file)
  .manifest_write_tsv(edge_module, edge_file)
  .manifest_write_tsv(intramodule, intramodule_file)
  data.frame(
    network = network,
    method = method,
    module_count = length(unique(node_module$module)),
    node_count = nrow(node_module),
    edge_count = nrow(edge_module),
    intramodule_edge_count = nrow(intramodule),
    node_module_file_abs = normalizePath(node_file, winslash = "/", mustWork = TRUE),
    edge_module_file_abs = normalizePath(edge_file, winslash = "/", mustWork = TRUE),
    intramodule_edge_file_abs = normalizePath(intramodule_file, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

.division_legacy_index <- function(base_input_path, networks, subtypes, strict) {
  root <- normalizePath(base_input_path, winslash = "/", mustWork = TRUE)
  rows <- list()
  row_index <- 0L
  for (subtype in subtypes) {
    subtype_dir <- file.path(root, subtype)
    candidates <- if (dir.exists(subtype_dir)) list.dirs(subtype_dir, recursive = FALSE, full.names = TRUE) else character()
    candidate_names <- basename(candidates)
    candidate_normalized <- vapply(candidate_names, function(x) {
      tryCatch(normalize_network_name(x, output = "display"), error = function(e) NA_character_)
    }, character(1))
    for (network in networks) {
      hit <- which(candidate_normalized == network)
      file <- if (length(hit)) file.path(candidates[hit[[1L]]], paste0("ppi_", subtype, ".txt")) else NA_character_
      if (is.na(file) || !file.exists(file)) {
        if (isTRUE(strict)) stop("Legacy PPI file missing for ", subtype, "/", network, ".", call. = FALSE)
        next
      }
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        network = network,
        subtype = subtype,
        ppi_file_abs = normalizePath(file, winslash = "/", mustWork = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) stop("No PPI networks were discovered for module division.", call. = FALSE)
  do.call(rbind, rows)
}

.division_manifest_index <- function(path) {
  manifest <- .upstream_read_tabular(path)
  required <- c("network", "subtype", "ppi_file")
  missing <- setdiff(required, names(manifest))
  if (length(missing)) stop("network_manifest.tsv is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  root <- dirname(normalizePath(path, winslash = "/", mustWork = TRUE))
  manifest$network <- normalize_network_name(manifest$network, output = "display")
  manifest$subtype <- .normalize_subtype(manifest$subtype)
  manifest$ppi_file_abs <- vapply(manifest$ppi_file, .drug_absolute_path, character(1), root = root)
  if (any(!file.exists(manifest$ppi_file_abs))) stop("network_manifest.tsv references missing PPI files.", call. = FALSE)
  manifest[, c("network", "subtype", "ppi_file_abs"), drop = FALSE]
}

#' Detect Louvain or WF modules in subtype-specific PPI networks
#'
#' Step 3 accepts the canonical step-2 network manifest or a legacy
#' Netconstruct_result root. It writes subnetDR-compatible node/edge module
#' files plus a canonical module-division manifest consumed by step 4.
#'
#' @param subtype_file Optional phenotype file used only to obtain subtype
#'   names for legacy directory input.
#' @param base_input_path Legacy Netconstruct_result root.
#' @param output_base_path New step-3 output directory.
#' @param network_method,module_method Requested networks and algorithms.
#' @param network_manifest_file Optional step-2 `network_manifest.tsv`.
#' @param subtypes Optional explicit subtype vector.
#' @param subtype_column Phenotype subtype column.
#' @param seed Reproducibility seed.
#' @param wf_pvalue_cutoff WF overlap threshold.
#' @param strict Require every requested subtype/network group.
#' @return Module-division manifest with output paths in attributes.
#' @export
module_division <- function(
  subtype_file = NULL,
  base_input_path = "./Netconstruct_result/",
  output_base_path = "./ModuleDivision/",
  network_method = c("String", "physicalPPIN", "chengF"),
  module_method = c("Louvain", "WF"),
  network_manifest_file = NULL,
  subtypes = NULL,
  subtype_column = "Subtype",
  seed = 123L,
  wf_pvalue_cutoff = 0.05,
  strict = TRUE
) {
  networks <- unique(normalize_network_name(network_method, output = "display"))
  methods <- unique(normalize_method_name(module_method, output = "display"))
  index <- if (!is.null(network_manifest_file)) {
    .division_manifest_index(network_manifest_file)
  } else {
    if (!is.null(subtype_file) && !is.null(subtypes)) {
      stop("Supply subtype_file or subtypes, not both.", call. = FALSE)
    }
    if (!is.null(subtype_file)) subtypes <- .selection_read_subtypes(subtype_file, subtype_column)
    if (is.null(subtypes)) {
      root_dirs <- list.dirs(base_input_path, recursive = FALSE, full.names = FALSE)
      subtypes <- root_dirs[grepl("^C[1-4]$", toupper(root_dirs))]
    }
    subtypes <- sort(unique(.normalize_subtype(subtypes)))
    .division_legacy_index(base_input_path, networks, subtypes, strict)
  }
  requested_subtypes <- if (is.null(subtypes)) sort(unique(index$subtype)) else
    sort(unique(.normalize_subtype(subtypes)))
  index <- index[index$network %in% networks, , drop = FALSE]
  index <- index[index$subtype %in% requested_subtypes, , drop = FALSE]
  if (!nrow(index)) stop("No step-3 network groups remain after filtering.", call. = FALSE)
  expected <- expand.grid(
    network = networks,
    subtype = requested_subtypes,
    stringsAsFactors = FALSE
  )
  observed_key <- paste(index$network, index$subtype, sep = "\r")
  expected_key <- paste(expected$network, expected$subtype, sep = "\r")
  missing_groups <- expected[!expected_key %in% observed_key, , drop = FALSE]
  if (nrow(missing_groups) && isTRUE(strict)) {
    stop("Step 3 input lacks requested network/subtype group(s).", call. = FALSE)
  }

  output_target <- .upstream_output_target(output_base_path, "Module-division")
  staging_dir <- tempfile(pattern = paste0(".", basename(output_target), "-"), tmpdir = dirname(output_target))
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create module-division staging directory.", call. = FALSE)
  }
  on.exit(if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
  manifest_rows <- list()
  row_index <- 0L
  for (source_index in seq_len(nrow(index))) {
    source <- index[source_index, , drop = FALSE]
    for (method in methods) {
      relative_dir <- file.path(source$network, method, source$subtype)
      group_dir <- file.path(staging_dir, relative_dir)
      result <- tryCatch(
        subtype_module(
          ppi_file_path = source$ppi_file_abs,
          output_base_path = group_dir,
          network_method = source$network,
          module_method = method,
          seed = seed,
          wf_pvalue_cutoff = wf_pvalue_cutoff
        ),
        error = function(error) {
          if (isTRUE(strict)) stop(error)
          warning(conditionMessage(error), call. = FALSE)
          NULL
        }
      )
      if (is.null(result)) next
      row_index <- row_index + 1L
      relative_node <- .upstream_relative(file.path(relative_dir, basename(result$node_module_file_abs)))
      relative_edge <- .upstream_relative(file.path(relative_dir, basename(result$edge_module_file_abs)))
      relative_intramodule <- .upstream_relative(file.path(relative_dir, basename(result$intramodule_edge_file_abs)))
      manifest_rows[[row_index]] <- data.frame(
        division_uid = paste(
          normalize_network_name(source$network, output = "slug"),
          normalize_method_name(method, output = "slug"),
          source$subtype,
          sep = "__"
        ),
        network = source$network,
        method = method,
        subtype = source$subtype,
        module_count = result$module_count,
        node_count = result$node_count,
        edge_count = result$edge_count,
        intramodule_edge_count = result$intramodule_edge_count,
        node_module_file = relative_node,
        edge_module_file = relative_edge,
        intramodule_edge_file = relative_intramodule,
        node_sha256 = digest::digest(result$node_module_file_abs, algo = "sha256", file = TRUE),
        edge_sha256 = digest::digest(result$edge_module_file_abs, algo = "sha256", file = TRUE),
        source_ppi_file = source$ppi_file_abs,
        seed = as.integer(seed),
        wf_pvalue_cutoff = wf_pvalue_cutoff,
        analysis_status = "completed",
        status_reason = "",
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(manifest_rows)) stop("No module-division group completed.", call. = FALSE)
  manifest <- do.call(rbind, manifest_rows)
  manifest <- manifest[order(manifest$network, manifest$method, manifest$subtype), , drop = FALSE]
  rownames(manifest) <- NULL
  manifest_file <- file.path(staging_dir, "module_division_manifest.tsv")
  .manifest_write_tsv(manifest, manifest_file)
  if (!file.rename(staging_dir, output_target)) {
    stop("Could not finalize module-division output directory.", call. = FALSE)
  }
  attr(manifest, "output_dir") <- output_target
  attr(manifest, "manifest_file") <- file.path(output_target, basename(manifest_file))
  manifest
}
