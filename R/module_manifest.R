# ModuleSelection adapter -----------------------------------------------------

.manifest_required_columns <- function(data, required, path) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Missing required column(s) in ", path, ": ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  data[, required, drop = FALSE]
}

.manifest_read_table <- function(path, sep, required) {
  data <- utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM",
    colClasses = "character"
  )
  .manifest_required_columns(data, required, path)
}

.manifest_relative_path <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  prefix <- paste0(root, "/")

  comparable_path <- if (.Platform$OS.type == "windows") tolower(path) else path
  comparable_prefix <- if (.Platform$OS.type == "windows") tolower(prefix) else prefix
  if (!startsWith(comparable_path, comparable_prefix)) {
    stop("Path is outside module_selection_dir: ", path, call. = FALSE)
  }
  substring(path, nchar(prefix) + 1L)
}

.manifest_empty <- function() {
  data.frame(
    module_uid = character(),
    legacy_module_id = character(),
    network = character(),
    method = character(),
    subtype = character(),
    module = character(),
    module_size = integer(),
    edge_count = integer(),
    self_loops_excluded = integer(),
    node_file = character(),
    edge_file = character(),
    prefilter_pass = logical(),
    prefilter_reason = character(),
    node_sha256 = character(),
    edge_sha256 = character(),
    source_module_file = character(),
    source_node_file = character(),
    source_edge_file = character(),
    stringsAsFactors = FALSE
  )
}

.manifest_write_tsv <- function(data, path) {
  utils::write.table(
    data,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "",
    eol = "\n"
  )
}

.manifest_sort <- function(manifest) {
  if (!nrow(manifest)) return(manifest)
  network_order <- match(manifest$network, c("String", "physicalPPIN", "chengF"))
  method_order <- match(manifest$method, c("Louvain", "WF"))
  module_order <- as.integer(sub("^M", "", manifest$module))
  manifest[
    order(network_order, method_order, manifest$subtype, module_order),
    ,
    drop = FALSE
  ]
}

#' Build a canonical module manifest from ModuleSelection outputs
#'
#' The expected input hierarchy is
#' `<network>/<method>/<subtype>/Module_select_<network>_<subtype>.csv`,
#' accompanied by the matching aggregate node and edge files. The adapter
#' validates the source tables, splits them into one node and edge file per
#' module, computes SHA256 checksums and writes `module_manifest.tsv` plus a
#' group-level QC table.
#'
#' @param module_selection_dir Root directory containing ModuleSelection output.
#' @param output_dir New output directory. It must not already exist.
#' @param min_size_exclusive Module sizes must be strictly greater than this
#'   value. The paper workflow uses 9.
#' @return The manifest data frame, with `manifest_file`, `qc_file` and
#'   `output_dir` attributes.
#' @export
build_module_manifest <- function(
  module_selection_dir,
  output_dir,
  min_size_exclusive = 9
) {
  if (length(module_selection_dir) != 1L || !dir.exists(module_selection_dir)) {
    stop("module_selection_dir must be one existing directory.", call. = FALSE)
  }
  if (length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one non-empty path.", call. = FALSE)
  }
  if (
    length(min_size_exclusive) != 1L ||
      !is.numeric(min_size_exclusive) ||
      is.na(min_size_exclusive) ||
      min_size_exclusive < 0
  ) {
    stop("min_size_exclusive must be one non-negative number.", call. = FALSE)
  }

  module_selection_dir <- normalizePath(
    module_selection_dir,
    winslash = "/",
    mustWork = TRUE
  )

  output_dir <- path.expand(output_dir)
  output_parent <- dirname(output_dir)
  if (!dir.exists(output_parent)) {
    if (!dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create output parent directory: ", output_parent, call. = FALSE)
    }
  }
  output_parent <- normalizePath(output_parent, winslash = "/", mustWork = TRUE)
  output_dir <- file.path(output_parent, basename(output_dir))
  if (file.exists(output_dir)) {
    stop(
      "output_dir already exists; choose a new directory to avoid overwriting: ",
      output_dir,
      call. = FALSE
    )
  }

  module_files <- sort(list.files(
    module_selection_dir,
    pattern = "^Module_select_.+[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  ))
  if (!length(module_files)) {
    stop("No Module_select_*.csv files found under module_selection_dir.", call. = FALSE)
  }

  staging_dir <- tempfile(
    pattern = paste0(".", basename(output_dir), "-"),
    tmpdir = output_parent
  )
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create staging directory: ", staging_dir, call. = FALSE)
  }
  on.exit({
    if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  manifest_rows <- list()
  qc_rows <- list()
  row_index <- 0L
  qc_index <- 0L

  for (module_file in module_files) {
    relative_module_file <- .manifest_relative_path(module_file, module_selection_dir)
    parts <- strsplit(relative_module_file, "/", fixed = TRUE)[[1]]
    if (length(parts) != 4L) {
      stop(
        "Expected <network>/<method>/<subtype>/<file>, found: ",
        relative_module_file,
        call. = FALSE
      )
    }

    network_source <- parts[1]
    method_source <- parts[2]
    subtype_source <- parts[3]
    network <- normalize_network_name(network_source, output = "display")
    method <- normalize_method_name(method_source, output = "display")
    subtype <- .normalize_subtype(subtype_source)
    source_dir <- dirname(module_file)

    expected_module_file <- file.path(
      source_dir,
      paste0("Module_select_", network_source, "_", subtype_source, ".csv")
    )
    node_file <- file.path(
      source_dir,
      paste0("node_Module_select_", network_source, "_", subtype_source, ".txt")
    )
    edge_file <- file.path(
      source_dir,
      paste0("edges_select_", network_source, "_", subtype_source, ".txt")
    )
    required_files <- c(expected_module_file, node_file, edge_file)
    missing_files <- required_files[!file.exists(required_files)]
    if (length(missing_files)) {
      stop(
        "Missing source file(s) for ", relative_module_file, ": ",
        paste(basename(missing_files), collapse = ", "),
        call. = FALSE
      )
    }
    if (
      normalizePath(module_file, winslash = "/", mustWork = TRUE) !=
        normalizePath(expected_module_file, winslash = "/", mustWork = TRUE)
    ) {
      stop("Unexpected module summary filename: ", module_file, call. = FALSE)
    }

    modules <- .manifest_read_table(module_file, ",", c("module", "count"))
    nodes <- .manifest_read_table(node_file, "\t", c("node", "module"))
    edges <- .manifest_read_table(edge_file, "\t", c("node1", "node2", "module"))

    modules$module <- .normalize_module(modules$module)
    counts <- suppressWarnings(as.numeric(modules$count))
    if (any(!is.finite(counts)) || any(counts < 1) || any(counts != floor(counts))) {
      stop("Module counts must be positive integers in: ", module_file, call. = FALSE)
    }
    modules$count <- as.integer(counts)
    if (anyDuplicated(modules$module)) {
      stop("Duplicate module labels in: ", module_file, call. = FALSE)
    }
    if (any(modules$count <= min_size_exclusive)) {
      failed <- modules$module[modules$count <= min_size_exclusive]
      stop(
        "ModuleSelection contains module(s) that fail size > ",
        min_size_exclusive, ": ", paste(failed, collapse = ", "),
        call. = FALSE
      )
    }

    if (nrow(nodes)) {
      nodes$node <- trimws(nodes$node)
      nodes$module <- .normalize_module(nodes$module)
      if (any(is.na(nodes$node)) || any(!nzchar(nodes$node))) {
        stop("Blank node identifier in: ", node_file, call. = FALSE)
      }
      if (anyDuplicated(nodes$node)) {
        duplicate_nodes <- unique(nodes$node[duplicated(nodes$node)])
        stop(
          "Node assigned more than once in ", node_file, ": ",
          paste(utils::head(duplicate_nodes, 10L), collapse = ", "),
          call. = FALSE
        )
      }
    }

    if (nrow(edges)) {
      edges$node1 <- trimws(edges$node1)
      edges$node2 <- trimws(edges$node2)
      edges$module <- toupper(trimws(edges$module))
      if (
        any(is.na(edges$node1)) || any(!nzchar(edges$node1)) ||
          any(is.na(edges$node2)) || any(!nzchar(edges$node2))
      ) {
        stop("Blank edge endpoint in: ", edge_file, call. = FALSE)
      }
      if (any(!grepl("^M[0-9]+$", edges$module))) {
        stop("Invalid edge module label in: ", edge_file, call. = FALSE)
      }
    }

    unexpected_node_modules <- setdiff(unique(nodes$module), modules$module)
    if (length(unexpected_node_modules)) {
      stop(
        "Node table contains unselected module(s) in ", node_file, ": ",
        paste(unexpected_node_modules, collapse = ", "),
        call. = FALSE
      )
    }
    unexpected_edge_modules <- setdiff(unique(edges$module), c("M0", modules$module))
    if (length(unexpected_edge_modules)) {
      stop(
        "Edge table contains unselected module(s) in ", edge_file, ": ",
        paste(unexpected_edge_modules, collapse = ", "),
        call. = FALSE
      )
    }

    source_module_path <- .manifest_relative_path(module_file, module_selection_dir)
    source_node_path <- .manifest_relative_path(node_file, module_selection_dir)
    source_edge_path <- .manifest_relative_path(edge_file, module_selection_dir)
    group_edge_total <- 0L

    if (nrow(modules)) {
      module_numbers <- as.integer(sub("^M", "", modules$module))
      modules <- modules[order(module_numbers), , drop = FALSE]
    }

    for (i in seq_len(nrow(modules))) {
      module <- modules$module[i]
      module_size <- modules$count[i]
      module_nodes <- nodes[nodes$module == module, c("node", "module"), drop = FALSE]
      module_edges <- edges[edges$module == module, c("node1", "node2", "module"), drop = FALSE]
      self_loop_count <- sum(module_edges$node1 == module_edges$node2)
      module_edges <- module_edges[module_edges$node1 != module_edges$node2, , drop = FALSE]

      if (nrow(module_nodes) != module_size) {
        stop(
          "Declared size does not match node count for ",
          paste(network, method, subtype, module, sep = "/"),
          ": declared=", module_size, ", observed=", nrow(module_nodes),
          call. = FALSE
        )
      }
      if (!nrow(module_edges)) {
        stop(
          "Selected module has no within-module edges: ",
          paste(network, method, subtype, module, sep = "/"),
          call. = FALSE
        )
      }

      outside_endpoints <- setdiff(
        unique(c(module_edges$node1, module_edges$node2)),
        module_nodes$node
      )
      if (length(outside_endpoints)) {
        stop(
          "Within-module edge contains node(s) outside ", module, ": ",
          paste(utils::head(outside_endpoints, 10L), collapse = ", "),
          call. = FALSE
        )
      }
      edge_keys <- paste(
        pmin(module_edges$node1, module_edges$node2),
        pmax(module_edges$node1, module_edges$node2),
        sep = "\r"
      )
      if (anyDuplicated(edge_keys)) {
        stop("Duplicate undirected edge in ", edge_file, " for ", module, call. = FALSE)
      }

      module_uid <- make_module_uid(network, method, subtype, module)
      legacy_module_id <- parse_module_uid(module_uid)$legacy_module_id
      relative_node_output <- file.path("modules", module_uid, "nodes.tsv")
      relative_edge_output <- file.path("modules", module_uid, "edges.tsv")
      module_output_dir <- file.path(staging_dir, "modules", module_uid)
      if (!dir.create(module_output_dir, recursive = TRUE, showWarnings = FALSE)) {
        stop("Could not create module output directory: ", module_output_dir, call. = FALSE)
      }

      module_nodes <- module_nodes[order(module_nodes$node), , drop = FALSE]
      module_edges <- module_edges[order(module_edges$node1, module_edges$node2), , drop = FALSE]
      output_node_file <- file.path(staging_dir, relative_node_output)
      output_edge_file <- file.path(staging_dir, relative_edge_output)
      .manifest_write_tsv(module_nodes, output_node_file)
      .manifest_write_tsv(module_edges, output_edge_file)

      row_index <- row_index + 1L
      group_edge_total <- group_edge_total + nrow(module_edges)
      manifest_rows[[row_index]] <- data.frame(
        module_uid = module_uid,
        legacy_module_id = legacy_module_id,
        network = network,
        method = method,
        subtype = subtype,
        module = module,
        module_size = as.integer(module_size),
        edge_count = as.integer(nrow(module_edges)),
        self_loops_excluded = as.integer(self_loop_count),
        node_file = gsub("\\\\", "/", relative_node_output),
        edge_file = gsub("\\\\", "/", relative_edge_output),
        prefilter_pass = TRUE,
        prefilter_reason = paste0("module_size > ", min_size_exclusive),
        node_sha256 = digest::digest(output_node_file, algo = "sha256", file = TRUE),
        edge_sha256 = digest::digest(output_edge_file, algo = "sha256", file = TRUE),
        source_module_file = source_module_path,
        source_node_file = source_node_path,
        source_edge_file = source_edge_path,
        stringsAsFactors = FALSE
      )
    }

    qc_index <- qc_index + 1L
    qc_rows[[qc_index]] <- data.frame(
      network = network,
      method = method,
      subtype = subtype,
      n_modules = as.integer(nrow(modules)),
      total_nodes = as.integer(sum(modules$count)),
      total_within_edges = as.integer(group_edge_total),
      excluded_self_loops = as.integer(sum(
        edges$module != "M0" & edges$node1 == edges$node2
      )),
      excluded_m0_edges = as.integer(sum(edges$module == "M0")),
      source_module_file = source_module_path,
      stringsAsFactors = FALSE
    )
  }

  manifest <- if (length(manifest_rows)) {
    do.call(rbind, manifest_rows)
  } else {
    .manifest_empty()
  }
  manifest <- .manifest_sort(manifest)
  rownames(manifest) <- NULL
  if (anyDuplicated(manifest$module_uid)) {
    duplicate_uids <- unique(manifest$module_uid[duplicated(manifest$module_uid)])
    stop(
      "Duplicate module_uid values across source directories: ",
      paste(duplicate_uids, collapse = ", "),
      call. = FALSE
    )
  }

  qc <- do.call(rbind, qc_rows)
  qc <- qc[
    order(
      match(qc$network, c("String", "physicalPPIN", "chengF")),
      match(qc$method, c("Louvain", "WF")),
      qc$subtype
    ),
    ,
    drop = FALSE
  ]
  rownames(qc) <- NULL

  .manifest_write_tsv(manifest, file.path(staging_dir, "module_manifest.tsv"))
  .manifest_write_tsv(qc, file.path(staging_dir, "module_manifest_qc.tsv"))

  if (!file.rename(staging_dir, output_dir)) {
    stop("Could not finalize output directory: ", output_dir, call. = FALSE)
  }

  attr(manifest, "output_dir") <- output_dir
  attr(manifest, "manifest_file") <- file.path(output_dir, "module_manifest.tsv")
  attr(manifest, "qc_file") <- file.path(output_dir, "module_manifest_qc.tsv")
  manifest
}

#' Read and validate a canonical module manifest
#'
#' Relative node and edge paths are resolved against the directory containing
#' `module_manifest.tsv`. The returned table retains the portable relative
#' columns and adds `node_file_abs` and `edge_file_abs` for downstream steps.
#'
#' @param path Path to `module_manifest.tsv`.
#' @param check_files Whether every referenced node and edge file must exist.
#' @param verify_hashes Whether to recompute and verify every SHA256 checksum.
#' @return A validated manifest data frame with resolved absolute file columns.
#' @export
read_module_manifest <- function(path, check_files = TRUE, verify_hashes = FALSE) {
  if (length(path) != 1L || !file.exists(path)) {
    stop("path must be one existing module_manifest.tsv file.", call. = FALSE)
  }
  if (!is.logical(check_files) || length(check_files) != 1L || is.na(check_files)) {
    stop("check_files must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verify_hashes) || length(verify_hashes) != 1L || is.na(verify_hashes)) {
    stop("verify_hashes must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(verify_hashes) && !isTRUE(check_files)) {
    stop("verify_hashes = TRUE requires check_files = TRUE.", call. = FALSE)
  }

  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  manifest <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  required <- names(.manifest_empty())
  manifest <- .manifest_required_columns(manifest, required, path)
  if (!nrow(manifest)) stop("Module manifest contains no modules: ", path, call. = FALSE)
  if (anyDuplicated(manifest$module_uid)) {
    stop("module_manifest.tsv contains duplicate module_uid values.", call. = FALSE)
  }

  identities <- parse_module_uid(manifest$module_uid)
  identity_fields <- c("network", "method", "subtype", "module", "legacy_module_id")
  for (field in identity_fields) {
    if (any(as.character(manifest[[field]]) != as.character(identities[[field]]))) {
      stop("Manifest identity field does not match module_uid: ", field, call. = FALSE)
    }
  }
  if (
    any(!is.finite(manifest$module_size)) ||
      any(manifest$module_size < 1) ||
      any(manifest$module_size != floor(manifest$module_size))
  ) {
    stop("module_size must contain positive integers.", call. = FALSE)
  }
  if (
    any(!is.finite(manifest$edge_count)) ||
      any(manifest$edge_count < 1) ||
      any(manifest$edge_count != floor(manifest$edge_count))
  ) {
    stop("edge_count must contain positive integers.", call. = FALSE)
  }

  portable_paths <- c(manifest$node_file, manifest$edge_file)
  absolute_pattern <- "^(?:[A-Za-z]:[/\\\\]|/|\\\\\\\\)"
  if (any(grepl(absolute_pattern, portable_paths))) {
    stop("Manifest node_file and edge_file values must be relative paths.", call. = FALSE)
  }
  manifest_root <- dirname(path)
  manifest$node_file_abs <- normalizePath(
    file.path(manifest_root, manifest$node_file),
    winslash = "/",
    mustWork = FALSE
  )
  manifest$edge_file_abs <- normalizePath(
    file.path(manifest_root, manifest$edge_file),
    winslash = "/",
    mustWork = FALSE
  )

  if (isTRUE(check_files)) {
    missing <- c(
      manifest$node_file_abs[!file.exists(manifest$node_file_abs)],
      manifest$edge_file_abs[!file.exists(manifest$edge_file_abs)]
    )
    if (length(missing)) {
      stop(
        "Manifest references missing module file(s): ",
        paste(utils::head(unique(missing), 10L), collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (isTRUE(verify_hashes)) {
    node_hashes <- vapply(
      manifest$node_file_abs,
      digest::digest,
      character(1),
      algo = "sha256",
      file = TRUE
    )
    edge_hashes <- vapply(
      manifest$edge_file_abs,
      digest::digest,
      character(1),
      algo = "sha256",
      file = TRUE
    )
    if (any(node_hashes != manifest$node_sha256)) {
      stop("Node SHA256 mismatch in module manifest.", call. = FALSE)
    }
    if (any(edge_hashes != manifest$edge_sha256)) {
      stop("Edge SHA256 mismatch in module manifest.", call. = FALSE)
    }
  }

  attr(manifest, "manifest_file") <- path
  attr(manifest, "manifest_root") <- manifest_root
  manifest
}
