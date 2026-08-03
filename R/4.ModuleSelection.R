# Step 4: module selection ----------------------------------------------------

.selection_source_modules <- function(x, allow_zero = FALSE, field = "module") {
  value <- toupper(trimws(as.character(x)))
  value <- ifelse(grepl("^M", value), value, paste0("M", value))
  pattern <- if (isTRUE(allow_zero)) "^M[0-9]+$" else "^M[1-9][0-9]*$"
  if (any(is.na(value)) || any(!grepl(pattern, value))) {
    stop("Invalid ", field, " label in ModuleDivision input.", call. = FALSE)
  }
  value
}

.selection_read_subtypes <- function(path, column) {
  if (!file.exists(path)) stop("Subtype file does not exist: ", path, call. = FALSE)
  extension <- tolower(tools::file_ext(path))
  data <- switch(
    extension,
    csv = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    tsv = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    txt = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    xlsx = {
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        stop("Reading .xlsx subtype files requires the openxlsx package.", call. = FALSE)
      }
      openxlsx::read.xlsx(path)
    },
    stop(
      "Unsupported subtype file extension: .", extension,
      ". Use .csv, .tsv, .txt or .xlsx.",
      call. = FALSE
    )
  )
  if (!column %in% names(data)) {
    stop("Subtype file is missing column: ", column, call. = FALSE)
  }
  unique(.normalize_subtype(data[[column]]))
}

.selection_discover_subtypes <- function(input_dir, networks, methods) {
  discovered <- character()
  for (network in networks) {
    for (method in methods) {
      group_dir <- file.path(input_dir, network, method)
      if (!dir.exists(group_dir)) next
      candidates <- list.dirs(group_dir, recursive = FALSE, full.names = FALSE)
      discovered <- c(discovered, candidates[grepl("^C[1-4]$", toupper(candidates))])
    }
  }
  unique(.normalize_subtype(discovered))
}

.selection_read_required <- function(path, required) {
  if (!file.exists(path)) stop("Input file does not exist: ", path, call. = FALSE)
  data <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Missing required column(s) in ", path, ": ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  data
}

.selection_output_target <- function(path) {
  path <- path.expand(path)
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    if (!dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create output parent directory: ", parent, call. = FALSE)
    }
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

.selection_plot_counts <- function(counts, selected, path) {
  colors <- ifelse(counts$module %in% selected, "lightslateblue", "gray")
  grDevices::tiff(path, width = 2500, height = 2500, res = 300)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::barplot(
    rev(counts$count),
    names.arg = rev(counts$module),
    horiz = TRUE,
    las = 1,
    col = rev(colors),
    border = NA,
    xlab = "Module size"
  )
}

#' Select ModuleDivision modules and create the canonical step-4 manifest
#'
#' This is a chaining-compatible replacement for subnetDR's
#' `module_selection()`. It preserves the original aggregate output layout and
#' additionally writes `standardized/module_manifest.tsv`, per-module node and
#' edge files, SHA256 checksums and QC. The returned manifest is the primary
#' input for steps 5, 6 and 6A.
#'
#' @param subtype_file Optional `.csv`, `.tsv`, `.txt` or `.xlsx` phenotype
#'   file used to obtain subtype names. If omitted, subtype directories are
#'   discovered from `base_input_path`.
#' @param base_input_path ModuleDivision root directory.
#' @param base_output_path New ModuleSelection output directory. It must not
#'   already exist.
#' @param network_method PPI networks to process.
#' @param module_method Module-detection methods to process.
#' @param numberCutoff Select modules with size strictly greater than this
#'   value.
#' @param subtypes Optional explicit subtype vector. It cannot be combined with
#'   `subtype_file`.
#' @param subtype_column Column containing subtype names when a subtype file is
#'   supplied.
#' @param strict When `TRUE`, missing ModuleDivision files stop the step. When
#'   `FALSE`, incomplete groups are skipped.
#' @param write_plots Whether to write the original module-size TIFF plots.
#' @return The canonical module manifest. Output paths are available from its
#'   `output_dir`, `manifest_file` and `qc_file` attributes.
#' @export
module_selection <- function(
  subtype_file = NULL,
  base_input_path = "./ModuleDivision/",
  base_output_path = "./ModuleSelection/",
  network_method = c("String", "physicalPPIN", "chengF"),
  module_method = c("Louvain", "WF"),
  numberCutoff = 9,
  subtypes = NULL,
  subtype_column = "Subtype",
  strict = TRUE,
  write_plots = FALSE
) {
  if (length(base_input_path) != 1L || !dir.exists(base_input_path)) {
    stop("base_input_path must be one existing ModuleDivision directory.", call. = FALSE)
  }
  if (
    length(numberCutoff) != 1L || !is.numeric(numberCutoff) ||
      is.na(numberCutoff) || numberCutoff < 0
  ) {
    stop("numberCutoff must be one non-negative number.", call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("strict must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(write_plots) || length(write_plots) != 1L || is.na(write_plots)) {
    stop("write_plots must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(subtype_file) && !is.null(subtypes)) {
    stop("Supply either subtype_file or subtypes, not both.", call. = FALSE)
  }

  input_dir <- normalizePath(base_input_path, winslash = "/", mustWork = TRUE)
  networks <- unique(normalize_network_name(network_method, output = "display"))
  methods <- unique(normalize_method_name(module_method, output = "display"))
  if (!is.null(subtype_file)) {
    subtypes <- .selection_read_subtypes(subtype_file, subtype_column)
  } else if (!is.null(subtypes)) {
    subtypes <- unique(.normalize_subtype(subtypes))
  } else {
    subtypes <- .selection_discover_subtypes(input_dir, networks, methods)
  }
  if (!length(subtypes)) {
    stop("No C1-C4 subtypes were supplied or discovered.", call. = FALSE)
  }
  subtypes <- sort(subtypes)

  output_dir <- .selection_output_target(base_output_path)
  staging_dir <- tempfile(
    pattern = paste0(".", basename(output_dir), "-"),
    tmpdir = dirname(output_dir)
  )
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create ModuleSelection staging directory.", call. = FALSE)
  }
  on.exit({
    if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  processed_groups <- 0L
  for (network in networks) {
    for (method in methods) {
      for (subtype in subtypes) {
        source_dir <- file.path(input_dir, network, method, subtype)
        node_file <- file.path(
          source_dir,
          paste0("node_Module_", network, "_", method, ".txt")
        )
        edge_file <- file.path(
          source_dir,
          paste0("edges_", network, "_", method, ".txt")
        )
        missing <- c(node_file, edge_file)[!file.exists(c(node_file, edge_file))]
        if (length(missing)) {
          message_text <- paste0(
            "Missing ModuleDivision file(s) for ", network, "/", method, "/", subtype,
            ": ", paste(basename(missing), collapse = ", ")
          )
          if (isTRUE(strict)) stop(message_text, call. = FALSE)
          warning(message_text, call. = FALSE)
          next
        }

        nodes <- .selection_read_required(node_file, c("node", "module"))
        edges <- .selection_read_required(edge_file, c("node1", "node2", "module"))
        nodes <- nodes[, c("node", "module"), drop = FALSE]
        edges <- edges[, c("node1", "node2", "module"), drop = FALSE]
        nodes$node <- trimws(as.character(nodes$node))
        nodes$module <- .selection_source_modules(nodes$module)
        edges$node1 <- trimws(as.character(edges$node1))
        edges$node2 <- trimws(as.character(edges$node2))
        edges$module <- .selection_source_modules(
          edges$module,
          allow_zero = TRUE,
          field = "edge module"
        )
        if (any(!nzchar(nodes$node)) || anyDuplicated(nodes$node)) {
          stop("ModuleDivision node identifiers must be non-empty and unique: ", node_file, call. = FALSE)
        }
        if (any(!nzchar(edges$node1)) || any(!nzchar(edges$node2))) {
          stop("ModuleDivision edge endpoints must be non-empty: ", edge_file, call. = FALSE)
        }

        counts <- as.data.frame(table(nodes$module), stringsAsFactors = FALSE)
        names(counts) <- c("module", "count")
        counts$count <- as.integer(counts$count)
        counts <- counts[order(-counts$count, counts$module), , drop = FALSE]
        selected <- counts[counts$count > numberCutoff, , drop = FALSE]
        selected_nodes <- nodes[nodes$module %in% selected$module, , drop = FALSE]
        selected_edges <- edges[
          edges$node1 %in% selected_nodes$node & edges$node2 %in% selected_nodes$node,
          ,
          drop = FALSE
        ]

        group_output <- file.path(staging_dir, network, method, subtype)
        if (!dir.create(group_output, recursive = TRUE, showWarnings = FALSE)) {
          stop("Could not create ModuleSelection group directory: ", group_output, call. = FALSE)
        }
        utils::write.csv(
          counts,
          file.path(group_output, paste0("Module_Number_", network, "_", subtype, ".csv")),
          row.names = FALSE,
          quote = FALSE
        )
        utils::write.csv(
          selected,
          file.path(group_output, paste0("Module_select_", network, "_", subtype, ".csv")),
          row.names = FALSE,
          quote = FALSE
        )
        .manifest_write_tsv(
          selected_nodes,
          file.path(group_output, paste0("node_Module_select_", network, "_", subtype, ".txt"))
        )
        .manifest_write_tsv(
          selected_edges,
          file.path(group_output, paste0("edges_select_", network, "_", subtype, ".txt"))
        )
        if (isTRUE(write_plots)) {
          .selection_plot_counts(
            counts,
            selected$module,
            file.path(group_output, paste0("Module_Number_", network, "_", subtype, ".tiff"))
          )
        }
        processed_groups <- processed_groups + 1L
      }
    }
  }
  if (!processed_groups) stop("No ModuleDivision groups were processed.", call. = FALSE)

  manifest <- build_module_manifest(
    module_selection_dir = staging_dir,
    output_dir = file.path(staging_dir, "standardized"),
    min_size_exclusive = numberCutoff
  )
  if (!file.rename(staging_dir, output_dir)) {
    stop("Could not finalize ModuleSelection output: ", output_dir, call. = FALSE)
  }

  standardized_dir <- file.path(output_dir, "standardized")
  attr(manifest, "output_dir") <- standardized_dir
  attr(manifest, "manifest_file") <- file.path(standardized_dir, "module_manifest.tsv")
  attr(manifest, "qc_file") <- file.path(standardized_dir, "module_manifest_qc.tsv")
  manifest
}
