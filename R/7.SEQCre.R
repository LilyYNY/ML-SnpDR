# Sequence and SMILES --------------------------------------------------------

.selected_read <- function(path, verify_files = TRUE) {
  selected <- .annotation_read_tabular(path)
  required <- c(
    "module_uid", "network", "method", "subtype", "module", "module_size",
    "node_file", "edge_file", "drn_file", "drn_info_file", "selection_rank"
  )
  missing <- setdiff(required, names(selected))
  if (length(missing)) stop("selected_modules.tsv is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyDuplicated(selected$module_uid)) stop("selected_modules.tsv contains duplicate module_uid values.", call. = FALSE)
  identities <- parse_module_uid(selected$module_uid)
  for (field in c("network", "method", "subtype", "module")) {
    if (any(as.character(selected[[field]]) != as.character(identities[[field]]))) {
      stop("selected_modules.tsv identity mismatch: ", field, call. = FALSE)
    }
  }
  root <- dirname(normalizePath(path, winslash = "/", mustWork = TRUE))
  for (column in c("node_file", "edge_file", "drn_file", "drn_info_file")) {
    selected[[paste0(column, "_abs")]] <- vapply(
      selected[[column]],
      .drug_absolute_path,
      character(1),
      root = root
    )
  }
  if (isTRUE(verify_files)) {
    absolute_columns <- paste0(c("node_file", "edge_file", "drn_file", "drn_info_file"), "_abs")
    missing_paths <- unlist(lapply(selected[, absolute_columns, drop = FALSE], function(paths) paths[!file.exists(paths)]))
    if (length(missing_paths)) {
      stop("selected_modules.tsv references missing handoff file(s): ", missing_paths[[1L]], call. = FALSE)
    }
  }
  attr(selected, "selected_file") <- normalizePath(path, winslash = "/", mustWork = TRUE)
  attr(selected, "selected_root") <- root
  selected
}

.seq_lookup <- function(path, type) {
  data <- .annotation_read_tabular(path)
  lowered <- tolower(names(data))
  if (type == "protein") {
    id_aliases <- c("node", "protein", "target", "gene", "gene_symbol")
    value_aliases <- c("sequence", "seq", "protein_sequence")
  } else {
    id_aliases <- c("node", "drug", "drug_name", "compound")
    value_aliases <- c("smiles", "canonical_smiles", "drug_smiles")
  }
  id_hit <- match(id_aliases, lowered, nomatch = 0L)
  value_hit <- match(value_aliases, lowered, nomatch = 0L)
  id_hit <- id_hit[id_hit > 0L]
  value_hit <- value_hit[value_hit > 0L]
  if (!length(id_hit) || !length(value_hit)) {
    stop(type, " lookup must contain an identifier and sequence/SMILES column.", call. = FALSE)
  }
  result <- data.frame(
    node = trimws(as.character(data[[id_hit[[1L]]]])),
    value = trimws(as.character(data[[value_hit[[1L]]]])),
    stringsAsFactors = FALSE
  )
  result <- result[nzchar(result$node) & nzchar(result$value), , drop = FALSE]
  if (anyDuplicated(result$node)) {
    conflicts <- vapply(split(result$value, result$node), function(x) length(unique(x)) > 1L, logical(1))
    if (any(conflicts)) stop(type, " lookup contains conflicting duplicate identifiers.", call. = FALSE)
    result <- result[!duplicated(result$node), , drop = FALSE]
  }
  result
}

#' Create protein-sequence and drug-SMILES files for selected modules only
#'
#' Step 7 reads `selected_modules.tsv` and its copied DRN info files. Protein
#' sequences and drug SMILES are joined from explicit lookup tables, and one
#' pair of files is written per selected module. No directory scanning occurs.
#'
#' @param input_base Path to step-6C `selected_modules.tsv` (kept as the first
#'   argument for compatibility with subnetDR's `run_SEQCre()`).
#' @param output_base New step-7 output directory.
#' @param protein_sequence_file Tabular lookup with protein/node and sequence.
#' @param drug_smiles_file Tabular lookup with drug/node and SMILES.
#' @param strict Require complete sequence and SMILES coverage.
#' @return `seq_smiles_manifest.tsv` with paths in attributes.
#' @export
run_SEQCre <- function(
  input_base,
  output_base,
  protein_sequence_file,
  drug_smiles_file,
  strict = TRUE
) {
  selected <- .selected_read(input_base, verify_files = TRUE)
  protein_lookup <- .seq_lookup(protein_sequence_file, "protein")
  drug_lookup <- .seq_lookup(drug_smiles_file, "drug")
  output_target <- .annotation_output_target(output_base)
  staging_dir <- tempfile(pattern = paste0(".", basename(output_target), "-"), tmpdir = dirname(output_target))
  if (!dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create sequence/SMILES staging directory.", call. = FALSE)
  }
  on.exit({
    if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  rows <- vector("list", nrow(selected))
  for (index in seq_len(nrow(selected))) {
    module <- selected[index, , drop = FALSE]
    info <- .drug_read_info(module$drn_info_file_abs)
    proteins <- sort(unique(info$node[info$type == "protein"]))
    drugs <- sort(unique(info$node[info$type == "drug"]))
    sequence <- protein_lookup[match(proteins, protein_lookup$node), , drop = FALSE]
    smiles <- drug_lookup[match(drugs, drug_lookup$node), , drop = FALSE]
    missing_proteins <- proteins[is.na(sequence$node)]
    missing_drugs <- drugs[is.na(smiles$node)]
    if (isTRUE(strict) && (length(missing_proteins) || length(missing_drugs))) {
      stop(
        "Sequence/SMILES coverage failed for ", module$module_uid,
        ": missing proteins=", length(missing_proteins),
        ", missing drugs=", length(missing_drugs),
        call. = FALSE
      )
    }
    sequence <- data.frame(
      node = proteins,
      sequence = sequence$value,
      sequence_status = ifelse(is.na(sequence$value), "missing", "available"),
      stringsAsFactors = FALSE
    )
    smiles <- data.frame(
      node = drugs,
      SMILES = smiles$value,
      smiles_status = ifelse(is.na(smiles$value), "missing", "available"),
      stringsAsFactors = FALSE
    )
    relative_dir <- file.path("selected", module$subtype, module$module_uid)
    module_dir <- file.path(staging_dir, relative_dir)
    if (!dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create sequence/SMILES module directory.", call. = FALSE)
    }
    seq_relative <- gsub("\\\\", "/", file.path(relative_dir, "protein_sequences.tsv"))
    smiles_relative <- gsub("\\\\", "/", file.path(relative_dir, "drug_smiles.tsv"))
    .manifest_write_tsv(sequence, file.path(staging_dir, seq_relative))
    .manifest_write_tsv(smiles, file.path(staging_dir, smiles_relative))
    complete <- !length(missing_proteins) && !length(missing_drugs)
    rows[[index]] <- data.frame(
      module_uid = module$module_uid,
      subtype = module$subtype,
      sequence_file = seq_relative,
      smiles_file = smiles_relative,
      required_protein_number = length(proteins),
      matched_protein_number = sum(sequence$sequence_status == "available"),
      required_drug_number = length(drugs),
      matched_drug_number = sum(smiles$smiles_status == "available"),
      analysis_status = if (complete) "completed" else "completed_with_missing_lookup",
      status_reason = if (complete) "" else paste0(
        "missing_proteins=", length(missing_proteins),
        ";missing_drugs=", length(missing_drugs)
      ),
      stringsAsFactors = FALSE
    )
  }
  manifest <- do.call(rbind, rows)
  rownames(manifest) <- NULL
  manifest_file <- file.path(staging_dir, "seq_smiles_manifest.tsv")
  .manifest_write_tsv(manifest, manifest_file)
  if (!file.rename(staging_dir, output_target)) {
    stop("Could not finalize sequence/SMILES output directory: ", output_target, call. = FALSE)
  }
  attr(manifest, "output_dir") <- output_target
  attr(manifest, "manifest_file") <- file.path(output_target, basename(manifest_file))
  manifest
}
