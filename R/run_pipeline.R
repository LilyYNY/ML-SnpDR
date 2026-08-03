#' Plan or run the ML-SnpDR chained pipeline
#'
#' @param config Path to a YAML configuration.
#' @param defaults Optional defaults YAML file.
#' @param stages Optional stage names to show or run. ML-SnpDR keeps subnetDR
#'   steps 1-9 and inserts `module_features`, `ml_scoring` and `module_triage`
#'   between the original steps 6 and 7.
#' @param dry_run When `TRUE`, validate configuration and return the stage plan.
#' @return A list containing configuration, plan, status and stage outputs.
#' @export
run_ml_snpdr_pipeline <- function(config, defaults = NULL, stages = NULL, dry_run = TRUE) {
  parsed <- read_mlsnpdr_config(config, defaults = defaults)
  plan <- mlsnpdr_stage_registry()
  if (!is.null(stages)) {
    unknown <- setdiff(stages, plan$name)
    if (length(unknown)) stop("Unknown stage(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    plan <- plan[plan$name %in% stages, , drop = FALSE]
  }

  if (isTRUE(dry_run)) {
    message("ML-SnpDR configuration is valid. Returning a dry-run plan.")
    return(list(config = parsed, plan = plan, status = "planned", outputs = list()))
  }

  unavailable <- plan$name[!plan$implemented]
  if (length(unavailable)) {
    stop(
      "Cannot execute because these requested stages are not implemented yet: ",
      paste(unavailable, collapse = ", "),
      ". Run dry_run = TRUE to inspect the complete chain.",
      call. = FALSE
    )
  }

  outputs <- list()
  if ("module_selection" %in% plan$name) {
    required_paths <- c("module_division_dir", "module_selection_dir")
    missing_paths <- required_paths[
      !vapply(required_paths, function(x) {
        value <- parsed$paths[[x]]
        !is.null(value) && length(value) == 1L && nzchar(value)
      }, logical(1))
    ]
    if (length(missing_paths)) {
      stop(
        "Missing path setting(s) required for module_selection: ",
        paste(missing_paths, collapse = ", "),
        call. = FALSE
      )
    }
    outputs$module_selection <- module_selection(
      subtype_file = NULL,
      base_input_path = parsed$paths$module_division_dir,
      base_output_path = parsed$paths$module_selection_dir,
      network_method = parsed$pipeline$network_methods,
      module_method = parsed$pipeline$module_methods,
      numberCutoff = parsed$module_prefilter$min_size_exclusive,
      subtypes = parsed$project$target_subtypes,
      strict = parsed$pipeline$strict_inputs,
      write_plots = parsed$pipeline$write_module_plots
    )
  }

  if ("module_annotation" %in% plan$name) {
    annotation_output <- parsed$paths$module_annotation_dir
    if (is.null(annotation_output) || length(annotation_output) != 1L ||
        is.na(annotation_output) || !nzchar(annotation_output)) {
      stop("Missing path setting required for module_annotation: module_annotation_dir", call. = FALSE)
    }
    manifest_file <- if (!is.null(outputs$module_selection)) {
      attr(outputs$module_selection, "manifest_file")
    } else {
      .annotation_manifest_path(NULL, parsed$paths$module_selection_dir)
    }
    outputs$module_annotation <- functional_annotation(
      subtype_file = NULL,
      base_input_path = parsed$paths$module_selection_dir,
      base_output_path = annotation_output,
      network_method = parsed$pipeline$network_methods,
      module_method = parsed$pipeline$module_methods,
      module_manifest_file = manifest_file,
      gene_set_file = parsed$paths$annotation_gene_set_file %||% NULL,
      background_gene_file = parsed$paths$annotation_background_gene_file %||% NULL,
      databases = parsed$annotation$databases,
      species = parsed$annotation$species,
      subtypes = parsed$project$target_subtypes,
      pvalueCutoff = parsed$annotation$pvalue_cutoff,
      pAdjustMethod = parsed$annotation$p_adjust_method,
      pAdjustCutoff = parsed$annotation$p_adjust_cutoff,
      qvalueCutoff = parsed$annotation$qvalue_cutoff,
      minGSSize = parsed$annotation$min_gene_set_size,
      maxGSSize = parsed$annotation$max_gene_set_size,
      top_n = parsed$annotation$top_n,
      strict = parsed$pipeline$strict_inputs,
      write_module_files = parsed$annotation$write_module_files
    )
  }

  if ("drug_response" %in% plan$name) {
    manifest_file <- if (!is.null(outputs$module_selection)) {
      attr(outputs$module_selection, "manifest_file")
    } else {
      .annotation_manifest_path(NULL, parsed$paths$module_selection_dir)
    }
    input_index <- parsed$paths$drug_response_input_index_file %||% NULL
    input_roots <- parsed$drug_response$input_roots %||% NULL
    if (!length(input_roots)) input_roots <- NULL
    if (!xor(is.null(input_index), is.null(input_roots))) {
      stop(
        "Configure exactly one of paths.drug_response_input_index_file or ",
        "drug_response.input_roots before running drug_response.",
        call. = FALSE
      )
    }
    outputs$drug_response <- drug_response_analysis(
      module_manifest_file = manifest_file,
      drug_response_path = parsed$paths$drug_response_dir,
      drug_response_roots = input_roots,
      drug_response_index_file = input_index,
      panels = parsed$drug_response$panels,
      network_methods = parsed$pipeline$network_methods,
      module_methods = parsed$pipeline$module_methods,
      subtypes = parsed$project$target_subtypes,
      significance_cutoff = parsed$drug_response$significance_cutoff,
      p_adjust_method = parsed$drug_response$p_adjust_method,
      strict = parsed$pipeline$strict_inputs
    )
  }

  if ("module_features" %in% plan$name) {
    manifest_file <- if (!is.null(outputs$module_selection)) {
      attr(outputs$module_selection, "manifest_file")
    } else {
      .annotation_manifest_path(NULL, parsed$paths$module_selection_dir)
    }
    feature_input <- parsed$paths$module_feature_input_file %||% NULL
    if (is.null(feature_input) || length(feature_input) != 1L ||
        is.na(feature_input) || !nzchar(feature_input)) {
      stop("Configure paths.module_feature_input_file before running module_features.", call. = FALSE)
    }
    outputs$module_features <- prepare_module_features(
      module_manifest_file = manifest_file,
      feature_file = feature_input,
      output_dir = parsed$paths$module_features_dir,
      feature_set = parsed$ml$feature_set,
      expected_feature_count = parsed$ml$expected_feature_count,
      strict = parsed$pipeline$strict_inputs
    )
  }

  if ("ml_scoring" %in% plan$name) {
    module_features_file <- if (!is.null(outputs$module_features)) {
      attr(outputs$module_features, "feature_file")
    } else {
      file.path(parsed$paths$module_features_dir, "module_features.tsv")
    }
    score_input <- parsed$paths$ml_score_input_file %||% NULL
    outputs$ml_scoring <- if (!is.null(score_input) && length(score_input) == 1L &&
                              !is.na(score_input) && nzchar(score_input)) {
      prepare_ml_scores(
        module_features_file = module_features_file,
        score_file = score_input,
        output_dir = parsed$paths$ml_scoring_dir,
        top_k = parsed$ml$top_k_per_subtype,
        require_predicted_subtype_match = parsed$ml$require_predicted_subtype_match,
        probability_min = parsed$ml$probability_min,
        margin_min = parsed$ml$margin_min,
        model_version = parsed$ml$model_version,
        strict = parsed$pipeline$strict_inputs
      )
    } else {
      feature_schema_file <- if (!is.null(outputs$module_features)) {
        attr(outputs$module_features, "schema_file")
      } else {
        file.path(parsed$paths$module_features_dir, "feature_schema.json")
      }
      run_nested_ml_scoring(
        module_features_file = module_features_file,
        feature_schema_file = feature_schema_file,
        output_dir = parsed$paths$ml_scoring_dir,
        python_executable = parsed$ml$python_executable,
        mode = if (parsed$project$mode == "paper_reproduction") "paper" else "fast",
        seed = parsed$project$seed,
        outer_splits = parsed$ml$outer_splits,
        outer_repeats = parsed$ml$outer_repeats,
        inner_splits = parsed$ml$inner_splits,
        top_k = parsed$ml$top_k_per_subtype,
        n_jobs = parsed$ml$n_jobs,
        require_predicted_subtype_match = parsed$ml$require_predicted_subtype_match
      )
    }
  }

  if ("module_triage" %in% plan$name) {
    manifest_file <- if (!is.null(outputs$module_selection)) {
      attr(outputs$module_selection, "manifest_file")
    } else {
      .annotation_manifest_path(NULL, parsed$paths$module_selection_dir)
    }
    ml_top_file <- if (!is.null(outputs$ml_scoring)) {
      attr(outputs$ml_scoring, "top_file")
    } else {
      file.path(
        parsed$paths$ml_scoring_dir,
        paste0("ml_top", parsed$ml$top_k_per_subtype, ".tsv")
      )
    }
    drug_summary_file <- if (!is.null(outputs$drug_response)) {
      attr(outputs$drug_response, "summary_file")
    } else {
      file.path(parsed$paths$drug_response_dir, "drug_response_summary.tsv")
    }
    survival_input <- parsed$paths$survival_input_file %||% NULL
    if (is.null(survival_input) || length(survival_input) != 1L ||
        is.na(survival_input) || !nzchar(survival_input)) {
      stop("Configure paths.survival_input_file before running module_triage.", call. = FALSE)
    }
    outputs$module_triage <- triage_modules(
      module_manifest_file = manifest_file,
      ml_top_file = ml_top_file,
      survival_file = survival_input,
      drug_response_summary_file = drug_summary_file,
      output_dir = parsed$paths$module_triage_dir,
      primary_drug_panel = parsed$candidate_selection$primary_drug_panel,
      min_module_size = parsed$candidate_selection$min_module_size,
      required_direction = parsed$survival$required_direction,
      logrank_p_max = parsed$survival$logrank_p_max,
      logrank_fdr_max = parsed$survival$logrank_fdr_max,
      sort_by = parsed$candidate_selection$sort_by,
      select_n_per_subtype = parsed$candidate_selection$select_n_per_subtype,
      strict = parsed$pipeline$strict_inputs
    )
  }

  if ("sequence_smiles" %in% plan$name) {
    selected_file <- if (!is.null(outputs$module_triage)) {
      attr(outputs$module_triage, "selected_file")
    } else {
      file.path(parsed$paths$module_triage_dir, "selected_modules.tsv")
    }
    protein_file <- parsed$paths$protein_sequence_file %||% NULL
    smiles_file <- parsed$paths$drug_smiles_file %||% NULL
    if (is.null(protein_file) || is.null(smiles_file)) {
      stop(
        "Configure paths.protein_sequence_file and paths.drug_smiles_file ",
        "before running sequence_smiles.",
        call. = FALSE
      )
    }
    outputs$sequence_smiles <- run_SEQCre(
      input_base = selected_file,
      output_base = parsed$paths$sequence_smiles_dir,
      protein_sequence_file = protein_file,
      drug_smiles_file = smiles_file,
      strict = parsed$pipeline$strict_inputs
    )
  }

  if ("binding_score" %in% plan$name) {
    selected_file <- if (!is.null(outputs$module_triage)) {
      attr(outputs$module_triage, "selected_file")
    } else {
      file.path(parsed$paths$module_triage_dir, "selected_modules.tsv")
    }
    seq_manifest_file <- if (!is.null(outputs$sequence_smiles)) {
      attr(outputs$sequence_smiles, "manifest_file")
    } else {
      file.path(parsed$paths$sequence_smiles_dir, "seq_smiles_manifest.tsv")
    }
    binding_input <- parsed$paths$binding_score_input_file %||% NULL
    if (is.null(binding_input)) {
      stop(
        "Configure paths.binding_score_input_file before running binding_score. ",
        "Use DeepPurpose output with module_uid/target/drug/binding_score columns.",
        call. = FALSE
      )
    }
    outputs$binding_score <- predict_BA(
      selected_modules_file = selected_file,
      seq_smiles_manifest_file = seq_manifest_file,
      output_base = parsed$paths$binding_score_dir,
      binding_score_file = binding_input,
      score_direction = parsed$post_selection$binding_score_direction,
      strict = parsed$pipeline$strict_inputs
    )
  }

  if ("perturbation_score" %in% plan$name) {
    selected_file <- if (!is.null(outputs$module_triage)) {
      attr(outputs$module_triage, "selected_file")
    } else {
      file.path(parsed$paths$module_triage_dir, "selected_modules.tsv")
    }
    binding_file <- if (!is.null(outputs$binding_score)) {
      attr(outputs$binding_score, "scores_file")
    } else {
      file.path(parsed$paths$binding_score_dir, "binding_scores.tsv")
    }
    sensitivity_input <- parsed$paths$sensitivity_input_file %||% NULL
    if (is.null(sensitivity_input)) {
      stop(
        "Configure paths.sensitivity_input_file before running perturbation_score. ",
        "Use ENM/PRS output with module_uid/target/sensitivity columns.",
        call. = FALSE
      )
    }
    outputs$perturbation_score <- process_prs_dti(
      selected_modules_file = selected_file,
      binding_scores_file = binding_file,
      output_base = parsed$paths$perturbation_score_dir,
      sensitivity_file = sensitivity_input,
      top_n = parsed$post_selection$final_top_n_per_subtype,
      strict = parsed$pipeline$strict_inputs
    )
  }

  list(config = parsed, plan = plan, status = "completed", outputs = outputs)
}
