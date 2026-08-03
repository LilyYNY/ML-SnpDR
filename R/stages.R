#' Return the ML-SnpDR stage registry
#'
#' @return A data frame describing the planned execution order and boundaries.
#' @export
mlsnpdr_stage_registry <- function() {
  data.frame(
    stage = sprintf("%02d", 1:15),
    name = c(
      "diff_expression", "network_construction", "module_detection",
      "module_prefilter", "module_manifest", "module_annotation",
      "module_features", "ml_scoring", "ml_topk_pool",
      "survival_and_drug_response", "candidate_selection",
      "sequence_smiles", "binding_affinity", "perturbation_score",
      "candidate_integration"
    ),
    scope = c(
      rep("cohort", 3), rep("all_prefiltered_modules", 5),
      "ml_ranked_modules", "configured_scope", "evidence_join",
      rep("selected_modules_only", 4)
    ),
    implemented = c(
      rep(FALSE, 4), TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE
    ),
    stringsAsFactors = FALSE
  )
}

