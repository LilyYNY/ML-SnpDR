#' Return the ML-SnpDR stage registry
#'
#' ML-SnpDR preserves subnetDR steps 1-9 and inserts the machine-learning
#' selection chain between the original steps 6 and 7. Every primary output is
#' the default primary input of the next dependent step.
#'
#' @return A data frame describing execution order, analysis scope, primary
#'   input/output contracts and implementation status.
#' @export
mlsnpdr_stage_registry <- function() {
  data.frame(
    stage = c("01", "02", "03", "04", "05", "06", "06A", "06B", "06C", "07", "08", "09"),
    name = c(
      "deps",
      "network_construction",
      "module_division",
      "module_selection",
      "module_annotation",
      "drug_response",
      "module_features",
      "ml_scoring",
      "module_triage",
      "sequence_smiles",
      "binding_score",
      "perturbation_score"
    ),
    scope = c(
      "cohort",
      "subtype",
      "subtype_network_method",
      rep("all_selected_by_size_modules", 5),
      "ml_top10_per_subtype",
      rep("one_best_module_per_subtype", 3)
    ),
    primary_input = c(
      "expression + subtype phenotype",
      "differential_expression.tsv",
      "network_manifest.tsv",
      "ModuleDivision node/edge tables",
      "module_manifest.tsv",
      "module_manifest.tsv + expression + drug training panels",
      "module_manifest.tsv + omics/network feature sources",
      "module_features.tsv + subtype labels/model",
      "ml_top10.tsv + survival + drug_response_summary.tsv",
      "selected_modules.tsv + selected DRN files",
      "seq_smiles_manifest.tsv + selected DRN files",
      "binding_scores.tsv + selected module edges"
    ),
    primary_output = c(
      "differential_expression.tsv",
      "network_manifest.tsv",
      "module_division_manifest.tsv",
      "module_manifest.tsv",
      "module_annotation.tsv",
      "drug_response_summary.tsv + per-module DRN files",
      "module_features.tsv",
      "ml_scores.tsv + ml_top10.tsv",
      "module_evidence.tsv + selected_modules.tsv",
      "seq_smiles_manifest.tsv",
      "binding_scores.tsv",
      "perturbation_scores.tsv + final_candidates.tsv"
    ),
    implemented = rep(TRUE, 12L),
    stringsAsFactors = FALSE
  )
}
