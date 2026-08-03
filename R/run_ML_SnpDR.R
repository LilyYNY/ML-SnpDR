#' Run the complete ML-SnpDR downstream chain
#'
#' This convenience runner mirrors subnetDR's `run_subnetDR.R` while exposing
#' explicit start/end stages. The default chain starts from ModuleDivision
#' outputs at step 4, runs annotation and all-module drug response/features,
#' creates ML Top-K and one best module per subtype, then runs steps 7-9 only
#' for those selected modules.
#'
#' @param config Path to an ML-SnpDR YAML configuration.
#' @param defaults Optional defaults YAML merged before `config`.
#' @param from First stage name; defaults to `module_selection`.
#' @param to Last stage name; defaults to `perturbation_score`.
#' @param dry_run Validate and show the selected chain without scientific
#'   computation.
#' @return The result of [run_ml_snpdr_pipeline()].
#' @export
run_ML_SnpDR <- function(
  config,
  defaults = NULL,
  from = "module_selection",
  to = "perturbation_score",
  dry_run = FALSE
) {
  registry <- mlsnpdr_stage_registry()
  from_index <- match(from, registry$name)
  to_index <- match(to, registry$name)
  if (is.na(from_index)) stop("Unknown from stage: ", from, call. = FALSE)
  if (is.na(to_index)) stop("Unknown to stage: ", to, call. = FALSE)
  if (from_index > to_index) stop("from stage must not occur after to stage.", call. = FALSE)
  selected <- registry$name[seq.int(from_index, to_index)]
  unavailable <- selected[!registry$implemented[match(selected, registry$name)]]
  if (length(unavailable) && !isTRUE(dry_run)) {
    stop(
      "The requested range includes stages without an ML-SnpDR implementation: ",
      paste(unavailable, collapse = ", "),
      call. = FALSE
    )
  }
  run_ml_snpdr_pipeline(
    config = config,
    defaults = defaults,
    stages = selected,
    dry_run = dry_run
  )
}
