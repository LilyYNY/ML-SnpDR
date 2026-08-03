#' Plan or run the ML-SnpDR pipeline
#'
#' @param config Path to a YAML configuration.
#' @param defaults Optional defaults YAML file.
#' @param stages Optional stage names to show or run.
#' @param dry_run When `TRUE`, validate configuration and return the stage plan.
#' @return A list containing configuration, plan and status for a dry run.
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
    return(list(config = parsed, plan = plan, status = "planned"))
  }

  stop(
    "Execution is intentionally disabled in bootstrap version 0.0.1. ",
    "Implement and validate each scientific stage before enabling non-dry runs.",
    call. = FALSE
  )
}

