# Explicit recipe-proxy builders for the course scenarios.
# These scale declared formulation inventories only; they do not read
# exchange measurements, direction summaries, or spent-medium data.
build_rpmi_recipe_proxy <- function(keibler_medium,
                                    formulations = utils::read.csv("config/medium_reference_formulations.csv", stringsAsFactors = FALSE),
                                    formulation_id = "Jain2012_NCI60_RPMI",
                                    keibler_glucose_mM = 10) {
  if (!is.list(keibler_medium) || !all(c("id", "source", "bounds") %in% names(keibler_medium))) {
    stop("keibler_medium must contain id, source, and bounds.", call. = FALSE)
  }
  bounds <- keibler_medium$bounds
  required <- c("reaction_id", "lower_bound", "upper_bound", "rationale")
  if (!is.data.frame(bounds) || !all(required %in% names(bounds)) || anyDuplicated(bounds$reaction_id)) {
    stop("Keibler medium bounds have an invalid schema.", call. = FALSE)
  }
  glucose_row <- match("R_EX_glc__D_e", bounds$reaction_id)
  if (is.na(glucose_row) || bounds$upper_bound[[glucose_row]] != 0 || bounds$lower_bound[[glucose_row]] >= 0 ||
      !is.finite(keibler_glucose_mM) || keibler_glucose_mM <= 0) {
    stop("Recipe proxy requires a one-sided Keibler glucose uptake ceiling.", call. = FALSE)
  }
  recipe <- formulations[formulations$formulation_id == formulation_id, , drop = FALSE]
  if (!nrow(recipe) || anyDuplicated(recipe$reaction_id) || !all(recipe$reaction_id %in% bounds$reaction_id) ||
      any(!is.finite(recipe$concentration_mM)) || any(recipe$concentration_mM < 0)) {
    stop("Formulation rows must map uniquely to listed exchanges with finite, non-negative mM values.", call. = FALSE)
  }
  # Same volume and assumed biomass-time as the Exercise 1 glucose cap (mmol/gDW/h per mM).
  inventory_factor <- -bounds$lower_bound[[glucose_row]] / keibler_glucose_mM

  out <- keibler_medium
  out$id <- "Jain_NCI60_RPMI_recipe_proxy_v2"
  out$source <- paste(
    "Jain NCI-60 complete RPMI-1640 + 2 mM glutamine + 5% FBS; standard RPMI-1640 reference formulation (Jain catalog unverified);",
    "every defined component capped with the Exercise 1 Keibler inventory factor (shared assumed biomass-time denominator);",
    "O2/CO2/H2O/H+ kept at the shared exercise setting; FBS not represented; Jain measured exchange outcomes held out"
  )
  rows <- match(recipe$reaction_id, out$bounds$reaction_id)
  out$bounds$lower_bound[rows] <- -recipe$concentration_mM * inventory_factor
  out$bounds$upper_bound[rows] <- 1000
  out$bounds$rationale[rows] <- ifelse(recipe$concentration_mM == 0,
    paste(recipe$component, "absent from formulation; uptake closed; secretion allowed"),
    paste0(recipe$component, " ", signif(recipe$concentration_mM, 4),
      " mM x shared inventory factor; complete-depletion ceiling; secretion allowed"))
  one_way <- out$bounds$reaction_id %in% c("R_EX_glc__D_e", "R_EX_gln__L_e")
  out$bounds$upper_bound[one_way] <- 0
  out$assumptions <- list(
    formulation_id = formulation_id,
    inventory_factor_mmol_per_gDW_h_per_mM = inventory_factor,
    shared_biomass_time_denominator = TRUE,
    exact_rpmi_catalog_resolved = FALSE,
    mcf7_biomass_time_resolved = FALSE,
    fbs_represented = FALSE,
    validation_measurements_used = FALSE
  )
  out
}
