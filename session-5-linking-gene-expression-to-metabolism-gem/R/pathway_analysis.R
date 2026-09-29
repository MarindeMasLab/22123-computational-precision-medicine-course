summarize_pathway_metrics <- function(reaction_fluxes,
                                     excluded_subsystems = c("Transport reactions", "Exchange/demand reactions", "Biomass", "Unassigned"),
                                     carrying_flux_tolerance = 1e-6) {
  required <- c("scenario", "reaction_id", "subsystem", "mapping_source", "all_human1_subsystems", "pfba_flux",
    "limited_by_expression", "biomass_objective", "medium_id", "objective_id", "objective_sense",
    "model_sha256", "pfba_status", "max_abs_sv")
  if (!is.data.frame(reaction_fluxes) || !all(required %in% names(reaction_fluxes))) {
    stop("Reaction pathway data are missing required columns.", call. = FALSE)
  }
  if (anyNA(reaction_fluxes[, required]) || anyDuplicated(reaction_fluxes[, c("scenario", "reaction_id")])) {
    stop("Reaction pathway data must be complete and unique by scenario/reaction.", call. = FALSE)
  }
  if (length(carrying_flux_tolerance) != 1L || !is.finite(carrying_flux_tolerance) || carrying_flux_tolerance < 0) {
    stop("carrying_flux_tolerance must be one non-negative finite number.", call. = FALSE)
  }
  if (any(!is.finite(reaction_fluxes$pfba_flux)) || any(!is.finite(reaction_fluxes$biomass_objective)) ||
      any(reaction_fluxes$biomass_objective <= 0)) {
    stop("pFBA fluxes must be finite and biomass objectives positive.", call. = FALSE)
  }
  if (any(reaction_fluxes$pfba_status != "optimal") || any(!is.finite(reaction_fluxes$max_abs_sv)) ||
      any(reaction_fluxes$max_abs_sv > 1e-7)) stop("All pFBA solves must be optimal and mass-balanced.", call. = FALSE)
  if (anyNA(reaction_fluxes$limited_by_expression) || any(!reaction_fluxes$limited_by_expression %in% c(TRUE, FALSE))) {
    stop("limited_by_expression must be logical and non-missing.", call. = FALSE)
  }

  scenario_ids <- unique(reaction_fluxes$scenario)
  if (length(scenario_ids) < 1L) stop("No scenarios supplied.", call. = FALSE)
  reference_ids <- sort(reaction_fluxes$reaction_id[reaction_fluxes$scenario == scenario_ids[[1L]]])
  reference_map <- reaction_fluxes[reaction_fluxes$scenario == scenario_ids[[1L]], c("reaction_id", "subsystem")]
  reference_map <- reference_map[order(reference_map$reaction_id), ]
  for (scenario_id in scenario_ids) {
    current <- reaction_fluxes[reaction_fluxes$scenario == scenario_id, ]
    context_fields <- c("medium_id", "objective_id", "objective_sense", "model_sha256", "biomass_objective", "pfba_status")
    if (any(vapply(current[context_fields], function(x) length(unique(x)) != 1L, logical(1)))) {
      stop("Context provenance/objective must be constant within each scenario.", call. = FALSE)
    }
    if (!scenario_id %in% scenario_ids[-1L]) next
    if (!identical(sort(current$reaction_id), reference_ids)) {
      stop("Scenarios must contain identical reaction ID sets.", call. = FALSE)
    }
    current_map <- current[match(reference_map$reaction_id, current$reaction_id), c("reaction_id", "subsystem")]
    if (!identical(as.character(current_map$subsystem), as.character(reference_map$subsystem))) {
      stop("Subsystem assignments differ between scenarios.", call. = FALSE)
    }
  }

  reaction_fluxes$pathway_class <- ifelse(reaction_fluxes$subsystem %in% excluded_subsystems,
    "coverage/excluded from intracellular pathway weights", "intracellular pathway")
  reaction_fluxes$carrying_flux <- abs(reaction_fluxes$pfba_flux) > carrying_flux_tolerance
  groups <- split(reaction_fluxes,
    list(reaction_fluxes$scenario, reaction_fluxes$subsystem, reaction_fluxes$pathway_class), drop = TRUE)
  summary <- do.call(rbind, lapply(groups, function(x) {
    data.frame(scenario = x$scenario[[1L]], subsystem = x$subsystem[[1L]], pathway_class = x$pathway_class[[1L]],
      mapping_source = paste(sort(unique(x$mapping_source)), collapse = ";"),
      all_human1_subsystems = paste(sort(unique(unlist(strsplit(x$all_human1_subsystems, ";", fixed = TRUE)))), collapse = ";"),
      reaction_count = nrow(x), carrying_flux_count = sum(x$carrying_flux),
      carrying_flux_fraction = mean(x$carrying_flux),
      limited_by_expression_count = sum(x$limited_by_expression),
      limited_by_expression_fraction = mean(x$limited_by_expression),
      abs_flux_sum = sum(abs(x$pfba_flux)),
      pFBA_biomass_objective = x$biomass_objective[[1L]],
      throughput_per_biomass_objective = sum(abs(x$pfba_flux)) / x$biomass_objective[[1L]],
      medium_id = x$medium_id[[1L]], objective_id = x$objective_id[[1L]],
      objective_sense = x$objective_sense[[1L]], model_sha256 = x$model_sha256[[1L]],
      solver_status = x$pfba_status[[1L]], stringsAsFactors = FALSE)
  }))
  internal_totals <- tapply(summary$abs_flux_sum[summary$pathway_class == "intracellular pathway"],
    summary$scenario[summary$pathway_class == "intracellular pathway"], sum)
  summary$intracellular_abs_flux_denominator <- unname(internal_totals[summary$scenario])
  if (any(!is.finite(summary$intracellular_abs_flux_denominator)) || any(summary$intracellular_abs_flux_denominator <= 0)) {
    stop("Each context requires positive internal pathway throughput for flux-share calculation.", call. = FALSE)
  }
  summary$intracellular_flux_share <- ifelse(summary$pathway_class == "intracellular pathway",
    summary$abs_flux_sum / summary$intracellular_abs_flux_denominator, NA_real_)
  summary$carrying_flux_percent <- 100 * summary$carrying_flux_fraction
  summary$limited_by_expression_percent <- 100 * summary$limited_by_expression_fraction
  summary$carrying_flux_tolerance <- carrying_flux_tolerance

  coverage_groups <- split(reaction_fluxes,
    list(reaction_fluxes$scenario, reaction_fluxes$pathway_class), drop = TRUE)
  coverage <- do.call(rbind, lapply(coverage_groups, function(x) data.frame(
    scenario = x$scenario[[1L]], pathway_class = x$pathway_class[[1L]],
    reaction_count = nrow(x), carrying_flux_count = sum(x$carrying_flux),
    limited_by_expression_count = sum(x$limited_by_expression), abs_flux_sum = sum(abs(x$pfba_flux)),
    model_reaction_count = length(reference_ids), stringsAsFactors = FALSE)))
  summary <- summary[order(summary$scenario, summary$pathway_class, summary$subsystem), ]
  rownames(summary) <- NULL
  coverage <- coverage[order(coverage$scenario, coverage$pathway_class), ]
  rownames(coverage) <- NULL
  list(summary = summary, coverage = coverage, reaction_fluxes = reaction_fluxes,
    carrying_flux_tolerance = carrying_flux_tolerance, excluded_subsystems = excluded_subsystems)
}

compare_pathway_metrics <- function(pathway_summary, scenario_a, scenario_b) {
  required <- c("scenario", "subsystem", "pathway_class", "mapping_source", "all_human1_subsystems", "reaction_count", "carrying_flux_count",
    "carrying_flux_fraction", "limited_by_expression_count", "limited_by_expression_fraction",
    "abs_flux_sum", "throughput_per_biomass_objective", "intracellular_flux_share")
  if (!is.data.frame(pathway_summary) || !all(required %in% names(pathway_summary))) {
    stop("Pathway summary is missing required comparison columns.", call. = FALSE)
  }
  a <- pathway_summary[pathway_summary$scenario == scenario_a, , drop = FALSE]
  b <- pathway_summary[pathway_summary$scenario == scenario_b, , drop = FALSE]
  if (!nrow(a) || !nrow(b)) stop("Both requested scenarios must occur in pathway_summary.", call. = FALSE)
  key <- c("subsystem", "pathway_class")
  if (anyDuplicated(a[, key]) || anyDuplicated(b[, key])) stop("Pathway summary has duplicate scenario/subsystem rows.", call. = FALSE)
  a <- a[order(a$subsystem, a$pathway_class), , drop = FALSE]
  b <- b[match(paste(a$subsystem, a$pathway_class), paste(b$subsystem, b$pathway_class)), , drop = FALSE]
  if (anyNA(b$scenario)) stop("Contexts do not have the same subsystem/category universe.", call. = FALSE)
  out <- data.frame(subsystem = a$subsystem, pathway_class = a$pathway_class,
    mapping_source_a = a$mapping_source, mapping_source_b = b$mapping_source,
    all_human1_subsystems_a = a$all_human1_subsystems, all_human1_subsystems_b = b$all_human1_subsystems,
    reaction_count = a$reaction_count,
    abs_flux_sum_a = a$abs_flux_sum, abs_flux_sum_b = b$abs_flux_sum,
    abs_flux_sum_delta_a_minus_b = a$abs_flux_sum - b$abs_flux_sum,
    throughput_per_biomass_a = a$throughput_per_biomass_objective,
    throughput_per_biomass_b = b$throughput_per_biomass_objective,
    throughput_per_biomass_delta_a_minus_b = a$throughput_per_biomass_objective - b$throughput_per_biomass_objective,
    flux_share_a = a$intracellular_flux_share, flux_share_b = b$intracellular_flux_share,
    flux_share_delta_a_minus_b = a$intracellular_flux_share - b$intracellular_flux_share,
    carrying_flux_count_a = a$carrying_flux_count, carrying_flux_count_b = b$carrying_flux_count,
    carrying_flux_fraction_a = a$carrying_flux_fraction, carrying_flux_fraction_b = b$carrying_flux_fraction,
    carrying_flux_delta_percentage_points_a_minus_b = 100 * (a$carrying_flux_fraction - b$carrying_flux_fraction),
    limited_by_expression_count_a = a$limited_by_expression_count,
    limited_by_expression_count_b = b$limited_by_expression_count,
    limited_by_expression_fraction_a = a$limited_by_expression_fraction,
    limited_by_expression_fraction_b = b$limited_by_expression_fraction,
    limited_by_expression_delta_percentage_points_a_minus_b =
      100 * (a$limited_by_expression_fraction - b$limited_by_expression_fraction),
    intracellular_abs_flux_denominator_a = a$intracellular_abs_flux_denominator,
    intracellular_abs_flux_denominator_b = b$intracellular_abs_flux_denominator,
    carrying_flux_tolerance_a = a$carrying_flux_tolerance,
    carrying_flux_tolerance_b = b$carrying_flux_tolerance,
    scenario_a = scenario_a, scenario_b = scenario_b,
    medium_id_a = a$medium_id, medium_id_b = b$medium_id,
    objective_id_a = a$objective_id, objective_id_b = b$objective_id,
    objective_sense_a = a$objective_sense, objective_sense_b = b$objective_sense,
    model_sha256_a = a$model_sha256, model_sha256_b = b$model_sha256,
    solver_status_a = a$solver_status, solver_status_b = b$solver_status,
    stringsAsFactors = FALSE)
  out[order(out$pathway_class, out$subsystem), ]
}
