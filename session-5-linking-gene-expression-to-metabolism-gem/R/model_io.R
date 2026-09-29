# Shared Session 5 SBML parser and baseline FBA API.
# Supported input: Level 3 SBML with FBC v2 bounds/objectives and nested AND/OR GPRs.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[1L])) y else x

.xml_attr_local <- function(node, local_name, required = TRUE) {
  attrs <- xml2::xml_attrs(node)
  hit <- which(sub("^.*:", "", names(attrs)) == local_name)
  if (length(hit) != 1L) {
    if (!required && !length(hit)) return(NA_character_)
    stop("Expected exactly one XML attribute '", local_name, "'.", call. = FALSE)
  }
  unname(attrs[[hit]])
}

.xml_children_local <- function(node, local_name) {
  xml2::xml_find_all(node, paste0("./*[local-name()='", local_name, "']"))
}

.unique_ids <- function(ids, label) {
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(label, " IDs must be present and unique.", call. = FALSE)
  }
  ids
}

.parse_gpr_node <- function(node, gene_ids) {
  kind <- sub("^.*:", "", xml2::xml_name(node))
  if (kind == "geneProductRef") {
    id <- .xml_attr_local(node, "geneProduct")
    if (!id %in% gene_ids) stop("GPR refers to unknown gene product: ", id, call. = FALSE)
    return(id)
  }
  if (!kind %in% c("and", "or")) stop("Unsupported GPR node: ", kind, call. = FALSE)
  children <- xml2::xml_children(node)
  if (length(children) < 2L) stop("FBC AND/OR GPR nodes need at least two children.", call. = FALSE)
  paste0("(", paste(vapply(children, .parse_gpr_node, character(1), gene_ids = gene_ids), collapse = paste0(" ", kind, " ")), ")")
}

.resolve_bound <- function(value, parameters, reaction_id, bound_name) {
  if (value %in% names(parameters)) return(unname(parameters[[value]]))
  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed)) stop("Unsupported or unresolved ", bound_name, " for ", reaction_id, ": ", value, call. = FALSE)
  parsed
}

.read_manifest <- function() {
  manifest <- utils::read.csv("../config/model_manifest.csv", stringsAsFactors = FALSE)
  if (!all(c("model_file", "sha256", "genes", "reactions", "metabolites", "active_objective") %in% names(manifest))) {
    stop("Model manifest has an invalid schema.", call. = FALSE)
  }
  manifest
}

# Load a source SBML without normalizing IDs. Run from exercises/R or pass a
# manifest_path explicitly. Every parse is checked against local source hashes/counts.
load_model <- function(path, manifest_path = "config/model_manifest.csv") {
  if (!requireNamespace("xml2", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE) || !requireNamespace("digest", quietly = TRUE)) {
    stop("load_model() requires xml2, Matrix, and digest.", call. = FALSE)
  }
  if (!file.exists(path)) stop("Model file does not exist: ", path, call. = FALSE)
  doc <- xml2::read_xml(path, options = c("NONET", "NOBLANKS"))
  model_node <- xml2::xml_find_first(doc, "//*[local-name()='model']")
  if (inherits(model_node, "xml_missing")) stop("SBML model element is missing.", call. = FALSE)
  reaction_nodes <- xml2::xml_find_all(doc, "//*[local-name()='listOfReactions']/*[local-name()='reaction']")
  species_nodes <- xml2::xml_find_all(doc, "//*[local-name()='listOfSpecies']/*[local-name()='species']")
  gene_nodes <- xml2::xml_find_all(doc, "//*[local-name()='listOfGeneProducts']/*[local-name()='geneProduct']")
  reaction_ids <- .unique_ids(vapply(reaction_nodes, .xml_attr_local, character(1), local_name = "id"), "Reaction")
  metabolite_ids <- .unique_ids(vapply(species_nodes, .xml_attr_local, character(1), local_name = "id"), "Metabolite")
  gene_ids <- .unique_ids(vapply(gene_nodes, .xml_attr_local, character(1), local_name = "id"), "Gene")
  reaction_names <- stats::setNames(vapply(reaction_nodes, .xml_attr_local, character(1), local_name = "name", required = FALSE), reaction_ids)
  metabolite_names <- stats::setNames(vapply(species_nodes, .xml_attr_local, character(1), local_name = "name", required = FALSE), metabolite_ids)
  missing_metabolite_names <- is.na(metabolite_names) | !nzchar(metabolite_names)
  metabolite_names[missing_metabolite_names] <- metabolite_ids[missing_metabolite_names]
  compartments <- stats::setNames(vapply(species_nodes, .xml_attr_local, character(1), local_name = "compartment"), metabolite_ids)
  parameter_nodes <- xml2::xml_find_all(doc, "//*[local-name()='listOfParameters']/*[local-name()='parameter']")
  parameter_ids <- vapply(parameter_nodes, .xml_attr_local, character(1), local_name = "id")
  parameter_values <- as.numeric(vapply(parameter_nodes, .xml_attr_local, character(1), local_name = "value"))
  if (any(!is.finite(parameter_values))) stop("Non-finite SBML parameter.", call. = FALSE)
  parameters <- stats::setNames(parameter_values, parameter_ids)

  lb <- ub <- stats::setNames(numeric(length(reaction_ids)), reaction_ids)
  gpr <- stats::setNames(rep("", length(reaction_ids)), reaction_ids)
  i <- j <- numeric(0); x <- numeric(0)
  for (k in seq_along(reaction_nodes)) {
    r <- reaction_nodes[[k]]; rid <- reaction_ids[[k]]
    lb[[rid]] <- .resolve_bound(.xml_attr_local(r, "lowerFluxBound"), parameters, rid, "lower bound")
    ub[[rid]] <- .resolve_bound(.xml_attr_local(r, "upperFluxBound"), parameters, rid, "upper bound")
    if (lb[[rid]] > ub[[rid]]) stop("Contradictory source bounds in ", rid, call. = FALSE)
    direct_children <- xml2::xml_children(r)
    child_names <- sub("^.*:", "", xml2::xml_name(direct_children))
    association_index <- which(child_names == "geneProductAssociation")
    if (length(association_index) > 1L) stop("Multiple FBC gene associations in ", rid, call. = FALSE)
    if (length(association_index) == 1L) {
      association_children <- xml2::xml_children(direct_children[[association_index]])
      if (length(association_children) != 1L) stop("Malformed FBC gene association in ", rid, call. = FALSE)
      gpr[[rid]] <- .parse_gpr_node(association_children[[1L]], gene_ids)
    }
    for (list_name in c("listOfReactants", "listOfProducts")) {
      list_index <- which(child_names == list_name)
      if (length(list_index) > 1L) stop("Duplicate ", list_name, " in ", rid, call. = FALSE)
      if (!length(list_index)) next
      side <- if (list_name == "listOfReactants") -1 else 1
      refs <- xml2::xml_children(direct_children[[list_index]])
      for (ref in refs) {
        mid <- .xml_attr_local(ref, "species")
        if (!mid %in% metabolite_ids) stop("Unresolved metabolite '", mid, "' in ", rid, call. = FALSE)
        stoichiometry_text <- .xml_attr_local(ref, "stoichiometry", required = FALSE)
        if (is.na(stoichiometry_text) || !nzchar(stoichiometry_text)) stoichiometry_text <- "1"
        coefficient <- suppressWarnings(as.numeric(stoichiometry_text))
        if (!is.finite(coefficient)) stop("Unsupported stoichiometry in ", rid, call. = FALSE)
        i <- c(i, match(mid, metabolite_ids)); j <- c(j, k); x <- c(x, side * coefficient)
      }
    }
  }
  S <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(length(metabolite_ids), length(reaction_ids)), dimnames = list(metabolite_ids, reaction_ids))
  objectives <- xml2::xml_find_all(doc, "//*[local-name()='listOfObjectives']/*[local-name()='objective']")
  active <- .xml_attr_local(xml2::xml_find_first(doc, "//*[local-name()='listOfObjectives']"), "activeObjective")
  objective_nodes <- objectives[vapply(objectives, .xml_attr_local, character(1), local_name = "id") == active]
  if (length(objective_nodes) != 1L) stop("Missing/ambiguous active FBC objective.", call. = FALSE)
  objective_node <- objective_nodes[[1L]]
  objective_id <- .xml_attr_local(objective_node, "id")
  objective_sense <- .xml_attr_local(objective_node, "type")
  if (!objective_sense %in% c("maximize", "minimize")) stop("Unsupported objective sense.", call. = FALSE)
  objective <- stats::setNames(numeric(length(reaction_ids)), reaction_ids)
  flux_objectives <- xml2::xml_find_all(objective_node, ".//*[local-name()='fluxObjective']")
  for (fo in flux_objectives) {
    rid <- .xml_attr_local(fo, "reaction")
    if (!rid %in% reaction_ids) stop("Objective refers to unknown reaction: ", rid, call. = FALSE)
    objective[[rid]] <- objective[[rid]] + as.numeric(.xml_attr_local(fo, "coefficient"))
  }
  if (!length(flux_objectives)) stop("Active objective has no flux terms.", call. = FALSE)
  sha <- digest::digest(file = path, algo = "sha256")
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  file_basename <- basename(path)
  expected <- manifest[manifest$model_file == file_basename, , drop = FALSE]
  if (nrow(expected) != 1L) stop("Model missing or duplicated in manifest: ", file_basename, call. = FALSE)
  actual <- c(sha256 = sha, genes = length(gene_ids), reactions = length(reaction_ids), metabolites = length(metabolite_ids), active_objective = objective_id)
  for (field in names(actual)) {
    if (!identical(as.character(expected[[field]][[1L]]), as.character(actual[[field]]))) {
      stop("Manifest mismatch for ", file_basename, " field ", field, ": expected ", expected[[field]][[1L]], ", got ", actual[[field]], call. = FALSE)
    }
  }
  flux_unit <- .xml_attr_local(xml2::xml_find_first(doc, "//*[local-name()='unitDefinition' and @id='mmol_per_gDW_per_hr']"), "id", required = FALSE)
  model <- list(S = S, lb = lb, ub = ub, objective = objective, objective_id = objective_id, objective_sense = objective_sense,
    reaction_names = reaction_names, metabolite_names = metabolite_names,
    gpr = gpr, reaction_ids = reaction_ids, metabolite_ids = metabolite_ids, gene_ids = gene_ids,
    exchange_ids = reaction_ids[grepl("^(R_)?EX_", reaction_ids)], model_path = normalizePath(path),
    source_sha256 = sha, flux_unit = if (is.na(flux_unit)) NA_character_ else "mmol_per_gDW_per_hr",
    medium_id = NA_character_, medium_source = NA_character_, medium_config = NULL,
    effective_exchange_bounds = NULL, scenario_id = NA_character_, compartments = compartments)
  class(model) <- "session5_model"
  model
}

# Apply an explicitly documented medium to a copy. In these Recon1/COBRA-style
# exchange reactions, negative flux is uptake and positive flux is secretion.
# Unlisted uptake on single-metabolite extracellular exchanges is closed.
set_medium <- function(model, medium) {
  stopifnot(inherits(model, "session5_model"))
  if (!is.list(medium) || !all(c("id", "source", "bounds") %in% names(medium)) ||
      length(medium$id) != 1L || is.na(medium$id) || !nzchar(medium$id) ||
      length(medium$source) != 1L || is.na(medium$source) || !nzchar(medium$source)) {
    stop("Medium requires nonempty id, source, and bounds.", call. = FALSE)
  }
  b <- medium$bounds
  required <- c("reaction_id", "lower_bound", "upper_bound", "rationale")
  if (!is.data.frame(b) || !all(required %in% names(b)) || anyDuplicated(b$reaction_id) || anyNA(b[, required])) {
    stop("Medium bounds need unique reaction IDs, finite explicit bounds, and rationale.", call. = FALSE)
  }
  if (!all(is.finite(b$lower_bound)) || !all(is.finite(b$upper_bound)) || any(b$lower_bound > b$upper_bound)) {
    stop("Invalid or contradictory medium bounds.", call. = FALSE)
  }
  if (any(!b$reaction_id %in% model$reaction_ids)) stop("Medium contains unknown reaction IDs.", call. = FALSE)
  if (any(!b$reaction_id %in% model$exchange_ids)) stop("Medium may only set exchange reaction bounds.", call. = FALSE)
  out <- model
  for (rid in model$exchange_ids) {
    col <- model$S[, rid]
    nz <- which(col != 0)
    external <- nz[model$compartments[rownames(model$S)[nz]] == "e"]
    if (length(nz) != 1L || length(external) != 1L) next
    if (col[external] < 0) out$lb[[rid]] <- max(out$lb[[rid]], 0) else out$ub[[rid]] <- min(out$ub[[rid]], 0)
  }
  for (row in seq_len(nrow(b))) {
    rid <- b$reaction_id[[row]]
    nz <- which(model$S[, rid] != 0)
    if (length(nz) != 1L || model$compartments[rownames(model$S)[nz]] != "e") {
      stop("Configured reaction is not a single-metabolite extracellular exchange: ", rid, call. = FALSE)
    }
    out$lb[[rid]] <- b$lower_bound[[row]]
    out$ub[[rid]] <- b$upper_bound[[row]]
  }
  out$medium_id <- medium$id
  out$medium_source <- medium$source
  out$medium_config <- medium
  out$effective_exchange_bounds <- data.frame(reaction_id = model$exchange_ids,
    lower_bound = unname(out$lb[model$exchange_ids]), upper_bound = unname(out$ub[model$exchange_ids]),
    stringsAsFactors = FALSE)
  out
}

# Shared HiGHS adapter for FBA and FVA; optional linear constraints are in
# original flux variables. lpSolve returned false "unbounded" on bounded Recon1 LPs.
.solve_lp <- function(model, objective, sense, tolerance, extra_constraint = NULL) {
  base <- list(status = "error", objective_value = NA_real_,
    fluxes = stats::setNames(rep(NA_real_, length(model$reaction_ids)), model$reaction_ids),
    max_abs_sv = NA_real_, solver = "highs", tolerance = tolerance,
    model_sha256 = model$source_sha256, medium_id = model$medium_id)
  if (!requireNamespace("highs", quietly = TRUE)) { base$error <- "highs package unavailable"; return(base) }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) { base$error <- "tolerance must be positive and finite"; return(base) }
  if (length(objective) != length(model$reaction_ids) || any(!is.finite(objective)) || !sense %in% c("maximize", "minimize")) {
    base$error <- "LP objective or sense is invalid"; return(base)
  }
  if (any(!is.finite(model$lb)) || any(!is.finite(model$ub)) || any(model$lb > model$ub)) {
    base$error <- "LP adapter requires finite, non-contradictory reaction bounds"; return(base)
  }
  n <- length(model$reaction_ids)
  A <- methods::as(model$S, "CsparseMatrix")
  lhs <- rhs <- rep(0, nrow(A))
  if (!is.null(extra_constraint)) {
    if (length(extra_constraint$coefficients) != n || any(!is.finite(extra_constraint$coefficients)) ||
        length(extra_constraint$rhs) != 1L || !is.finite(extra_constraint$rhs) ||
        !extra_constraint$direction %in% c("<=", ">=", "=")) {
      base$error <- "Additional LP constraint is invalid"; return(base)
    }
    A <- rbind(A, Matrix::Matrix(extra_constraint$coefficients, nrow = 1L, sparse = TRUE))
    lhs <- c(lhs, if (extra_constraint$direction == "<=") -Inf else extra_constraint$rhs)
    rhs <- c(rhs, if (extra_constraint$direction == ">=") Inf else extra_constraint$rhs)
  }
  fit <- tryCatch(highs::highs_solve(L = unname(objective), lower = unname(model$lb), upper = unname(model$ub),
    A = A, lhs = lhs, rhs = rhs, maximum = identical(sense, "maximize"),
    control = highs::highs_control(threads = 1L, log_to_console = FALSE,
      primal_feasibility_tolerance = tolerance, dual_feasibility_tolerance = tolerance)), error = identity)
  if (inherits(fit, "error")) { base$error <- conditionMessage(fit); return(base) }
  base$solver_code <- as.integer(fit$status)
  message <- tolower(as.character(fit$status_message %||% ""))
  base$status <- if (identical(message, "optimal")) "optimal" else if (grepl("infeasible", message)) "infeasible" else
    if (grepl("unbounded", message)) "unbounded" else if (grepl("time limit", message)) "time_limit" else "error"
  if (base$status == "optimal") {
    fluxes <- stats::setNames(fit$primal_solution, model$reaction_ids)
    base$fluxes <- fluxes
    base$objective_value <- sum(objective * fluxes)
    base$max_abs_sv <- max(abs(as.numeric(model$S %*% fluxes)))
  } else {
    base$error <- paste("HiGHS status:", fit$status_message)
  }
  base
}

solve_fba <- function(model, tolerance = 1e-7) {
  if (!inherits(model, "session5_model")) stop("model must inherit from session5_model.", call. = FALSE)
  result <- .solve_lp(model, model$objective, model$objective_sense, tolerance)
  result$objective_id <- model$objective_id
  result
}

# Limited-panel FVA: one min and max LP per verified reaction, conditional on
# retaining the requested fraction of the active objective optimum.
run_fva <- function(model, reaction_ids, fraction = 0.9, tolerance = 1e-7) {
  if (!inherits(model, "session5_model")) stop("model must inherit from session5_model.", call. = FALSE)
  if (!is.character(reaction_ids) || !length(reaction_ids) || anyNA(reaction_ids) ||
      any(!nzchar(reaction_ids)) || anyDuplicated(reaction_ids)) stop("reaction_ids must be unique, nonempty reaction IDs.", call. = FALSE)
  unknown <- setdiff(reaction_ids, model$reaction_ids)
  if (length(unknown)) stop("Unknown reaction ID(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  if (length(fraction) != 1L || !is.numeric(fraction) || !is.finite(fraction) || fraction < 0 || fraction > 1) {
    stop("fraction must be one finite value in [0, 1].", call. = FALSE)
  }
  if (length(tolerance) != 1L || !is.numeric(tolerance) || !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be one positive finite value.", call. = FALSE)
  }
  baseline <- solve_fba(model, tolerance)
  n <- length(reaction_ids)
  intervals <- data.frame(reaction_id = reaction_ids, minimum = rep(NA_real_, n), maximum = rep(NA_real_, n),
    fba_flux = rep(NA_real_, n), min_status = rep("not_run", n), max_status = rep("not_run", n),
    min_solver_code = rep(NA_integer_, n), max_solver_code = rep(NA_integer_, n), stringsAsFactors = FALSE)
  result <- list(status = baseline$status, baseline_fba_status = baseline$status,
    baseline_objective = baseline$objective_value, objective_id = model$objective_id,
    objective_sense = model$objective_sense, objective = model$objective,
    model_sha256 = model$source_sha256, model_path = model$model_path,
    medium_id = model$medium_id,
    medium_source = if (is.null(model$medium_source) || !length(model$medium_source)) NA_character_ else model$medium_source,
    scenario_id = if (is.null(model$scenario_id) || !length(model$scenario_id)) NA_character_ else model$scenario_id,
    imat_settings_id = if (is.null(model$imat_settings_id) || !length(model$imat_settings_id)) NA_character_ else model$imat_settings_id,
    fraction = fraction, tolerance = tolerance, reaction_ids = reaction_ids,
    fba_fluxes = baseline$fluxes, intervals = intervals, solver = baseline$solver,
    effective_exchange_bounds = if (is.null(model$effective_exchange_bounds)) NULL else model$effective_exchange_bounds)
  if (baseline$status != "optimal") {
    result$intervals$min_status <- baseline$status
    result$intervals$max_status <- baseline$status
    result$error <- paste("FVA requires optimal baseline FBA; got", baseline$status)
    return(result)
  }
  slack <- (1 - fraction) * abs(baseline$objective_value)
  retained <- list(coefficients = unname(model$objective),
    direction = if (model$objective_sense == "maximize") ">=" else "<=",
    rhs = if (model$objective_sense == "maximize") baseline$objective_value - slack else baseline$objective_value + slack)
  for (i in seq_along(reaction_ids)) {
    objective <- stats::setNames(numeric(length(model$reaction_ids)), model$reaction_ids)
    objective[[reaction_ids[[i]]]] <- 1
    lower <- .solve_lp(model, objective, "minimize", tolerance, retained)
    upper <- .solve_lp(model, objective, "maximize", tolerance, retained)
    intervals$min_status[[i]] <- lower$status
    intervals$max_status[[i]] <- upper$status
    intervals$min_solver_code[[i]] <- if (is.null(lower$solver_code)) NA_integer_ else lower$solver_code
    intervals$max_solver_code[[i]] <- if (is.null(upper$solver_code)) NA_integer_ else upper$solver_code
    intervals$fba_flux[[i]] <- baseline$fluxes[[reaction_ids[[i]]]]
    if (lower$status == "optimal") intervals$minimum[[i]] <- lower$objective_value
    if (upper$status == "optimal") intervals$maximum[[i]] <- upper$objective_value
  }
  result$intervals <- intervals
  valid <- intervals$min_status == "optimal" & intervals$max_status == "optimal"
  result$status <- if (all(valid)) "optimal" else if (any(valid)) "partial" else "error"
  if (result$status != "optimal") result$error <- "One or more FVA endpoint LPs did not return optimal status."
  result
}

# Parsimonious FBA: the minimal-total-flux vector retaining `fraction` of the optimum.
solve_pfba <- function(model, fraction = 1, tolerance = 1e-7) {
  if (!inherits(model, "session5_model")) stop("model must inherit from session5_model.", call. = FALSE)
  if (length(fraction) != 1L || !is.finite(fraction) || fraction < 0 || fraction > 1) stop("fraction must be one finite value in [0, 1].", call. = FALSE)
  baseline <- solve_fba(model, tolerance)
  result <- list(status = baseline$status, objective_id = model$objective_id, objective_value = NA_real_,
    baseline_objective = baseline$objective_value, fraction = fraction, total_flux = NA_real_,
    fluxes = baseline$fluxes * NA_real_, max_abs_sv = NA_real_, solver = "highs", tolerance = tolerance,
    model_sha256 = model$source_sha256, medium_id = model$medium_id)
  if (baseline$status != "optimal") { result$error <- paste("pFBA requires optimal FBA; got", baseline$status); return(result) }
  n <- length(model$reaction_ids)
  S <- methods::as(model$S, "CsparseMatrix")
  I <- Matrix::Diagonal(n)
  zero <- Matrix::Matrix(0, nrow(S), n, sparse = TRUE)
  objective_row <- Matrix::Matrix(unname(model$objective), nrow = 1L, sparse = TRUE)
  slack <- (1 - fraction) * abs(baseline$objective_value)
  retained <- if (model$objective_sense == "maximize") c(baseline$objective_value - slack, Inf) else c(-Inf, baseline$objective_value + slack)
  A <- rbind(cbind(S, zero), cbind(objective_row, Matrix::Matrix(0, 1L, n, sparse = TRUE)), cbind(I, I), cbind(-I, I))
  lhs <- c(rep(0, nrow(S)), retained[[1L]], rep(0, 2L * n))
  rhs <- c(rep(0, nrow(S)), retained[[2L]], rep(Inf, 2L * n))
  fit <- tryCatch(highs::highs_solve(L = c(rep(0, n), rep(1, n)), lower = c(unname(model$lb), rep(0, n)),
    upper = c(unname(model$ub), rep(Inf, n)), A = A, lhs = lhs, rhs = rhs, maximum = FALSE,
    control = highs::highs_control(threads = 1L, log_to_console = FALSE,
      primal_feasibility_tolerance = tolerance, dual_feasibility_tolerance = tolerance)), error = identity)
  if (inherits(fit, "error") || !identical(tolower(fit$status_message), "optimal")) {
    result$status <- "error"; result$error <- if (inherits(fit, "error")) conditionMessage(fit) else fit$status_message; return(result)
  }
  fluxes <- stats::setNames(fit$primal_solution[seq_len(n)], model$reaction_ids)
  result$fluxes <- fluxes
  result$objective_value <- sum(model$objective * fluxes)
  result$total_flux <- sum(abs(fluxes))
  result$max_abs_sv <- max(abs(as.numeric(model$S %*% fluxes)))
  result
}

# Apply an O2 uptake upper bound; this does not prescribe actual oxygen consumption or OXPHOS flux.
apply_oxygen_bound_scenario <- function(model, settings) {
  stopifnot(inherits(model, "session5_model"))
  required <- c("id", "o2_uptake_capacity", "fix_glucose_at_cap", "closed_secretions")
  if (!is.list(settings) || !all(required %in% names(settings))) stop("Respiratory-capacity settings are incomplete.", call. = FALSE)
  capacity <- settings$o2_uptake_capacity
  if (length(capacity) != 1L || !is.finite(capacity) || capacity < 0) stop("o2_uptake_capacity must be one non-negative number.", call. = FALSE)
  closed <- if (is.character(settings$closed_secretions)) trimws(strsplit(settings$closed_secretions, ",")[[1L]]) else character()
  closed <- closed[nzchar(closed)]
  if (!all(c("R_EX_o2_e", "R_EX_glc__D_e", closed) %in% model$exchange_ids)) stop("Capacity settings refer to unknown exchanges.", call. = FALSE)
  out <- model
  out$lb[["R_EX_o2_e"]] <- max(out$lb[["R_EX_o2_e"]], -capacity)
  if (isTRUE(settings$fix_glucose_at_cap)) out$ub[["R_EX_glc__D_e"]] <- out$lb[["R_EX_glc__D_e"]]
  out$ub[closed] <- 0
  if (any(out$lb > out$ub)) stop("Respiratory-capacity settings contradict existing bounds.", call. = FALSE)
  out$capacity_settings_id <- settings$id
  out
}

# Read the flat scalar settings used by the shared iMAT-like MILP. This avoids
# introducing a YAML package dependency for a small, versioned settings file.
read_imat_settings <- function(path = "config/imat_settings.yml") {
  if (!file.exists(path)) stop("iMAT settings file not found: ", path, call. = FALSE)
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(sub("#.*$", "", lines))
  lines <- lines[nzchar(lines)]
  split <- regmatches(lines, regexec("^([A-Za-z_][A-Za-z0-9_]*):\\s*(.*?)\\s*$", lines))
  if (any(lengths(split) != 3L)) stop("Unsupported iMAT settings YAML; expected flat key: scalar entries.", call. = FALSE)
  keys <- vapply(split, `[[`, character(1), 2L)
  values <- vapply(split, `[[`, character(1), 3L)
  if (anyDuplicated(keys)) stop("Duplicate iMAT settings key.", call. = FALSE)
  parsed <- lapply(values, function(value) {
    value <- gsub("^['\"]|['\"]$", "", value)
    if (tolower(value) %in% c("true", "false")) return(tolower(value) == "true")
    numeric_value <- suppressWarnings(as.numeric(value))
    if (!is.na(numeric_value)) return(numeric_value)
    value
  })
  stats::setNames(parsed, keys)
}

# Load a single named medium from the versioned shared CSV. Extra metadata
# columns are retained in the file; set_medium() receives its stable contract.
read_medium <- function(path = "config/medium_bounds.csv", medium_id) {
  if (length(medium_id) != 1L || is.na(medium_id) || !nzchar(medium_id)) stop("medium_id must be one nonempty string.", call. = FALSE)
  if (!file.exists(path)) stop("Medium configuration not found: ", path, call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("medium_id", "reaction_id", "lower_bound", "upper_bound", "rationale", "source")
  if (!all(required %in% names(x))) stop("Medium configuration has an invalid schema.", call. = FALSE)
  x <- x[x$medium_id == medium_id, , drop = FALSE]
  if (!nrow(x)) stop("No rows for medium_id: ", medium_id, call. = FALSE)
  sources <- unique(x$source)
  if (length(sources) != 1L || is.na(sources) || !nzchar(sources)) stop("A medium must have exactly one nonempty source.", call. = FALSE)
  list(id = medium_id, source = sources[[1L]], bounds = x[, c("reaction_id", "lower_bound", "upper_bound", "rationale"), drop = FALSE])
}

# GPR evidence policy shared by all Session 5 exercises. Unknown values remain
# NA except where low AND or high OR evidence is decisive.
reaction_evidence_from_gpr <- function(gpr, gene_evidence) {
  if (length(gpr) != 1L || is.na(gpr)) stop("gpr must be one non-missing string.", call. = FALSE)
  if (!is.numeric(gene_evidence) || is.null(names(gene_evidence)) ||
      anyNA(names(gene_evidence)) || any(!nzchar(names(gene_evidence))) || anyDuplicated(names(gene_evidence))) {
    stop("gene_evidence must be a uniquely named numeric vector.", call. = FALSE)
  }
  if (any(!is.na(gene_evidence) & !gene_evidence %in% c(-1, 0, 1))) stop("gene_evidence values must be -1, 0, 1, or NA.", call. = FALSE)
  if (!nzchar(trimws(gpr))) return(NA_real_)

  separated <- gsub("([()])", " \\1 ", trimws(gpr), perl = TRUE)
  tokens <- strsplit(trimws(separated), "\\s+", perl = TRUE)[[1L]]
  position <- 1L
  peek <- function() if (position <= length(tokens)) tokens[[position]] else ""
  combine <- function(operator, values) {
    if (operator == "and") {
      if (any(!is.na(values) & values == -1)) return(-1)
      if (anyNA(values)) return(NA_real_)
      return(min(values))
    }
    if (any(!is.na(values) & values == 1)) return(1)
    if (anyNA(values)) return(NA_real_)
    max(values)
  }
  parse_or <- NULL
  parse_primary <- function() {
    token <- peek()
    if (!nzchar(token)) stop("Unexpected end of GPR.", call. = FALSE)
    if (token == "(") {
      position <<- position + 1L
      value <- parse_or()
      if (peek() != ")") stop("Unclosed parenthesis in GPR.", call. = FALSE)
      position <<- position + 1L
      return(value)
    }
    if (token %in% c(")", "and", "or")) stop("Expected a gene product in GPR.", call. = FALSE)
    position <<- position + 1L
    if (token %in% names(gene_evidence)) unname(gene_evidence[[token]]) else NA_real_
  }
  parse_and <- function() {
    values <- parse_primary()
    while (peek() == "and") {
      position <<- position + 1L
      values <- c(values, parse_primary())
    }
    if (length(values) == 1L) values[[1L]] else combine("and", values)
  }
  parse_or <- function() {
    values <- parse_and()
    while (peek() == "or") {
      position <<- position + 1L
      values <- c(values, parse_and())
    }
    if (length(values) == 1L) values[[1L]] else combine("or", values)
  }
  result <- parse_or()
  if (position <= length(tokens)) stop("Unexpected token in GPR: ", peek(), call. = FALSE)
  as.numeric(result)
}

# Preindex GPR gene tokens once for larger screens; immutable GPR equality is
# checked by knock_out_genes() before any cached index is reused.
ko_gpr_index <- function(model) {
  if (!inherits(model, "session5_model") || !identical(names(model$gpr), model$reaction_ids))
    stop("Model GPRs must align with reaction IDs.", call. = FALSE)
  rules <- model$gpr[nzchar(trimws(model$gpr))]
  by_gene <- stats::setNames(vector("list", length(model$gene_ids)), model$gene_ids)
  for (rid in names(rules)) {
    tokens <- setdiff(unique(strsplit(trimws(gsub("[()]", " ", rules[[rid]])), "\\s+")[[1L]]), c("and", "or"))
    if (any(!tokens %in% model$gene_ids)) stop("GPR contains unknown model gene-product IDs.", call. = FALSE)
    for (g in tokens) by_gene[[g]] <- c(by_gene[[g]], rid)
  }
  list(gpr = model$gpr, gene_ids = model$gene_ids, by_gene = by_gene)
}

# Boolean deletion on a fresh WT copy. Explicitly unavailable genes are a
# caller-supplied context assumption, NOT inferred from continuous expression.
# A KO closes only reactions whose unchanged GPR changes from true to false.
knock_out_genes <- function(model, gene_ids, index = NULL, unavailable_gene_ids = character()) {
  if (!inherits(model, "session5_model")) stop("model must inherit from session5_model.", call. = FALSE)
  if (!is.character(gene_ids) || !length(gene_ids) || anyNA(gene_ids) ||
      any(!nzchar(gene_ids)) || anyDuplicated(gene_ids)) {
    stop("gene_ids must be unique, nonempty gene-product IDs.", call. = FALSE)
  }
  if (any(!gene_ids %in% model$gene_ids)) stop("Unknown model gene-product ID(s): ",
    paste(setdiff(gene_ids, model$gene_ids), collapse = ", "), call. = FALSE)
  if (!is.character(unavailable_gene_ids) || anyNA(unavailable_gene_ids) ||
      anyDuplicated(unavailable_gene_ids) || any(!unavailable_gene_ids %in% model$gene_ids))
    stop("unavailable_gene_ids must contain unique, known model gene IDs.", call. = FALSE)
  if (!identical(names(model$gpr), model$reaction_ids) ||
      !identical(names(model$lb), model$reaction_ids) || !identical(names(model$ub), model$reaction_ids)) {
    stop("Model GPR/bounds must align with reaction IDs.", call. = FALSE)
  }
  if (is.null(index)) index <- ko_gpr_index(model)
  if (!identical(index$gpr, model$gpr) || !identical(index$gene_ids, model$gene_ids))
    stop("Stale GPR index does not match the input model.", call. = FALSE)
  uncovered <- gene_ids[lengths(index$by_gene[gene_ids]) == 0L]
  if (length(uncovered)) stop("Gene ID(s) without GPR coverage: ", paste(uncovered, collapse = ", "), call. = FALSE)
  evidence <- stats::setNames(rep(1, length(model$gene_ids)), model$gene_ids)
  evidence[unavailable_gene_ids] <- -1
  if (length(unavailable_gene_ids)) {
    preblocked <- vapply(model$gpr[nzchar(model$gpr)], function(rule)
      identical(reaction_evidence_from_gpr(rule, evidence), -1), logical(1))
    already_false <- names(preblocked)[preblocked]
    if (any(model$lb[already_false] != 0 | model$ub[already_false] != 0))
      stop("Context-unavailable genes require their already-false GPR reactions to be closed in WT.", call. = FALSE)
  }
  possible <- unique(unlist(index$by_gene[gene_ids], use.names = FALSE))
  active_before <- vapply(model$gpr[possible], function(rule)
    identical(reaction_evidence_from_gpr(rule, evidence), 1), logical(1))
  evidence[gene_ids] <- -1
  inactive <- active_before & vapply(model$gpr[possible], function(rule)
    identical(reaction_evidence_from_gpr(rule, evidence), -1), logical(1))
  affected <- model$reaction_ids[model$reaction_ids %in% possible[inactive]]
  if (any(model$lb[affected] > 0 | model$ub[affected] < 0)) {
    stop("Cannot close an affected reaction with a compulsory nonzero flux without reopening its bounds.", call. = FALSE)
  }
  out <- model
  out$lb[affected] <- pmax(model$lb[affected], 0)
  out$ub[affected] <- pmin(model$ub[affected], 0)
  list(model = out, deleted_genes = gene_ids, affected_reactions = affected)
}

# A tested, reversible-aware iMAT-like MILP. High evidence rewards a flux whose
# signed magnitude exceeds epsilon; low evidence rewards inactivity within
# epsilon. Intermediate/unknown evidence is neutral. Biomass is constrained
# to a configured fraction of the same-medium FBA optimum.
reconstruct_context <- function(model, gene_evidence, medium, config) {
  stopifnot(inherits(model, "session5_model"))
  required_config <- c("id", "epsilon", "growth_fraction", "flux_bound_cap", "time_limit_seconds", "unknown_evidence")
  if (!is.list(config) || !all(required_config %in% names(config))) stop("iMAT config is incomplete.", call. = FALSE)
  if (!identical(config$unknown_evidence, "neutral")) stop("Only neutral unknown-evidence handling is supported.", call. = FALSE)
  if (is.null(names(gene_evidence)) || anyDuplicated(names(gene_evidence)) || any(!names(gene_evidence) %in% model$gene_ids)) stop("Gene evidence must be uniquely named with model gene-product IDs.", call. = FALSE)
  if (any(!is.na(gene_evidence) & !gene_evidence %in% c(-1, 0, 1))) stop("Gene evidence values must be -1, 0, 1, or NA.", call. = FALSE)
  if (!is.finite(config$epsilon) || config$epsilon <= 0 || !is.finite(config$growth_fraction) || config$growth_fraction <= 0 || config$growth_fraction > 1 ||
      !is.finite(config$flux_bound_cap) || config$flux_bound_cap <= config$epsilon || !is.finite(config$time_limit_seconds) || config$time_limit_seconds <= 0) stop("Invalid iMAT numeric settings.", call. = FALSE)
  active_objective_reactions <- names(model$objective)[model$objective != 0]
  objective_is_biomass <- grepl("biomass", model$objective_id, ignore.case = TRUE) || "R_biomass_reaction" %in% active_objective_reactions
  if (!identical(model$objective_sense, "maximize") || !objective_is_biomass) stop("Context reconstruction requires the active biomass-maximization objective.", call. = FALSE)
  solver <- if (is.null(config$solver)) "highs" else config$solver
  if (!solver %in% c("highs", "lpSolve")) stop("Unsupported configured MILP solver: ", solver, call. = FALSE)
  if (solver == "highs" && !requireNamespace("highs", quietly = TRUE)) stop("Configured HiGHS backend requires the CRAN 'highs' package.", call. = FALSE)
  if (solver == "lpSolve" && !requireNamespace("lpSolve", quietly = TRUE)) stop("Configured lpSolve backend requires the 'lpSolve' package.", call. = FALSE)
  runtime_solver_version <- if (solver == "highs") as.character(utils::packageVersion("highs")) else as.character(utils::packageVersion("lpSolve"))
  if (!is.null(config$solver_version) && !identical(as.character(config$solver_version), runtime_solver_version)) stop(solver, " package version mismatch: settings require ", config$solver_version, ", runtime has ", runtime_solver_version, call. = FALSE)
  runtime_engine_version <- if (solver == "highs") {
    probe_model <- highs::highs_model(L = 0, lower = 0)
    highs::hi_solver_version(highs::hi_new_solver(probe_model))
  } else runtime_solver_version
  if (!is.null(config$solver_engine_version) && !identical(as.character(config$solver_engine_version), as.character(runtime_engine_version))) stop(solver, " engine version mismatch: settings require ", config$solver_engine_version, ", runtime has ", runtime_engine_version, call. = FALSE)
  if (!is.null(config$objective_id) && !identical(config$objective_id, model$objective_id)) stop("Configured objective ID does not match the loaded model.", call. = FALSE)
  if (!is.null(config$objective_sense) && !identical(config$objective_sense, model$objective_sense)) stop("Configured objective sense does not match the loaded model.", call. = FALSE)
  if (!is.numeric(gene_evidence)) stop("gene_evidence must be numeric.", call. = FALSE)
  high_weight <- if (is.null(config$high_weight)) 1 else config$high_weight
  low_weight <- if (is.null(config$low_weight)) 1 else config$low_weight
  if (!is.numeric(high_weight) || length(high_weight) != 1L || !is.finite(high_weight) || high_weight <= 0 ||
      !is.numeric(low_weight) || length(low_weight) != 1L || !is.finite(low_weight) || low_weight <= 0) stop("Evidence weights must be positive finite scalars.", call. = FALSE)

  bounded_model <- set_medium(model, medium)
  bounded_model$lb <- pmax(bounded_model$lb, -config$flux_bound_cap)
  bounded_model$ub <- pmin(bounded_model$ub, config$flux_bound_cap)
  if (any(bounded_model$lb > bounded_model$ub)) {
    stop("Medium bounds conflict with the configured flux-bound cap.", call. = FALSE)
  }
  bounded_model$effective_exchange_bounds <- data.frame(
    reaction_id = bounded_model$exchange_ids,
    lower_bound = unname(bounded_model$lb[bounded_model$exchange_ids]),
    upper_bound = unname(bounded_model$ub[bounded_model$exchange_ids]),
    stringsAsFactors = FALSE
  )
  baseline <- solve_fba(bounded_model)
  if (baseline$status != "optimal" || !is.finite(baseline$objective_value) || baseline$objective_value <= 0) stop("Same-medium biomass preflight failed (status=", baseline$status, ", objective=", baseline$objective_value, ").", call. = FALSE)

  reaction_evidence <- stats::setNames(vapply(model$reaction_ids, function(rid) reaction_evidence_from_gpr(model$gpr[[rid]], gene_evidence), numeric(1)), model$reaction_ids)
  high_ids <- which(reaction_evidence == 1 & !is.na(reaction_evidence))
  low_ids <- which(reaction_evidence == -1 & !is.na(reaction_evidence))
  # Omit impossible high-evidence directions at the variable-generation stage.
  # One-sided reactions need one binary; only reversible reactions need both.
  high_always_active_ids <- high_ids[bounded_model$lb[high_ids] >= config$epsilon | bounded_model$ub[high_ids] <= -config$epsilon]
  high_positive_ids <- high_ids[bounded_model$ub[high_ids] >= config$epsilon & bounded_model$lb[high_ids] < config$epsilon]
  high_negative_ids <- high_ids[bounded_model$lb[high_ids] <= -config$epsilon & bounded_model$ub[high_ids] > -config$epsilon]
  high_both_ids <- intersect(high_positive_ids, high_negative_ids)
  low_all_ids <- low_ids
  low_always_inactive_ids <- low_all_ids[bounded_model$lb[low_all_ids] >= -config$epsilon & bounded_model$ub[low_all_ids] <= config$epsilon]
  low_ids <- setdiff(low_all_ids, low_always_inactive_ids)
  n_flux <- length(model$reaction_ids); n_positive <- length(high_positive_ids); n_negative <- length(high_negative_ids); n_both <- length(high_both_ids); n_low <- length(low_ids)
  pos_start <- n_flux + 1L
  neg_start <- pos_start + n_positive
  low_start <- neg_start + n_negative
  n_var <- n_flux + n_positive + n_negative + n_low
  row_count <- nrow(model$S) + n_flux + 1L + n_positive + n_negative + n_both + 2L * n_low
  constraints <- Matrix::Matrix(0, nrow = row_count, ncol = n_var, sparse = TRUE)
  rhs <- numeric(row_count)
  direction <- rep("<=", row_count)
  row <- 0L

  row_idx <- seq_len(nrow(model$S))
  constraints[row_idx, seq_len(n_flux)] <- model$S
  rhs[row_idx] <- -as.numeric(model$S %*% bounded_model$lb)
  direction[row_idx] <- "="
  row <- nrow(model$S)

  bound_rows <- row + seq_len(n_flux)
  constraints[cbind(bound_rows, seq_len(n_flux))] <- 1
  rhs[bound_rows] <- bounded_model$ub - bounded_model$lb
  row <- max(bound_rows)

  growth_row <- row + 1L
  constraints[growth_row, seq_len(n_flux)] <- unname(model$objective)
  rhs[growth_row] <- config$growth_fraction * baseline$objective_value - sum(model$objective * bounded_model$lb)
  direction[growth_row] <- ">="
  row <- growth_row

  for (j in seq_along(high_positive_ids)) {
    ridx <- high_positive_ids[[j]]; p <- pos_start + j - 1L
    big_m <- max(abs(bounded_model$lb[[ridx]]), abs(bounded_model$ub[[ridx]])) + config$epsilon
    row <- row + 1L
    constraints[row, ridx] <- 1; constraints[row, p] <- -big_m
    rhs[[row]] <- config$epsilon - big_m - bounded_model$lb[[ridx]]; direction[[row]] <- ">="
  }
  for (j in seq_along(high_negative_ids)) {
    ridx <- high_negative_ids[[j]]; n <- neg_start + j - 1L
    big_m <- max(abs(bounded_model$lb[[ridx]]), abs(bounded_model$ub[[ridx]])) + config$epsilon
    row <- row + 1L
    constraints[row, ridx] <- 1; constraints[row, n] <- big_m
    rhs[[row]] <- -config$epsilon + big_m - bounded_model$lb[[ridx]]
  }
  for (ridx in high_both_ids) {
    p <- pos_start + match(ridx, high_positive_ids) - 1L
    n <- neg_start + match(ridx, high_negative_ids) - 1L
    row <- row + 1L
    constraints[row, p] <- 1; constraints[row, n] <- 1
    rhs[[row]] <- 1
  }
  for (j in seq_along(low_ids)) {
    ridx <- low_ids[[j]]; z <- low_start + j - 1L
    big_m <- max(abs(bounded_model$lb[[ridx]]), abs(bounded_model$ub[[ridx]])) + config$epsilon
    row <- row + 1L
    constraints[row, ridx] <- 1; constraints[row, z] <- big_m
    rhs[[row]] <- config$epsilon + big_m - bounded_model$lb[[ridx]]
    row <- row + 1L
    constraints[row, ridx] <- 1; constraints[row, z] <- -big_m
    rhs[[row]] <- -config$epsilon - big_m - bounded_model$lb[[ridx]]; direction[[row]] <- ">="
  }
  if (row != row_count) stop("Internal iMAT constraint-count mismatch.", call. = FALSE)

  objective <- numeric(n_var)
  high_positive <- if (n_positive) seq.int(pos_start, pos_start + n_positive - 1L) else integer()
  high_negative <- if (n_negative) seq.int(neg_start, neg_start + n_negative - 1L) else integer()
  low_inactive <- if (n_low) seq.int(low_start, n_var) else integer()
  if (length(high_positive)) objective[high_positive] <- high_weight
  if (length(high_negative)) objective[high_negative] <- high_weight
  if (length(low_inactive)) objective[low_inactive] <- low_weight
  binary <- c(high_positive, high_negative, low_inactive)
  fit <- tryCatch({
    if (solver == "highs") {
      variable_upper <- rep(Inf, n_var)
      variable_upper[binary] <- 1
      variable_types <- rep("C", n_var)
      variable_types[binary] <- "I"
      lhs <- ifelse(direction == ">=", rhs, ifelse(direction == "=", rhs, -Inf))
      solver_rhs <- ifelse(direction == "<=", rhs, ifelse(direction == "=", rhs, Inf))
      highs::highs_solve(
        L = objective, lower = rep(0, n_var), upper = variable_upper,
        A = constraints, lhs = lhs, rhs = solver_rhs, types = variable_types,
        maximum = TRUE,
        control = highs::highs_control(
          threads = if (is.null(config$threads)) 1L else as.integer(config$threads),
          time_limit = config$time_limit_seconds,
          log_to_console = FALSE,
          primal_feasibility_tolerance = if (is.null(config$feasibility_tolerance)) 1e-7 else config$feasibility_tolerance,
          dual_feasibility_tolerance = if (is.null(config$feasibility_tolerance)) 1e-7 else config$feasibility_tolerance,
          mip_feasibility_tolerance = if (is.null(config$mip_feasibility_tolerance)) 1e-7 else config$mip_feasibility_tolerance
        )
      )
    } else {
      triplets <- as.matrix(Matrix::summary(constraints))
      solver_triplets <- cbind(triplets[, "i"], triplets[, "j"], triplets[, "x"])
      lp_call <- function() lpSolve::lp(direction = "max", objective.in = objective,
        const.dir = direction, const.rhs = rhs, dense.const = solver_triplets,
        binary.vec = binary, timeout = as.integer(config$time_limit_seconds),
        presolve = 1, scale = 1)
      if (.Platform$OS.type == "windows") lp_call() else {
        worker <- parallel::mcparallel(tryCatch(list(ok = TRUE, result = lp_call()),
          error = function(e) list(ok = FALSE, error = conditionMessage(e))), silent = TRUE, mc.set.seed = FALSE)
        completed <- parallel::mccollect(worker, wait = TRUE, timeout = config$time_limit_seconds)
        if (is.null(completed)) {
          try(tools::pskill(worker$pid), silent = TRUE)
          list(status = 7L, solution = numeric(), timeout_enforced = TRUE)
        } else {
          worker_result <- completed[[1L]]
          if (!isTRUE(worker_result$ok)) stop(worker_result$error, call. = FALSE)
          worker_result$result
        }
      }
    }
  }, error = identity)
  if (inherits(fit, "error")) stop("iMAT MILP solver error: ", conditionMessage(fit), call. = FALSE)
  if (solver == "highs") {
    status_message <- if (is.null(fit$status_message) || !length(fit$status_message) || is.na(fit$status_message[[1L]])) "" else tolower(as.character(fit$status_message[[1L]]))
    status <- if (grepl("optimal", status_message)) "optimal" else if (grepl("time.limit|time limit|timelimit", status_message)) "time_limit" else if (grepl("infeasible", status_message)) "infeasible" else if (grepl("unbounded", status_message)) "unbounded" else "error"
    solver_code <- fit$status
    solver_solution <- fit$primal_solution
    timeout_enforced <- status == "time_limit"
  } else {
    status <- switch(as.character(fit$status), `0` = "optimal", `2` = "infeasible", `3` = "unbounded", `7` = "time_limit", "error")
    status_message <- paste("lpSolve status", fit$status)
    solver_code <- fit$status
    solver_solution <- fit$solution
    timeout_enforced <- identical(as.integer(fit$status), 7L)
  }
  qc <- data.frame(
    reactions_total = n_flux,
    with_gpr = sum(nzchar(model$gpr)),
    evaluated = sum(!is.na(reaction_evidence)),
    unknown = sum(is.na(reaction_evidence)),
    positive = sum(reaction_evidence == 1, na.rm = TRUE),
    intermediate = sum(reaction_evidence == 0, na.rm = TRUE),
    negative = sum(reaction_evidence == -1, na.rm = TRUE),
    positive_active = NA_integer_, high_always_active = length(high_always_active_ids),
    negative_inactive = NA_integer_, low_always_inactive = length(low_always_inactive_ids),
    stringsAsFactors = FALSE
  )
  result <- list(model = bounded_model, scenario_model = bounded_model,
    scenario_exchange_bounds = bounded_model$effective_exchange_bounds,
    status = status, objective_value = NA_real_, incumbent_biomass_value = NA_real_, imat_objective = NA_real_,
    fluxes = stats::setNames(rep(NA_real_, n_flux), model$reaction_ids), evidence_qc = qc,
    reaction_evidence = reaction_evidence, settings = config, source_sha256 = model$source_sha256, medium = medium,
    baseline = baseline, max_abs_sv = NA_real_, solver = solver, solver_version = runtime_solver_version,
    solver_engine_version = as.character(runtime_engine_version), solver_code = solver_code,
    solver_status_message = status_message, incumbent_feasible = FALSE, solver_incumbent_objective = NA_real_,
    solver_timeout_enforced = timeout_enforced)
  if (!status %in% c("optimal", "time_limit") || length(solver_solution) != n_var || any(!is.finite(solver_solution))) return(result)
  candidate <- solver_solution
  activity <- as.numeric(constraints %*% candidate)
  feasibility_tolerance <- if (is.null(config$feasibility_tolerance)) 1e-7 else config$feasibility_tolerance
  feasible <- all(vapply(seq_along(direction), function(k) {
    scale <- max(1, abs(rhs[[k]]))
    if (direction[[k]] == "=") abs(activity[[k]] - rhs[[k]]) <= feasibility_tolerance * scale
    else if (direction[[k]] == "<=") activity[[k]] <= rhs[[k]] + feasibility_tolerance * scale
    else activity[[k]] >= rhs[[k]] - feasibility_tolerance * scale
  }, logical(1)))
  binary_deviation <- if (length(binary)) max(abs(candidate[binary] - round(candidate[binary]))) else 0
  feasible <- feasible && binary_deviation <= feasibility_tolerance && all(candidate[binary] >= -feasibility_tolerance & candidate[binary] <= 1 + feasibility_tolerance)
  if (!feasible) return(result)
  result$incumbent_feasible <- TRUE
  result$solver_incumbent_objective <- sum(objective * candidate)

  flux <- candidate[seq_len(n_flux)] + bounded_model$lb
  names(flux) <- model$reaction_ids
  integrated <- bounded_model
  pos <- if (n_positive) candidate[pos_start:(pos_start + n_positive - 1L)] else numeric()
  neg <- if (n_negative) candidate[neg_start:(neg_start + n_negative - 1L)] else numeric()
  inactive <- if (n_low) candidate[low_start:n_var] else numeric()
  positive_active_ids <- high_positive_ids[pos > 0.5]
  negative_active_ids <- high_negative_ids[neg > 0.5]
  # As in iMAT, only MILP-inactive low-evidence reactions constrain the context; active ones are not forced.
  inactive_ids <- low_ids[inactive > 0.5]
  if (length(inactive_ids)) {
    for (ridx in inactive_ids) {
      rid <- model$reaction_ids[[ridx]]
      integrated$lb[[rid]] <- max(integrated$lb[[rid]], -config$epsilon)
      integrated$ub[[rid]] <- min(integrated$ub[[rid]], config$epsilon)
    }
  }
  class(integrated) <- "session5_model"
  postfit <- solve_fba(integrated)
  result$model <- integrated
  result$incumbent_biomass_value <- sum(model$objective * flux)
  result$objective_value <- if (status == "optimal") result$incumbent_biomass_value else NA_real_
  solver_objective_value <- if (solver == "highs") fit$objective_value else fit$objval
  result$imat_objective <- if (status == "optimal") solver_objective_value else NA_real_
  result$fluxes <- flux
  result$max_abs_sv <- max(abs(as.numeric(model$S %*% flux)))
  result$evidence_qc$positive_active <- length(positive_active_ids) + length(negative_active_ids)
  result$evidence_qc$negative_inactive <- length(inactive_ids) + length(low_always_inactive_ids)
  result$reconstructed_fba_status <- postfit$status
  result$reconstructed_fba_objective_value <- postfit$objective_value
  result$reconstructed_fba_max_abs_sv <- postfit$max_abs_sv
  result
}
