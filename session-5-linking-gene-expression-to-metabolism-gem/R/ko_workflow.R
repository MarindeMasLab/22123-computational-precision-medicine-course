# Exercise 4 workflow. Source R/model_io.R first; run from exercises/.
ko_context_preflight <- function() {
  qc <- utils::read.csv("Exercise2_expression_integration/outputs/reconstruction_qc.csv", stringsAsFactors = FALSE)
  bounds <- utils::read.csv("Exercise2_expression_integration/outputs/effective_exchange_bounds.csv", stringsAsFactors = FALSE)
  manifest <- utils::read.csv("config/model_manifest.csv", stringsAsFactors = FALSE)
  settings <- read_imat_settings("config/imat_settings.yml")
  expected <- c(MCF7_ACH_000019_Jain = "Jain_NCI60_RPMI_recipe_proxy_v2",
    GTEx_breast_Keibler = "Keibler_184A1_EV_EGF5_growth_adjusted_v2")
  sha <- manifest$sha256[manifest$model_file == "RECON1_biomass.xml"]
  if (length(sha) != 1L || digest::digest(file = "datasets_and_models/models/RECON1_biomass.xml", algo = "sha256") != sha ||
      digest::digest(file = "config/imat_settings.yml", algo = "sha256") != unique(qc$settings_sha256)[[1L]] ||
      !setequal(qc$scenario_id, names(expected)) || anyDuplicated(qc$scenario_id)) stop("Source/settings/QC provenance mismatch.")
  contexts <- list(); rows <- list()
  for (id in names(expected)) {
    rds_path <- paste0("Exercise2_expression_integration/outputs/context_", id, ".rds")
    if (!file.exists(paste0("Exercise2_expression_integration/outputs/context_", id, ".xml"))) {
      stop("Missing matching canonical XML: ", id)
    }
    record <- readRDS(rds_path)
    rds_sha256 <- digest::digest(file = rds_path, algo = "sha256")
    model <- record$model; row <- qc[qc$scenario_id == id, , drop = FALSE]; p <- record$provenance
    if (nrow(row) != 1L || !inherits(model, "session5_model") ||
        !identical(p$scenario_id, id) || !identical(p$profile_id, row$profile_id) ||
        !identical(p$source_model_sha256, sha) || !identical(record$source_sha256, sha) ||
        !identical(model$source_sha256, sha) || !identical(row$source_model_sha256, sha) ||
        !identical(model$medium_id, expected[[id]]) || !identical(p$medium_id, model$medium_id) ||
        !identical(record$medium$id, model$medium_id) || !identical(row$medium_id, model$medium_id) ||
        !identical(record$settings, settings) || !identical(row$settings_id, settings$id) ||
        !identical(p$settings_sha256, row$settings_sha256) ||
        !identical(model$objective_id, "obj_biomass") || !identical(row$objective_id, model$objective_id) ||
        !identical(model$objective_sense, "maximize") || !identical(row$objective_sense, model$objective_sense) ||
        !identical(unname(model$objective[["R_biomass_reaction"]]), 1) || sum(model$objective != 0) != 1L ||
        !identical(record$status, "optimal") || !isTRUE(record$incumbent_feasible) ||
        !identical(record$reconstructed_fba_status, "optimal") || !isTRUE(row$incumbent_feasible) ||
        !identical(row$milp_status, "optimal") || !identical(row$reconstructed_fba_status, "optimal") ||
        !identical(record$solver, "highs") || !identical(record$solver_version, as.character(utils::packageVersion("highs"))) ||
        !identical(record$solver_engine_version, as.character(settings$solver_engine_version)) ||
        !identical(record$solver_version, row$solver_version)) stop("Reconstruction gate failed: ", id)
    for (field in c("evidence", "threshold", "crosswalk", "settings", "medium")) {
      file <- p[[paste0(field, "_file")]]; hash <- p[[paste0(field, "_sha256")]]
      if (length(file) != 1L || !file.exists(file) ||
          !identical(digest::digest(file = file, algo = "sha256"), hash)) stop("Input hash mismatch: ", id, "/", field)
    }
    b <- bounds[bounds$context == id, c("reaction_id", "lower_bound", "upper_bound")]
    m <- model$effective_exchange_bounds
    if (nrow(b) != nrow(m) || anyDuplicated(b$reaction_id) ||
        !identical(b$reaction_id, m$reaction_id) ||
        max(abs(b$lower_bound - m$lower_bound), abs(b$upper_bound - m$upper_bound)) > 1e-8) {
      stop("Effective exchange bounds differ: ", id)
    }
    wt <- solve_fba(model)
    if (wt$status != "optimal" || !is.finite(wt$objective_value) || wt$objective_value <= 0 ||
        !is.finite(wt$max_abs_sv) || wt$max_abs_sv > wt$tolerance ||
        abs(wt$objective_value - row$reconstructed_fba_objective_value) > 1e-6 ||
        !identical(wt$model_sha256, sha) || !identical(wt$medium_id, expected[[id]])) stop("Fresh WT FBA gate failed: ", id)
    contexts[[id]] <- list(model = model, wt = wt, reconstruction = record, rds_sha256 = rds_sha256)
    rows[[id]] <- data.frame(context = id, medium_id = model$medium_id,
      source_model_sha256 = sha, context_rds_sha256 = rds_sha256,
      settings_id = settings$id, settings_sha256 = p$settings_sha256,
      objective_id = model$objective_id, objective_sense = model$objective_sense,
      milp_status = record$status, incumbent_feasible = record$incumbent_feasible,
      reconstructed_fba_status = record$reconstructed_fba_status,
      wt_status = wt$status, wt_objective = wt$objective_value, wt_max_abs_sv = wt$max_abs_sv,
      solver = wt$solver, tolerance = wt$tolerance, stringsAsFactors = FALSE)
  }
  if (!identical(contexts[[1L]]$model$objective, contexts[[2L]]$model$objective) ||
      !identical(contexts[[1L]]$model$gpr, contexts[[2L]]$model$gpr) ||
      !identical(contexts[[1L]]$model$S, contexts[[2L]]$model$S) ||
      identical(contexts[[1L]]$model$medium_id, contexts[[2L]]$model$medium_id)) {
    stop("Source network, GPR, objective mismatch or native media unexpectedly identical.")
  }
  list(contexts = contexts, qc = do.call(rbind, rows))
}

# Input is a documented shortlist, NOT automatically an instructor-validated panel.
ko_panel_audit <- function(model, shortlist) {
  required <- c("gene_id", "source", "rationale")
  if (!is.data.frame(shortlist) || !all(required %in% names(shortlist)) ||
      anyNA(shortlist[, required]) || anyDuplicated(shortlist$gene_id)) stop("Shortlist needs unique gene IDs and source/rationale.")
  x <- shortlist[, required, drop = FALSE]
  x$gpr_reactions <- rep("", nrow(x)); x$affected_reactions <- rep("", nrow(x))
  x$coverage_count <- rep(0L, nrow(x))
  x$eligibility <- rep("excluded", nrow(x)); x$exclusion_reason <- rep("", nrow(x))
  for (i in seq_len(nrow(x))) {
    g <- x$gene_id[[i]]
    if (!nzchar(g) || !nzchar(x$source[[i]]) || !nzchar(x$rationale[[i]])) {
      x$exclusion_reason[[i]] <- "missing identifier/source/rationale"; next
    }
    if (!g %in% model$gene_ids) { x$exclusion_reason[[i]] <- "unknown model gene ID"; next }
    tokens <- lapply(model$gpr, function(rule) {
      if (!nzchar(rule)) return(character())
      setdiff(strsplit(trimws(gsub("[()]", " ", rule)), "\\s+")[[1L]], c("and", "or"))
    })
    covered <- names(tokens)[vapply(tokens, function(t) g %in% t, logical(1))]
    x$gpr_reactions[[i]] <- paste(covered, collapse = ";")
    x$coverage_count[[i]] <- length(covered)
    if (!length(covered)) { x$exclusion_reason[[i]] <- "no GPR coverage"; next }
    ko <- tryCatch(knock_out_genes(model, g), error = identity)
    if (inherits(ko, "error")) { x$exclusion_reason[[i]] <- conditionMessage(ko); next }
    x$affected_reactions[[i]] <- paste(ko$affected_reactions, collapse = ";")
    x$eligibility[[i]] <- "GPR-covered (not biologically validated)"
  }
  x
}

ko_screen <- function(contexts, shortlist, approved = FALSE) {
  if (!isTRUE(approved)) stop("Instructor review of shortlist and medium-confounded interpretation required before native screening.")
  audit <- ko_panel_audit(contexts[[1L]]$model, shortlist)
  rows <- list(); n <- 0L
  for (id in names(contexts)) {
    entry <- contexts[[id]]; wt <- solve_fba(entry$model)
    for (i in seq_len(nrow(audit))) {
      n <- n + 1L; g <- audit$gene_id[[i]]; excluded <- audit$eligibility[[i]] == "excluded"
      ko <- if (!excluded) tryCatch(knock_out_genes(entry$model, g,
        unavailable_gene_ids = entry$unavailable_gene_ids %||% character()), error = identity) else NULL
      failed <- inherits(ko, "error")
      fit <- if (!excluded && !failed) solve_fba(ko$model) else NULL
      status <- if (excluded) "excluded" else if (failed) "error" else fit$status
      value <- if (!is.null(fit)) fit$objective_value else NA_real_
      ratio <- if (wt$status == "optimal" && is.finite(wt$objective_value) && wt$objective_value > 0 &&
                   status == "optimal" && is.finite(value) && is.finite(fit$max_abs_sv) && fit$max_abs_sv <= fit$tolerance)
        value / wt$objective_value else NA_real_
      rows[[n]] <- data.frame(context = id, gene_id = g, source = audit$source[[i]],
        rationale = audit$rationale[[i]], gpr_reactions = audit$gpr_reactions[[i]],
        affected_reactions = if (excluded || failed) audit$affected_reactions[[i]] else paste(ko$affected_reactions, collapse = ";"),
        eligibility = if (excluded) "excluded" else if (failed) "unclassified" else if (is.finite(ratio)) "ratio available" else "unclassified",
        classification = if (is.finite(ratio)) "eligible for paired assessment" else "unclassified",
        reason = if (excluded) audit$exclusion_reason[[i]] else if (failed) conditionMessage(ko) else if (is.na(ratio)) status else "",
        deleted_genes = if (excluded || failed) "" else paste(ko$deleted_genes, collapse = ";"),
        wt_status = wt$status, ko_status = status, wt_objective = wt$objective_value,
        ko_objective = value, wt_max_abs_sv = wt$max_abs_sv,
        ko_max_abs_sv = if (!is.null(fit)) fit$max_abs_sv else NA_real_,
        ratio_to_wt = ratio, model_sha256 = entry$model$source_sha256,
        objective_id = entry$model$objective_id, objective_sense = entry$model$objective_sense,
        medium_id = entry$model$medium_id, solver = wt$solver, tolerance = wt$tolerance)
    }
  }
  if (!length(rows)) return(data.frame(context = character(), gene_id = character(),
    eligibility = character(), classification = character(), ko_status = character(), ratio_to_wt = numeric()))
  do.call(rbind, rows)
}

ko_selectivity <- function(rows, cancer, healthy, cancer_cutoff = 0.10, healthy_floor = 0.50) {
  eligible <- rows[!is.na(rows$ratio_to_wt) & rows$eligibility == "ratio available", , drop = FALSE]
  a <- eligible[eligible$context == cancer, c("gene_id", "ratio_to_wt")]
  b <- eligible[eligible$context == healthy, c("gene_id", "ratio_to_wt")]
  names(a)[2L] <- "cancer_ratio"; names(b)[2L] <- "healthy_ratio"
  paired <- merge(a, b, by = "gene_id")
  paired[paired$cancer_ratio < cancer_cutoff & paired$healthy_ratio >= healthy_floor, , drop = FALSE]
}

# Two parallel transporters with immutable rules; healthy has both routes,
# disease has the G_A-only route closed and G_A explicitly unavailable.
plot_toy_ko_network <- function(healthy, disease, path) {
  ids <- c("EX_toy_e", "T_toy_1", "T_toy_2", "BIOMASS_toy")
  if (!identical(healthy$S, disease$S) || !identical(healthy$gpr, disease$gpr) ||
      !identical(unname(healthy$gpr[ids]), c("", "G_A", "G_A or G_B", "G_C")) ||
      !identical(unname(as.numeric(healthy$S[, ids])), c(-1, 0, -1, 1, -1, 1, 0, -1)) ||
      !identical(disease$lb[["T_toy_1"]], 0) || !identical(disease$ub[["T_toy_1"]], 0))
    stop("Diagram requires the shared-GPR, two-branch synthetic network.")
  contexts <- list(healthy = list(model = healthy, unavailable = character()),
    disease = list(model = disease, unavailable = "G_A"))
  grDevices::png(path, width = 2200, height = 1270, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(mfrow = c(2, 4), mar = c(2.4, 1.2, 3.1, 1.2))
  on.exit(graphics::par(old), add = TRUE)
  for (name in names(contexts)) for (deleted in c("WT", "G_A", "G_B", "G_C")) {
    context <- contexts[[name]]
    ko <- if (deleted == "WT") NULL else knock_out_genes(context$model, deleted,
      unavailable_gene_ids = context$unavailable)
    model <- if (is.null(ko)) context$model else ko$model
    fit <- solve_fba(model); wt <- solve_fba(context$model)
    closed <- model$lb[ids] == 0 & model$ub[ids] == 0
    preclosed <- context$model$lb[ids] == 0 & context$model$ub[ids] == 0
    graphics::plot.new(); graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
    graphics::title(main = paste(name, deleted, sep = " | "), cex.main = 1.1)
    graphics::rect(.02, .08, .98, .98, border = "#cbd5df")
    graphics::points(c(.28, .76), c(.5, .5), pch = 21, bg = "#ecf4f7", cex = 1.7)
    graphics::text(c(.28, .76), c(.40, .40), c("substrate [e]", "precursor [c]"), cex = .77)
    paths <- list(c(.06, .5, .26, .5), c(.30, .5, .74, .5),
      c(.30, .5, .74, .5), c(.78, .5, .96, .5))
    labels <- c("EX / no GPR", "T1 / G_A", "T2 / G_A or G_B", "biomass / G_C")
    for (j in seq_along(ids)) {
      p <- paths[[j]]; color <- if (preclosed[[j]]) "#778899" else if (closed[[j]]) "#c44848" else "#197857"
      if (j %in% c(2L, 3L)) {
        bend <- if (j == 2L) .70 else .30
        graphics::segments(.30, .5, .51, bend, col = color, lwd = 3)
        graphics::arrows(.51, bend, .74, .5, col = color, lwd = 3, length = .09)
      } else graphics::arrows(p[1], p[2], p[3], p[4], col = color, lwd = 3, length = .09)
      if (closed[[j]]) graphics::text(if (j %in% c(2L, 3L)) .48 else mean(p[c(1, 3)]),
        if (j == 2L) .67 else if (j == 3L) .33 else .5, "X", col = color, cex = 1.5, font = 2)
      y <- if (j == 2L) .82 else if (j == 3L) .18 else .63
      graphics::text(mean(p[c(1, 3)]), y, labels[[j]], cex = .72, font = 2)
    }
    graphics::text(.5, .10, sprintf("unavailable: %s | KO/WT: %.1f",
      if (length(context$unavailable)) paste(context$unavailable, collapse = ",") else "none",
      fit$objective_value / wt$objective_value), cex = .74)
  }
  invisible(NULL)
}

# Exploratory census of all model gene products; NOT an instructor-reviewed
# dependency panel. Precompute exact GPR coverage, then solve only KOs that
# actually change a bound in that context. Every row has an auditable method.
ko_full_census <- function(verified, output_path = NULL, tolerance = 1e-7) {
  contexts <- verified$contexts
  if (length(contexts) != 2L || !identical(contexts[[1L]]$model$gene_ids, contexts[[2L]]$model$gene_ids))
    stop("Full census requires both verified contexts with identical gene IDs.")
  reference <- contexts[[1L]]$model
  gene_ids <- reference$gene_ids
  gpr_indexes <- lapply(contexts, function(entry) ko_gpr_index(entry$model))
  results <- vector("list", length(gene_ids) * length(contexts))
  index <- 0L
  for (id in names(contexts)) {
    m <- contexts[[id]]$model
    coverage <- gpr_indexes[[id]]$by_gene
    wt <- solve_fba(m, tolerance)
    if (wt$status != "optimal" || !is.finite(wt$objective_value) || wt$objective_value <= 0 ||
        wt$max_abs_sv > tolerance || abs(wt$objective_value - contexts[[id]]$wt$objective_value) > 1e-6)
      stop("Fresh WT gate failed during KO census: ", id)
    for (g in gene_ids) {
      index <- index + 1L
      covered <- coverage[[g]]
      ko <- if (length(covered)) tryCatch(knock_out_genes(m, g, index = gpr_indexes[[id]]), error = identity) else NULL
      error <- inherits(ko, "error")
      affected <- if (error || is.null(ko)) character() else ko$affected_reactions
      changed <- affected[m$lb[affected] != 0 | m$ub[affected] != 0]
      fit <- if (!error && length(changed)) solve_fba(ko$model, tolerance) else NULL
      status <- if (!length(covered)) "no_gpr_coverage" else if (error) "error" else
        if (!length(changed)) "optimal" else fit$status
      method <- if (!length(covered)) "not_tested" else if (error) "failed" else
        if (!length(changed)) "WT_reused_identical_bounds" else "HiGHS_KO_LP"
      objective <- if (status != "optimal") NA_real_ else if (!length(changed)) wt$objective_value else fit$objective_value
      residual <- if (status != "optimal") NA_real_ else if (!length(changed)) wt$max_abs_sv else fit$max_abs_sv
      valid <- status == "optimal" && is.finite(objective) && is.finite(residual) && residual <= tolerance
      results[[index]] <- data.frame(context = id, gene_id = g,
        gpr_coverage_count = length(covered), gpr_reactions = paste(covered, collapse = ";"),
        affected_count = length(affected), affected_reactions = paste(affected, collapse = ";"),
        changed_bound_count = length(changed), changed_bound_reactions = paste(changed, collapse = ";"),
        deleted_genes = if (length(covered) && !error) g else "", method = method,
        wt_status = wt$status, ko_status = status, wt_objective = wt$objective_value,
        ko_objective = objective, wt_max_abs_sv = wt$max_abs_sv, ko_max_abs_sv = residual,
        ratio_to_wt = if (valid) objective / wt$objective_value else NA_real_,
        eligibility = if (valid) "computational_only" else "unclassified",
        reason = if (!length(covered)) "gene not present in any GPR" else if (error) conditionMessage(ko) else
          if (!valid) status else if (!length(changed)) "no reaction bound changed" else "",
        medium_id = m$medium_id, objective_id = m$objective_id, objective_sense = m$objective_sense,
        model_sha256 = m$source_sha256,
        context_rds_sha256 = contexts[[id]]$rds_sha256 %||% NA_character_,
        settings_sha256 = contexts[[id]]$reconstruction$provenance$settings_sha256 %||% NA_character_,
        solver = wt$solver, tolerance = tolerance,
        stringsAsFactors = FALSE)
    }
    if (!identical(serialize(m, NULL), serialize(contexts[[id]]$model, NULL)))
      stop("WT model was mutated during census: ", id)
    if (!is.null(output_path)) {
      utils::write.csv(do.call(rbind, results[seq_len(index)]), output_path, row.names = FALSE)
    }
  }
  do.call(rbind, results)
}

# Native scenarios use different media: differential labels are *conditional*
# model contrasts, not disease-selective or expression-specific dependencies.
ko_compare_census <- function(rows, cancer = "MCF7_ACH_000019_Jain", healthy = "GTEx_breast_Keibler",
                              tolerance = 1e-6, cancer_cutoff = 0.10, healthy_floor = 0.50) {
  a <- rows[rows$context == cancer, , drop = FALSE]
  b <- rows[rows$context == healthy, , drop = FALSE]
  if (anyDuplicated(a$gene_id) || anyDuplicated(b$gene_id) || !setequal(a$gene_id, b$gene_id))
    stop("Paired census must have each gene exactly once in both contexts.")
  b <- b[match(a$gene_id, b$gene_id), , drop = FALSE]
  valid <- is.finite(a$ratio_to_wt) & is.finite(b$ratio_to_wt)
  out <- data.frame(gene_id = a$gene_id, cancer_ratio = a$ratio_to_wt, healthy_ratio = b$ratio_to_wt,
    cancer_status = a$ko_status, healthy_status = b$ko_status,
    cancer_method = a$method, healthy_method = b$method,
    cancer_affected = a$affected_reactions, healthy_affected = b$affected_reactions,
    cancer_changed_bounds = a$changed_bound_reactions, healthy_changed_bounds = b$changed_bound_reactions,
    paired_eligible = valid, delta_cancer_minus_healthy = ifelse(valid, a$ratio_to_wt - b$ratio_to_wt, NA_real_),
    classification = rep("unclassified", nrow(a)), stringsAsFactors = FALSE)
  out$classification[valid] <- "no_material_difference"
  out$classification[valid & abs(out$delta_cancer_minus_healthy) > tolerance] <- "different_native_scenarios"
  out$classification[valid & a$ratio_to_wt < cancer_cutoff & b$ratio_to_wt >= healthy_floor] <- "MCF7_lower_teaching_cutoffs"
  out$classification[valid & b$ratio_to_wt < cancer_cutoff & a$ratio_to_wt >= healthy_floor] <- "GTEx_lower_teaching_cutoffs"
  out$classification[a$ko_status == "no_gpr_coverage" | b$ko_status == "no_gpr_coverage"] <- "no_gpr_coverage"
  out
}

ko_pathway_analysis <- function(rows, comparisons, map_path = "config/recon1_subsystems.csv") {
  map <- utils::read.csv(map_path, stringsAsFactors = FALSE)
  if (anyDuplicated(map$reaction_id) || !all(c("reaction_id", "subsystem", "mapping_source") %in% names(map)))
    stop("Invalid subsystem crosswalk.")
  selected <- comparisons$gene_id[comparisons$classification %in% c(
    "MCF7_lower_teaching_cutoffs", "GTEx_lower_teaching_cutoffs", "different_native_scenarios")]
  subset <- rows[rows$gene_id %in% selected, , drop = FALSE]
  long <- lapply(seq_len(nrow(subset)), function(i) {
    ids <- strsplit(subset$changed_bound_reactions[[i]], ";", fixed = TRUE)[[1L]]
    ids <- ids[nzchar(ids)]
    if (!length(ids)) return(NULL)
    found <- match(ids, map$reaction_id)
    data.frame(context = subset$context[[i]], gene_id = subset$gene_id[[i]], reaction_id = ids,
      subsystem = ifelse(is.na(found) | is.na(map$subsystem[found]) | !nzchar(map$subsystem[found]),
        "Unassigned", map$subsystem[found]),
      mapping_source = ifelse(is.na(found), "unmapped", map$mapping_source[found]),
      ratio_to_wt = subset$ratio_to_wt[[i]], stringsAsFactors = FALSE)
  })
  long <- Filter(Negate(is.null), long)
  if (!length(long)) return(data.frame(context = character(), gene_id = character(), reaction_id = character(),
    subsystem = character(), mapping_source = character(), ratio_to_wt = numeric()))
  do.call(rbind, long)
}