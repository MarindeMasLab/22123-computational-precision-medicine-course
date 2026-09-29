# Exercise 5 workflow. Run from hands-on_session/ after sourcing R/model_io.R
# and R/ko_workflow.R.

exercise5_preflight <- function(
    ko_path = "Exercise4_disease_selective_knockouts/outputs/native_full_gene_ko.csv",
    comparison_path = "Exercise4_disease_selective_knockouts/outputs/native_full_gene_comparison.csv",
    hard_zero_path = NULL,
    manifest_path = "config/model_manifest.csv") {
  verified <- ko_context_preflight()
  native_ko <- utils::read.csv(ko_path, stringsAsFactors = FALSE, check.names = FALSE)
  native_comparison <- utils::read.csv(comparison_path, stringsAsFactors = FALSE, check.names = FALSE)
  hard_zero <- if (!is.null(hard_zero_path) && file.exists(hard_zero_path)) {
    utils::read.csv(hard_zero_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame(context = "GTEx_breast_Keibler", original_wt_biomass = NA_real_,
      variant_wt_status = "not_distributed", variant_wt_biomass = NA_real_,
      variant_status = "not_distributed", stringsAsFactors = FALSE)
  }
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  required_ko <- c("context", "gene_id", "gpr_coverage_count", "affected_reactions",
    "wt_status", "ko_status", "wt_objective", "ko_objective", "wt_max_abs_sv",
    "ko_max_abs_sv", "ratio_to_wt", "medium_id", "objective_id", "objective_sense",
    "model_sha256", "context_rds_sha256", "settings_sha256", "solver", "tolerance")
  required_comparison <- c("gene_id", "cancer_ratio", "healthy_ratio", "cancer_status",
    "healthy_status", "paired_eligible", "classification")
  if (!all(required_ko %in% names(native_ko)) || !all(required_comparison %in% names(native_comparison)))
    stop("Exercise 4 census schema mismatch.", call. = FALSE)
  if (nrow(native_ko) != 3810L || nrow(native_comparison) != 1905L ||
      anyDuplicated(native_comparison$gene_id)) stop("Exercise 4 census dimensions mismatch.", call. = FALSE)
  for (id in names(verified$contexts)) {
    rows <- native_ko[native_ko$context == id, , drop = FALSE]
    entry <- verified$contexts[[id]]
    if (nrow(rows) != 1905L || anyDuplicated(rows$gene_id) ||
        !identical(sort(rows$gene_id), sort(entry$model$gene_ids)) ||
        !identical(unique(rows$context_rds_sha256), entry$rds_sha256) ||
        !identical(unique(rows$model_sha256), entry$model$source_sha256) ||
        !identical(unique(rows$medium_id), entry$model$medium_id) ||
        !identical(unique(rows$objective_id), entry$model$objective_id) ||
        !identical(unique(rows$objective_sense), entry$model$objective_sense) ||
        !identical(unique(rows$settings_sha256), entry$reconstruction$provenance$settings_sha256) ||
        any(rows$wt_status != "optimal") || any(rows$ko_status != "optimal") ||
        any(!is.finite(rows$ratio_to_wt)) ||
        max(rows$wt_max_abs_sv, na.rm = TRUE) > unique(rows$tolerance) ||
        max(rows$ko_max_abs_sv, na.rm = TRUE) > unique(rows$tolerance) ||
        max(abs(rows$wt_objective - entry$wt$objective_value)) > 1e-10) {
      stop("Exercise 4 census provenance/QC mismatch: ", id, call. = FALSE)
    }
  }
  expected_source <- manifest$sha256[manifest$model_file == "RECON1_biomass.xml"]
  if (length(expected_source) != 1L ||
      !identical(expected_source, verified$contexts[[1L]]$model$source_sha256))
    stop("Current manifest does not match canonical contexts.", call. = FALSE)
  if (!is.null(hard_zero_path) && file.exists(hard_zero_path)) {
    hard_healthy <- hard_zero[hard_zero$context == "GTEx_breast_Keibler", , drop = FALSE]
    if (nrow(hard_healthy) != 1L || hard_healthy$variant_wt_status != "optimal" ||
        !identical(hard_healthy$variant_wt_biomass, 0) ||
        hard_healthy$variant_status != "blocked_nonpositive_or_failed_wt")
      stop("Hard-zero healthy gate no longer matches its quarantined result.", call. = FALSE)
  }
  list(verified = verified, native_ko = native_ko, native_comparison = native_comparison,
    hard_zero = hard_zero, manifest = manifest,
    manifest_sha256 = digest::digest(file = manifest_path, algo = "sha256"),
    native_ko_sha256 = digest::digest(file = ko_path, algo = "sha256"),
    native_comparison_sha256 = digest::digest(file = comparison_path, algo = "sha256"))
}

exercise5_safe_ratio <- function(wt, knockout, tolerance = 1e-7) {
  valid_wt <- identical(wt$status, "optimal") && is.finite(wt$objective_value) &&
    wt$objective_value > 0 && is.finite(wt$max_abs_sv) && wt$max_abs_sv <= tolerance
  valid_ko <- identical(knockout$status, "optimal") && is.finite(knockout$objective_value) &&
    is.finite(knockout$max_abs_sv) && knockout$max_abs_sv <= tolerance
  if (!valid_wt || !valid_ko) return(NA_real_)
  knockout$objective_value / wt$objective_value
}

exercise5_toy_model <- function() {
  S <- Matrix::sparseMatrix(i = c(1, 1, 2, 1, 2, 2), j = c(1, 2, 2, 3, 3, 4),
    x = c(-1, -1, 1, -1, 1, -1), dims = c(2, 4),
    dimnames = list(c("substrate_e", "precursor_c"),
      c("EX_substrate_e", "PATH_A", "PATH_B", "BIOMASS")))
  model <- list(S = S,
    lb = c(EX_substrate_e = -10, PATH_A = 0, PATH_B = 0, BIOMASS = 0),
    ub = c(EX_substrate_e = 0, PATH_A = 10, PATH_B = 10, BIOMASS = 100),
    objective = c(EX_substrate_e = 0, PATH_A = 0, PATH_B = 0, BIOMASS = 1),
    objective_id = "toy_biomass", objective_sense = "maximize",
    gpr = c(EX_substrate_e = "", PATH_A = "G_A and (G_H1 or G_H2)",
      PATH_B = "G_B", BIOMASS = "G_C"),
    reaction_ids = c("EX_substrate_e", "PATH_A", "PATH_B", "BIOMASS"),
    metabolite_ids = c("substrate_e", "precursor_c"),
    gene_ids = c("G_A", "G_B", "G_C", "G_H1", "G_H2"),
    exchange_ids = "EX_substrate_e", model_path = "synthetic_fixture",
    source_sha256 = "synthetic_two_route_v2", flux_unit = "arbitrary_fixture_units",
    medium_id = "synthetic_two_route_medium_v2", medium_source = "invented teaching fixture",
    medium_config = NULL, compartments = c(substrate_e = "e", precursor_c = "c"))
  class(model) <- "session5_model"
  model
}

.exercise5_fit_row <- function(model, wt, label, deletion_type, genes, tolerance = 1e-7) {
  if (!length(genes)) {
    fit <- wt
    affected <- character()
  } else {
    deleted <- knock_out_genes(model, genes)
    fit <- solve_fba(deleted$model, tolerance)
    affected <- deleted$affected_reactions
  }
  ratio <- exercise5_safe_ratio(wt, fit, tolerance)
  data.frame(context = "synthetic_two_route", condition = label,
    deletion_type = deletion_type, deleted_genes = paste(genes, collapse = ";"),
    affected_reactions = paste(affected, collapse = ";"), status = fit$status,
    objective_value = fit$objective_value, max_abs_sv = fit$max_abs_sv,
    ratio_to_wt = ratio, wt_status = wt$status, wt_objective = wt$objective_value,
    wt_max_abs_sv = wt$max_abs_sv, objective_id = model$objective_id,
    objective_sense = model$objective_sense, medium_id = model$medium_id,
    model_sha256 = model$source_sha256, solver = fit$solver, tolerance = tolerance,
    evidence_type = "synthetic software fixture; not a native model prediction",
    stringsAsFactors = FALSE)
}

exercise5_screen_pairs <- function(model, pairs, single_rows, context,
    context_rds_sha256 = NA_character_, tolerance = 1e-7,
    single_floor = 0.50, pair_cutoff = 0.10) {
  pair_fields <- c("pair_id", "gene_a", "gene_b", "selection_source", "selection_rationale")
  single_fields <- c("gene_id", "status", "objective_value", "max_abs_sv", "ratio_to_wt",
    "affected_reactions", "objective_id", "medium_id", "model_sha256", "context_rds_sha256")
  if (!is.data.frame(pairs) || !all(pair_fields %in% names(pairs)) ||
      anyNA(pairs[, pair_fields]) || anyDuplicated(pairs$pair_id) ||
      any(!nzchar(pairs$selection_source)) || any(!nzchar(pairs$selection_rationale)) ||
      any(pairs$gene_a == pairs$gene_b)) stop("Pairs require unique IDs, two genes, and predeclared source/rationale.", call. = FALSE)
  if (!is.data.frame(single_rows) || !all(single_fields %in% names(single_rows)) ||
      anyDuplicated(single_rows$gene_id)) stop("Single-control table schema mismatch.", call. = FALSE)
  genes <- unique(c(pairs$gene_a, pairs$gene_b))
  if (any(!genes %in% model$gene_ids)) stop("Pair contains an unknown model gene ID.", call. = FALSE)
  if (any(!genes %in% single_rows$gene_id)) stop("Pair is missing a single-KO control.", call. = FALSE)
  if (any(single_rows$objective_id != model$objective_id) ||
      any(single_rows$medium_id != model$medium_id) ||
      any(single_rows$model_sha256 != model$source_sha256) ||
      any(single_rows$context_rds_sha256 != context_rds_sha256))
    stop("Single-control provenance does not match the pair context.", call. = FALSE)
  wt <- solve_fba(model, tolerance)
  original <- serialize(model, NULL)
  index <- ko_gpr_index(model)
  output <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    a <- single_rows[single_rows$gene_id == pairs$gene_a[[i]], , drop = FALSE]
    b <- single_rows[single_rows$gene_id == pairs$gene_b[[i]], , drop = FALSE]
    valid_single <- function(row) row$status == "optimal" && is.finite(row$objective_value) &&
      is.finite(row$max_abs_sv) && row$max_abs_sv <= tolerance &&
      is.finite(row$ratio_to_wt) && row$ratio_to_wt >= single_floor
    eligible <- valid_single(a) && valid_single(b)
    deletion <- fit <- NULL
    if (eligible) {
      deletion <- tryCatch(knock_out_genes(model, c(pairs$gene_a[[i]], pairs$gene_b[[i]]),
        index = index), error = identity)
      fit <- if (!inherits(deletion, "error")) solve_fba(deletion$model, tolerance) else NULL
    }
    pair_status <- if (!eligible) "not_run" else if (inherits(deletion, "error")) "error" else fit$status
    pair_ratio <- if (!is.null(fit)) exercise5_safe_ratio(wt, fit, tolerance) else NA_real_
    classification <- if (!eligible) "excluded_single_control" else if (is.na(pair_ratio))
      "unclassified" else if (pair_ratio < pair_cutoff) "teaching_synthetic_lethal" else "tested_not_hit"
    output[[i]] <- data.frame(pair_id = pairs$pair_id[[i]], gene_a = pairs$gene_a[[i]],
      gene_b = pairs$gene_b[[i]], selection_source = pairs$selection_source[[i]],
      selection_rationale = pairs$selection_rationale[[i]],
      single_a_status = a$status, single_a_objective = a$objective_value,
      single_a_max_abs_sv = a$max_abs_sv, single_a_ratio = a$ratio_to_wt,
      single_a_affected_reactions = a$affected_reactions,
      single_b_status = b$status, single_b_objective = b$objective_value,
      single_b_max_abs_sv = b$max_abs_sv, single_b_ratio = b$ratio_to_wt,
      single_b_affected_reactions = b$affected_reactions,
      eligible_single_controls = eligible, pair_status = pair_status,
      pair_objective = if (!is.null(fit)) fit$objective_value else NA_real_,
      pair_max_abs_sv = if (!is.null(fit)) fit$max_abs_sv else NA_real_, pair_ratio = pair_ratio,
      pair_affected_reactions = if (!is.null(deletion) && !inherits(deletion, "error"))
        paste(deletion$affected_reactions, collapse = ";") else "",
      classification = classification, context = context, wt_status = wt$status,
      wt_objective = wt$objective_value, wt_max_abs_sv = wt$max_abs_sv,
      objective_id = model$objective_id, objective_sense = model$objective_sense,
      medium_id = model$medium_id, model_sha256 = model$source_sha256,
      context_rds_sha256 = context_rds_sha256, solver = wt$solver, tolerance = tolerance,
      stringsAsFactors = FALSE)
  }
  if (!identical(original, serialize(model, NULL))) stop("WT model was mutated during pair screen.", call. = FALSE)
  do.call(rbind, output)
}

exercise5_native_panel <- function(preflight, panel_size = NULL, single_floor = 0.50) {
  if (!is.null(panel_size) && (length(panel_size) != 1L || !is.finite(panel_size) || panel_size < 1L))
    stop("panel_size must be NULL or one positive integer.", call. = FALSE)
  cancer_id <- "MCF7_ACH_000019_Jain"
  model <- preflight$verified$contexts[[cancer_id]]$model
  controls <- preflight$native_ko[preflight$native_ko$context == cancer_id &
    preflight$native_ko$ko_status == "optimal" &
    is.finite(preflight$native_ko$ratio_to_wt) &
    preflight$native_ko$ratio_to_wt >= single_floor, , drop = FALSE]
  pure_or <- model$gpr[nzchar(model$gpr) & grepl(" or ", model$gpr, fixed = TRUE) &
    !grepl(" and ", model$gpr, fixed = TRUE)]
  candidates <- lapply(names(pure_or), function(reaction_id) {
    genes <- unique(setdiff(strsplit(trimws(gsub("[()]", " ", pure_or[[reaction_id]])), "\\s+")[[1L]],
      c("and", "or")))
    if (length(genes) != 2L || any(!genes %in% controls$gene_id)) return(NULL)
    genes <- sort(genes)
    data.frame(gene_a = genes[[1L]], gene_b = genes[[2L]], reaction_id = reaction_id,
      gpr = pure_or[[reaction_id]], stringsAsFactors = FALSE)
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) stop("No exact two-gene OR candidates pass the single-KO gate.", call. = FALSE)
  candidates <- unique(do.call(rbind, candidates))
  groups <- split(seq_len(nrow(candidates)), paste(candidates$gene_a, candidates$gene_b, sep = "|"))
  panel <- do.call(rbind, lapply(groups, function(rows) {
    reaction_ids <- sort(unique(candidates$reaction_id[rows]))
    rules <- unique(candidates$gpr[rows])
    data.frame(pair_id = paste(candidates$gene_a[[rows[[1L]]]], candidates$gene_b[[rows[[1L]]]], sep = "__"),
      gene_a = candidates$gene_a[[rows[[1L]]]], gene_b = candidates$gene_b[[rows[[1L]]]],
      selection_source = "Recon1 exact two-gene pure-OR GPR topology plus Exercise 4 MCF7 single-KO controls",
      selection_rationale = paste0(length(reaction_ids),
        " shared pure-OR reaction(s); both MCF7 single ratios >= ", format(single_floor, nsmall = 2)),
      shared_or_reaction_count = length(reaction_ids),
      shared_or_reactions = paste(reaction_ids, collapse = ";"),
      shared_or_gprs = paste(rules, collapse = ";"), stringsAsFactors = FALSE)
  }))
  panel <- panel[order(-panel$shared_or_reaction_count, panel$gene_a, panel$gene_b), , drop = FALSE]
  panel$selection_rank <- seq_len(nrow(panel))
  selected_count <- if (is.null(panel_size)) nrow(panel) else min(as.integer(panel_size), nrow(panel))
  panel$selected <- panel$selection_rank <= selected_count
  panel$exclusion_reason <- ifelse(panel$selected, "", "outside requested panel_size")
  panel
}

.exercise5_native_single_rows <- function(preflight, context) {
  x <- preflight$native_ko[preflight$native_ko$context == context, , drop = FALSE]
  data.frame(gene_id = x$gene_id, status = x$ko_status, objective_value = x$ko_objective,
    max_abs_sv = x$ko_max_abs_sv, ratio_to_wt = x$ratio_to_wt,
    affected_reactions = x$affected_reactions, objective_id = x$objective_id,
    medium_id = x$medium_id, model_sha256 = x$model_sha256,
    context_rds_sha256 = x$context_rds_sha256, stringsAsFactors = FALSE)
}

exercise5_run_native_panel <- function(preflight, panel_size = NULL, tolerance = 1e-7) {
  audit <- exercise5_native_panel(preflight, panel_size = panel_size)
  selected <- audit[audit$selected, c("pair_id", "gene_a", "gene_b", "selection_source",
    "selection_rationale"), drop = FALSE]
  contexts <- c("MCF7_ACH_000019_Jain", "GTEx_breast_Keibler")
  results <- lapply(contexts, function(context) {
    entry <- preflight$verified$contexts[[context]]
    exercise5_screen_pairs(entry$model, selected,
      .exercise5_native_single_rows(preflight, context), context,
      context_rds_sha256 = entry$rds_sha256, tolerance = tolerance,
      single_floor = if (context == "MCF7_ACH_000019_Jain") 0.50 else -Inf)
  })
  results <- do.call(rbind, results)
  cancer <- results[results$context == contexts[[1L]], , drop = FALSE]
  healthy <- results[results$context == contexts[[2L]], , drop = FALSE]
  healthy <- healthy[match(cancer$pair_id, healthy$pair_id), , drop = FALSE]
  contrast <- data.frame(pair_id = cancer$pair_id, gene_a = cancer$gene_a, gene_b = cancer$gene_b,
    selection_source = cancer$selection_source, selection_rationale = cancer$selection_rationale,
    mcf7_single_a_ratio = cancer$single_a_ratio, mcf7_single_b_ratio = cancer$single_b_ratio,
    mcf7_pair_status = cancer$pair_status, mcf7_pair_ratio = cancer$pair_ratio,
    gtex_single_a_ratio = healthy$single_a_ratio, gtex_single_b_ratio = healthy$single_b_ratio,
    gtex_pair_status = healthy$pair_status, gtex_pair_ratio = healthy$pair_ratio,
    scenario_conditional_classification = ifelse(
      cancer$pair_status == "optimal" & is.finite(cancer$pair_ratio) & cancer$pair_ratio < 0.10 &
      healthy$pair_status == "optimal" & is.finite(healthy$pair_ratio) & healthy$pair_ratio >= 0.50,
      "MCF7_low_GTEx_retained_teaching_cutoffs",
      ifelse(cancer$pair_status == "optimal" & is.finite(cancer$pair_ratio) & cancer$pair_ratio < 0.10,
        "MCF7_low_not_GTEx_retained", ifelse(is.finite(cancer$pair_ratio), "tested_not_hit", "unclassified"))),
    interpretation_scope = "exploratory native scenarios; unmatched profiles and different medium proxies",
    stringsAsFactors = FALSE)
  list(audit = audit, results = results, contrast = contrast)
}

exercise5_run_synthetic <- function(tolerance = 1e-7) {
  model <- exercise5_toy_model()
  original <- serialize(model, NULL)
  wt <- solve_fba(model, tolerance)
  conditions <- list(WT = character(), single_G_A = "G_A", single_G_B = "G_B",
    single_G_C = "G_C", pair_G_A_G_B = c("G_A", "G_B"))
  types <- c("WT", "single", "single", "single", "simultaneous_pair")
  rows <- do.call(rbind, Map(function(label, genes, type)
    .exercise5_fit_row(model, wt, label, type, genes, tolerance),
    names(conditions), conditions, types))
  if (!identical(original, serialize(model, NULL))) stop("Synthetic WT was mutated.", call. = FALSE)
  get_row <- function(label) rows[rows$condition == label, , drop = FALSE]
  singles <- do.call(rbind, lapply(c(G_A = "single_G_A", G_B = "single_G_B", G_C = "single_G_C"),
    function(label) get_row(label)))
  singles$gene_id <- c("G_A", "G_B", "G_C")
  singles$context_rds_sha256 <- "synthetic_not_rds"
  panel <- data.frame(pair_id = c("toy_parallel_routes", "toy_contains_essential_single"),
    gene_a = c("G_A", "G_A"), gene_b = c("G_B", "G_C"),
    selection_source = "predeclared synthetic teaching fixture",
    selection_rationale = c("genes support separate alternative routes",
      "negative control containing an essential single"),
    stringsAsFactors = FALSE)
  result <- exercise5_screen_pairs(model, panel, singles, "synthetic_two_route",
    context_rds_sha256 = "synthetic_not_rds", tolerance = tolerance)
  result$classification[result$pair_id == "toy_contains_essential_single"] <-
    "excluded_already_essential_single"
  result$evidence_type <- "synthetic fixture only"
  list(model = model, rows = rows, pairs = result)
}

exercise5_write_outputs <- function(output_dir = "Exercise5_synthetic_lethality/outputs") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stale_solutions <- file.path(output_dir, c("native_pair_context_results.csv",
    "native_pair_contrasts.csv", "native_pair_context_ratios.png"))
  unlink(stale_solutions[file.exists(stale_solutions)])
  preflight <- exercise5_preflight()
  synthetic <- exercise5_run_synthetic()
  native_audit <- exercise5_native_panel(preflight)
  context_qc <- preflight$verified$qc
  context_qc$manifest_sha256 <- preflight$manifest_sha256
  context_qc$native_ko_sha256 <- preflight$native_ko_sha256
  context_qc$native_comparison_sha256 <- preflight$native_comparison_sha256
  context_qc$comparison_scope <- c(
    "MCF7 expression plus RPMI recipe-proxy assumptions",
    "GTEx breast expression plus Keibler/DMEM proxy assumptions")
  utils::write.csv(context_qc, file.path(output_dir, "canonical_context_preflight.csv"), row.names = FALSE)
  utils::write.csv(synthetic$rows, file.path(output_dir, "synthetic_wt_single_pair.csv"), row.names = FALSE)
  utils::write.csv(synthetic$pairs, file.path(output_dir, "synthetic_pair_classification.csv"), row.names = FALSE)
  utils::write.csv(native_audit, file.path(output_dir, "native_pair_panel_audit.csv"), row.names = FALSE)
  panel_status <- data.frame(status = "student_live_screen_not_precomputed",
    pair_count = sum(native_audit$selected), double_ko_runs = 0L,
    reason = "All exact two-gene pure-OR GPR pairs after valid MCF7 single-KO ratio >= 0.50; students run both contexts live and no native solution checkpoint is distributed.",
    native_ko_rows_verified = nrow(preflight$native_ko),
    native_comparison_rows_verified = nrow(preflight$native_comparison),
    native_teaching_single_selectivity_hits = sum(preflight$native_comparison$classification == "MCF7_lower_teaching_cutoffs"),
    hard_zero_healthy_wt = preflight$hard_zero$variant_wt_biomass[preflight$hard_zero$context == "GTEx_breast_Keibler"],
    hard_zero_use = "excluded_from_default_workflow", stringsAsFactors = FALSE)
  utils::write.csv(panel_status, file.path(output_dir, "native_pair_panel_status.csv"), row.names = FALSE)
  plot_exercise5_synthetic_network(synthetic$model, file.path(output_dir, "synthetic_pair_network.png"))
  plot_exercise5_synthetic_biomass(synthetic$rows, file.path(output_dir, "synthetic_pair_biomass.png"))
  invisible(list(preflight = preflight, synthetic = synthetic, native_audit = native_audit,
    panel_status = panel_status))
}

plot_exercise5_native_pairs <- function(contrast, path) {
  order_index <- order(contrast$mcf7_pair_ratio, contrast$gtex_pair_ratio, contrast$pair_id,
    na.last = TRUE)
  x <- seq_along(order_index)
  grDevices::png(path, width = 1800, height = 1050, res = 170)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(5, 5, 4, 2))
  graphics::plot(x, contrast$mcf7_pair_ratio[order_index], pch = 16, cex = .72,
    col = "#b43a3a", ylim = c(0, 1.05), xlab = "Pair rank after sorting by MCF7 ratio",
    ylab = "Pair biomass objective / context-specific WT",
    main = "Full restricted native pair screen (170 pairs × 2 saved scenarios)")
  graphics::points(x, contrast$gtex_pair_ratio[order_index], pch = 1, cex = .8,
    col = "#147a63", lwd = 1.4)
  graphics::abline(h = c(.10, .50), lty = 2, col = c("#b43a3a", "#59636e"))
  graphics::legend("bottomright", c("MCF7 + RPMI proxy", "GTEx/Keibler + DMEM proxy"),
    pch = c(16, 1), col = c("#b43a3a", "#147a63"), bty = "n")
  graphics::mtext("Model predictions under different media and unmatched profiles; not observed growth",
    side = 1, line = 3.4, cex = .8)
}

plot_exercise5_synthetic_network <- function(model, path) {
  states <- list(WT = character(), `single G_A` = "G_A", `single G_B` = "G_B",
    `simultaneous G_A + G_B` = c("G_A", "G_B"))
  grDevices::png(path, width = 2200, height = 900, res = 170)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(mfrow = c(1, 4), mar = c(2, 1, 3, 1))
  on.exit(graphics::par(old), add = TRUE)
  for (label in names(states)) {
    genes <- states[[label]]
    changed <- if (length(genes)) knock_out_genes(model, genes)$affected_reactions else character()
    graphics::plot.new(); graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
    graphics::title(label, cex.main = .95)
    graphics::points(c(.18, .82), c(.5, .5), pch = 21, bg = "#eef4f2", cex = 1.8)
    graphics::text(c(.18, .82), c(.39, .39), c("substrate", "precursor"), cex = .72)
    for (route in c("PATH_A", "PATH_B")) {
      upper <- route == "PATH_A"; y <- if (upper) .70 else .30
      color <- if (route %in% changed) "#b43a3a" else "#147a63"
      graphics::segments(.21, .5, .48, y, col = color, lwd = 3)
      graphics::arrows(.48, y, .79, .5, col = color, lwd = 3, length = .09)
      graphics::text(.5, y + if (upper) .10 else -.10,
        if (upper) "PATH_A\nG_A and (G_H1 or G_H2)" else "PATH_B\nG_B", cex = .62)
      if (route %in% changed) graphics::text(.49, y, "X", col = color, cex = 1.5, font = 2)
    }
    graphics::arrows(.85, .5, .97, .5, col = "#46515d", lwd = 3, length = .09)
    graphics::text(.9, .64, "BIOMASS\nG_C", cex = .65)
  }
}

plot_exercise5_synthetic_biomass <- function(rows, path) {
  shown <- rows[rows$condition %in% c("WT", "single_G_A", "single_G_B", "pair_G_A_G_B"), ]
  labels <- c("WT", "single G_A", "single G_B", "G_A + G_B")
  grDevices::png(path, width = 1300, height = 850, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  colors <- c("#44556b", "#2a7f78", "#2a7f78", "#b43a3a")
  bars <- graphics::barplot(shown$ratio_to_wt, names.arg = labels, col = colors,
    ylim = c(0, 1.12), ylab = "Biomass objective ratio to synthetic WT",
    main = "Synthetic fixture: singles tolerated, simultaneous pair disruptive")
  graphics::abline(h = c(.10, .50), lty = 2, col = c("#b43a3a", "#59636e"))
  graphics::text(bars, shown$ratio_to_wt + .05, labels = format(shown$ratio_to_wt, digits = 2))
}