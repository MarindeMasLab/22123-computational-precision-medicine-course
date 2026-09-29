exercise6_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

exercise6_hash <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("The digest package is required.", call. = FALSE)
  digest::digest(file = path, algo = "sha256")
}

exercise6_read_settings <- function(path) {
  exercise6_assert(file.exists(path), paste("Missing validation settings:", path))
  lines <- trimws(readLines(path, warn = FALSE))
  get_scalar <- function(key) {
    hit <- grep(paste0("^", key, "\\s*:"), lines, value = TRUE)
    exercise6_assert(length(hit) == 1L, paste("Settings key is missing or duplicated:", key))
    value <- trimws(sub(paste0("^", key, "\\s*:\\s*"), "", hit))
    sub("^['\"]|['\"]$", "", value)
  }
  out <- list(
    version = get_scalar("version"),
    observation_unit = get_scalar("observation_unit"),
    model_unit = get_scalar("model_unit"),
    observation_zero_tolerance = as.numeric(get_scalar("observation_zero_tolerance")),
    model_zero_tolerance = as.numeric(get_scalar("model_zero_tolerance")),
    replicate_policy = get_scalar("replicate_policy"),
    seed = as.integer(get_scalar("seed"))
  )
  exercise6_assert(is.finite(out$observation_zero_tolerance) && out$observation_zero_tolerance >= 0,
    "Observation tolerance must be finite and non-negative.")
  exercise6_assert(is.finite(out$model_zero_tolerance) && out$model_zero_tolerance >= 0,
    "Model tolerance must be finite and non-negative.")
  exercise6_assert(out$observation_unit == "fmol/cell/h" && out$model_unit == "mmol/gDW/h",
    "Observation and model units must remain separate.")
  out
}

exercise6_normalize_exchange <- function(ids) {
  ids <- as.character(ids)
  ifelse(grepl("^R_", ids), ids, paste0("R_", ids))
}

exercise6_load_preflight <- function(root, settings) {
  model_path <- file.path(root, "Exercise2_expression_integration/outputs/context_MCF7_ACH_000019_Jain.rds")
  healthy_path <- file.path(root, "Exercise2_expression_integration/outputs/context_GTEx_breast_Keibler.rds")
  fva_path <- file.path(root, "Exercise3_pathway_and_flux_comparison/outputs/native_context_fva_comparison.csv")
  fva_provenance_path <- file.path(root, "Exercise3_pathway_and_flux_comparison/outputs/native_context_fva_provenance.csv")
  qc_path <- file.path(root, "Exercise2_expression_integration/outputs/reconstruction_qc.csv")
  bounds_path <- file.path(root, "Exercise2_expression_integration/outputs/effective_exchange_bounds.csv")
  required <- c(model_path, healthy_path, fva_path, fva_provenance_path, qc_path, bounds_path)
  exercise6_assert(all(file.exists(required)), "Canonical context/FVA/QC artifact is missing.")

  context <- readRDS(model_path)
  healthy <- readRDS(healthy_path)
  fva <- utils::read.csv(fva_path, stringsAsFactors = FALSE, check.names = FALSE)
  fva_provenance <- utils::read.csv(fva_provenance_path, stringsAsFactors = FALSE, check.names = FALSE)
  qc <- utils::read.csv(qc_path, stringsAsFactors = FALSE, check.names = FALSE)
  bounds <- utils::read.csv(bounds_path, stringsAsFactors = FALSE, check.names = FALSE)
  expected_hash <- "fd87a09e24e5b61e1b3c9c69178484792735a55123a32f5dec57be08cb66fe17"
  exercise6_assert(identical(context$source_sha256, expected_hash) && identical(healthy$source_sha256, expected_hash),
    "Canonical source model hash mismatch.")
  exercise6_assert(context$provenance$scenario_id == "MCF7_ACH_000019_Jain" &&
    healthy$provenance$scenario_id == "GTEx_breast_Keibler", "Canonical context scenario mismatch.")
  exercise6_assert(context$provenance$objective_id == "obj_biomass" && context$provenance$objective_sense == "maximize",
    "MCF7 objective provenance mismatch.")
  exercise6_assert(context$status == "optimal" && context$reconstructed_fba_status == "optimal" &&
    isTRUE(context$incumbent_feasible) && context$reconstructed_fba_objective_value > 0 &&
    context$reconstructed_fba_max_abs_sv < 1e-6, "MCF7 WT preflight is not optimal, positive, and residual-valid.")
  exercise6_assert(healthy$status == "optimal" && healthy$reconstructed_fba_status == "optimal" &&
    isTRUE(healthy$incumbent_feasible) && healthy$reconstructed_fba_objective_value > 0 &&
    healthy$reconstructed_fba_max_abs_sv < 1e-6, "GTEx WT preflight is not optimal, positive, and residual-valid.")
  exercise6_assert(all(fva$model_sha256 == expected_hash) && all(fva_provenance$model_sha256 == expected_hash),
    "FVA source model hash mismatch.")
  exercise6_assert(all(fva$objective_id == "obj_biomass") && all(fva$objective_sense == "maximize") &&
    all(fva$fva_fraction == 0.9) && all(fva$fva_min_status == "optimal") && all(fva$fva_max_status == "optimal"),
    "FVA status/objective/fraction provenance mismatch.")
  exercise6_assert(all(c("R_EX_ala__L_e", "R_EX_glc__D_e") %in% fva$reaction_id), "Expected FVA exchange panel is incomplete.")
  exercise6_assert(all(qc$source_model_sha256 == expected_hash) && all(bounds$reaction_id != ""), "Upstream QC/bounds provenance mismatch.")

  preflight <- data.frame(
    artifact = c("MCF7 context", "GTEx context", "Exercise 3 FVA comparison", "Exercise 3 FVA provenance", "reconstruction QC", "effective exchange bounds", "validation settings"),
    path = c(model_path, healthy_path, fva_path, fva_provenance_path, qc_path, bounds_path, settings$path),
    sha256 = c(exercise6_hash(model_path), exercise6_hash(healthy_path), exercise6_hash(fva_path), exercise6_hash(fva_provenance_path), exercise6_hash(qc_path), exercise6_hash(bounds_path), exercise6_hash(settings$path)),
    status = "verified", stringsAsFactors = FALSE
  )
  list(context = context, healthy = healthy, fva = fva, fva_provenance = fva_provenance,
    qc = qc, bounds = bounds, preflight = preflight)
}

exercise6_read_observations <- function(root, settings) {
  source_path <- file.path(root, "datasets_and_models/treated/cancer/metabolomics/Jain2012_Database_S1.csv")
  crosswalk_path <- file.path(root, "datasets_and_models/treated/cancer/transcriptomics/DepMap_Jain_NCI60_crosswalk.csv")
  source <- utils::read.csv(source_path, stringsAsFactors = FALSE, check.names = FALSE)
  crosswalk <- utils::read.csv(crosswalk_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_source <- c("cell_line", "replicate", "method", "metabolite", "calibrated", "observed_rate", "rate_unit", "exchange_reaction_id", "rate_in_cobra_exchange_sign", "mapping_status", "depmap_sample_id")
  exercise6_assert(all(required_source %in% names(source)), "Jain source schema is incomplete.")
  source$source_row_id <- with(source, paste("JAIN", cell_line, replicate, method, metabolite, sep = "::"))
  exercise6_assert(!anyDuplicated(source$source_row_id), "Jain source surrogate row IDs are not unique.")
  exercise6_assert(any(crosswalk$depmap_sample_id == "ACH_000019"), "MCF7 ACH_000019 crosswalk match is missing.")
  eligible <- source$depmap_sample_id == "ACH_000019" & source$calibrated == 1 &
    is.finite(source$observed_rate) & source$rate_unit == settings$observation_unit &
    !is.na(source$rate_in_cobra_exchange_sign) & is.finite(source$rate_in_cobra_exchange_sign) &
    !is.na(source$exchange_reaction_id) & nzchar(source$exchange_reaction_id) &
    source$mapping_status %in% c("exact_name_or_id", "curated_alias")
  observations <- source[eligible, , drop = FALSE]
  observations$normalized_exchange_id <- exercise6_normalize_exchange(observations$exchange_reaction_id)
  observations$replicate_direction <- ifelse(observations$rate_in_cobra_exchange_sign < -settings$observation_zero_tolerance, "uptake",
    ifelse(observations$rate_in_cobra_exchange_sign > settings$observation_zero_tolerance, "secretion", "near_zero_or_ambiguous"))
  group_ids <- unique(observations$normalized_exchange_id)
  group_label <- vapply(group_ids, function(id) {
    signs <- observations$replicate_direction[observations$normalized_exchange_id == id]
    if (length(signs) == 0L || any(signs == "near_zero_or_ambiguous") || length(unique(signs)) != 1L) "near_zero_or_ambiguous" else signs[[1L]]
  }, character(1))
  observations$observed_group_direction <- group_label[match(observations$normalized_exchange_id, group_ids)]
  list(source = source, crosswalk = crosswalk, observations = observations, source_path = source_path,
    crosswalk_path = crosswalk_path, group_ids = group_ids, group_label = stats::setNames(group_label, group_ids))
}

exercise6_model_features <- function(preflight, observation_data, settings) {
  fva <- preflight$fva
  mcf7_fva <- fva[fva$scenario == "MCF7_ACH_000019_Jain", , drop = FALSE]
  groups <- observation_data$group_ids
  features <- data.frame(exchange_reaction_id = groups, stringsAsFactors = FALSE)
  features <- merge(features, mcf7_fva[, c("reaction_id", "fba_flux", "fva_minimum", "fva_maximum", "fva_width", "fba_status", "fva_min_status", "fva_max_status", "medium_id", "objective_id", "model_sha256")],
    by.x = "exchange_reaction_id", by.y = "reaction_id", all.x = TRUE, sort = FALSE)
  features$zero_feasible <- ifelse(is.na(features$fva_minimum), NA, features$fva_minimum <= settings$model_zero_tolerance & features$fva_maximum >= -settings$model_zero_tolerance)
  features$can_uptake <- ifelse(is.na(features$fva_minimum), NA, features$fva_minimum < -settings$model_zero_tolerance)
  features$can_secrete <- ifelse(is.na(features$fva_maximum), NA, features$fva_maximum > settings$model_zero_tolerance)
  features$model_allowed_direction <- ifelse(is.na(features$fva_minimum), "unsupported",
    ifelse(!features$can_uptake & !features$can_secrete, "near_zero_or_ambiguous",
      ifelse(features$can_uptake & features$can_secrete, "uptake_or_secretion",
        ifelse(features$can_uptake, "uptake_only", "secretion_only"))))
  features$observed_group_direction <- unname(observation_data$group_label[features$exchange_reaction_id])
  features
}

exercise6_read_keibler <- function(root, settings) {
  source_path <- file.path(root, "datasets_and_models/treated/healthy/metabolomics/Keibler2021_Supporting_File_2.csv")
  source <- utils::read.csv(source_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("trial", "condition", "replicate", "metabolite_or_process", "observed_rate", "rate_unit", "exchange_reaction_id", "rate_in_cobra_exchange_sign", "mapping_status")
  exercise6_assert(all(required %in% names(source)), "Keibler source schema is incomplete.")
  source$source_row_id <- with(source, paste("KEIBLER", seq_len(nrow(source)), trial, condition, replicate, metabolite_or_process, sep = "::"))
  exercise6_assert(!anyDuplicated(source$source_row_id), "Keibler source surrogate row IDs are not unique.")
  eligible <- source$mapping_status %in% c("exact_name_or_id", "curated_alias") &
    !is.na(source$exchange_reaction_id) & nzchar(source$exchange_reaction_id) &
    !is.na(source$rate_in_cobra_exchange_sign) & is.finite(source$rate_in_cobra_exchange_sign)
  observations <- source[eligible, , drop = FALSE]
  observations$normalized_exchange_id <- exercise6_normalize_exchange(observations$exchange_reaction_id)
  observations$replicate_direction <- ifelse(observations$rate_in_cobra_exchange_sign < 0, "uptake",
    ifelse(observations$rate_in_cobra_exchange_sign > 0, "secretion", "near_zero_or_ambiguous"))
  group_key <- paste(observations$trial, observations$condition, observations$normalized_exchange_id, sep = "::")
  groups <- unique(group_key)
  labels <- vapply(groups, function(key) {
    signs <- observations$replicate_direction[group_key == key]
    if (any(signs == "near_zero_or_ambiguous") || length(unique(signs)) != 1L) "near_zero_or_ambiguous" else signs[[1L]]
  }, character(1))
  observations$measurement_group <- group_key
  observations$observed_group_direction <- labels[match(group_key, groups)]
  list(source = source, observations = observations, source_path = source_path,
    group_ids = unique(observations$normalized_exchange_id), groups = groups,
    group_label = stats::setNames(labels, groups))
}

exercise6_run_full_fva <- function(preflight, reaction_ids, scenario_name, fraction = 0.9, tolerance = 1e-7) {
  context <- if (scenario_name == "MCF7_ACH_000019_Jain") preflight$context else preflight$healthy
  ids <- sort(unique(reaction_ids))
  ids <- ids[ids %in% context$model$exchange_ids]
  exercise6_assert(length(ids) > 0L, paste("No eligible exchanges are present in", scenario_name))
  result <- run_fva(context$model, ids, fraction = fraction, tolerance = tolerance)
  exercise6_assert(result$status == "optimal", paste("Full FVA failed for", scenario_name, result$error %||% ""))
  intervals <- result$intervals
  intervals$scenario_id <- scenario_name
  intervals$medium_id <- context$model$medium_id
  intervals$objective_id <- context$model$objective_id
  intervals$objective_sense <- context$model$objective_sense
  intervals$model_sha256 <- context$model$source_sha256
  intervals$fva_fraction <- fraction
  intervals$model_unit <- "mmol/gDW/h"
  intervals
}

exercise6_fixed_category <- function(observed, feature, tolerance) {
  if (is.na(feature$fva_minimum) || feature$fba_status != "optimal" || feature$fva_min_status != "optimal" || feature$fva_max_status != "optimal") return("unsupported")
  if (observed == "near_zero_or_ambiguous") return("ambiguous")
  can_observed <- if (observed == "uptake") feature$fva_minimum < -tolerance else feature$fva_maximum > tolerance
  can_other <- if (observed == "uptake") feature$fva_maximum > tolerance else feature$fva_minimum < -tolerance
  can_zero <- feature$fva_minimum <= tolerance && feature$fva_maximum >= -tolerance
  if (!can_observed && !can_other && can_zero) return("blocked")
  if (!can_observed) return("contradicted")
  if (!can_other && !can_zero) return("required")
  "possible"
}

exercise6_fit_classifier <- function(training, feature_columns) {
  training <- training[training$observed_group_direction %in% c("uptake", "secretion") & complete.cases(training[, feature_columns, drop = FALSE]), , drop = FALSE]
  classes <- sort(unique(training$observed_group_direction))
  if (length(classes) < 2L) return(list(status = "unscored", reason = "single_class_training", n_training = nrow(training), state = NULL))
  center <- vapply(feature_columns, function(x) mean(training[[x]]), numeric(1))
  scale <- vapply(feature_columns, function(x) stats::sd(training[[x]]), numeric(1))
  scale[!is.finite(scale) | scale == 0] <- 1
  scaled <- sweep(sweep(as.matrix(training[, feature_columns, drop = FALSE]), 2L, center, "-"), 2L, scale, "/")
  centroids <- do.call(rbind, lapply(classes, function(label) colMeans(scaled[training$observed_group_direction == label, , drop = FALSE])))
  rownames(centroids) <- classes
  list(status = "fitted", reason = "", n_training = nrow(training), state = list(feature_columns = feature_columns, center = center, scale = scale, centroids = centroids))
}

exercise6_predict_classifier <- function(fit, row) {
  if (fit$status != "fitted") return(NA_character_)
  x <- as.numeric(row[fit$state$feature_columns])
  if (any(!is.finite(x))) return(NA_character_)
  x <- (x - fit$state$center) / fit$state$scale
  distances <- rowSums((fit$state$centroids - matrix(x, nrow = nrow(fit$state$centroids), ncol = length(x), byrow = TRUE))^2)
  rownames(fit$state$centroids)[which.min(distances)]
}

exercise6_metrics <- function(prediction, observed, method) {
  keep <- !is.na(prediction) & observed %in% c("uptake", "secretion")
  pred <- prediction[keep]; obs <- observed[keep]
  if (!length(obs)) return(data.frame(method = method, total = length(prediction), scoreable = 0L, accuracy = NA_real_, balanced_accuracy = NA_real_, macro_f1 = NA_real_, stringsAsFactors = FALSE))
  classes <- c("uptake", "secretion")
  cm <- table(factor(obs, levels = classes), factor(pred, levels = classes))
  recalls <- diag(cm) / rowSums(cm)
  recalls[!is.finite(recalls)] <- NA_real_
  precision <- diag(cm) / colSums(cm)
  precision[!is.finite(precision)] <- NA_real_
  f1 <- 2 * precision * recalls / (precision + recalls)
  f1[!is.finite(f1)] <- NA_real_
  data.frame(method = method, total = length(prediction), scoreable = length(obs), accuracy = mean(pred == obs),
    balanced_accuracy = mean(recalls, na.rm = TRUE), macro_f1 = mean(f1, na.rm = TRUE), stringsAsFactors = FALSE)
}

exercise6_make_confusion <- function(prediction, observed) {
  keep <- !is.na(prediction) & observed %in% c("uptake", "secretion")
  out <- as.data.frame.matrix(table(factor(observed[keep], levels = c("uptake", "secretion")), factor(prediction[keep], levels = c("uptake", "secretion"))))
  out$observed <- rownames(out)
  rownames(out) <- NULL
  out[, c("observed", "uptake", "secretion"), drop = FALSE]
}

exercise6_run_loocv <- function(features, feature_columns) {
  fold_rows <- list()
  for (i in seq_len(nrow(features))) {
    heldout <- features$exchange_reaction_id[[i]]
    train <- features[-i, , drop = FALSE]
    fit <- exercise6_fit_classifier(train, feature_columns)
    labels <- train$observed_group_direction[train$observed_group_direction %in% c("uptake", "secretion")]
    majority <- if (length(labels)) names(sort(table(labels), decreasing = TRUE))[1L] else NA_character_
    observed <- features$observed_group_direction[[i]]
    scored <- observed %in% c("uptake", "secretion") && fit$status == "fitted" && complete.cases(features[i, feature_columns, drop = FALSE])
    reason <- if (scored) "" else if (observed == "near_zero_or_ambiguous") "ambiguous_heldout_label" else fit$reason
    fold_rows[[i]] <- data.frame(fold_id = i, heldout_exchange_id = heldout,
      training_group_count = nrow(train), fitted_training_count = fit$n_training,
      heldout_observed_direction = observed,
      classifier_prediction = if (scored) exercise6_predict_classifier(fit, features[i, , drop = FALSE]) else NA_character_,
      majority_baseline_prediction = if (scored) majority else NA_character_, scoreable = scored,
      unscored_reason = reason,
      fitted_state_hash = if (fit$status == "fitted") digest::digest(fit$state, algo = "sha256") else NA_character_,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, fold_rows)
}

exercise6_plot_confusion <- function(confusion_tables, output_dir) {
  png(file.path(output_dir, "model_confusion_matrices.png"), width = 1400, height = 650, res = 140)
  par(mfrow = c(1, length(confusion_tables)), mar = c(5, 5, 4, 2) + 0.1)
  for (name in names(confusion_tables)) {
    cm <- confusion_tables[[name]]
    values <- as.matrix(cm[, c("uptake", "secretion"), drop = FALSE])
    image(x = 1:2, y = 1:2, z = t(values[2:1, , drop = FALSE]), axes = FALSE,
      col = colorRampPalette(c("#edf6f9", "#2a9d8f"))(20), main = name,
      xlab = "Predicted direction", ylab = "Observed direction")
    axis(1, at = 1:2, labels = c("uptake", "secretion"))
    axis(2, at = 1:2, labels = rev(c("uptake", "secretion")), las = 1)
    for (row in seq_len(nrow(values))) for (col in seq_len(ncol(values))) text(col, 3 - row, values[row, col], cex = 1.4, font = 2)
  }
  dev.off()
}

exercise6_run <- function(root = if (dir.exists("datasets_and_models")) "." else "..") {
  if (!exists("run_fva", mode = "function")) source(file.path(root, "R", "model_io.R"))
  output_dir <- file.path(root, "Exercise6_validation", "outputs")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  settings_path <- file.path(root, "config", "validation_settings.yml")
  settings <- exercise6_read_settings(settings_path); settings$path <- settings_path
  preflight <- exercise6_load_preflight(root, settings)
  observations <- exercise6_read_observations(root, settings)
  keibler <- exercise6_read_keibler(root, settings)
  full_mcf7_fva <- exercise6_run_full_fva(preflight, observations$group_ids, "MCF7_ACH_000019_Jain")
  full_gtex_fva <- exercise6_run_full_fva(preflight, keibler$group_ids, "GTEx_breast_Keibler")
  full_fva <- rbind(full_mcf7_fva, full_gtex_fva)
  full_mcf7_fva$fva_minimum <- full_mcf7_fva$minimum
  full_mcf7_fva$fva_maximum <- full_mcf7_fva$maximum
  full_mcf7_fva$fva_width <- full_mcf7_fva$maximum - full_mcf7_fva$minimum
  full_mcf7_fva$fba_status <- "optimal"
  full_mcf7_fva$fva_min_status <- full_mcf7_fva$min_status
  full_mcf7_fva$fva_max_status <- full_mcf7_fva$max_status
  full_mcf7_fva$scenario <- full_mcf7_fva$scenario_id
  preflight$fva <- full_mcf7_fva
  utils::write.csv(full_fva, file.path(output_dir, "all_available_exchange_fva.csv"), row.names = FALSE)
  utils::write.csv(keibler$observations, file.path(output_dir, "eligible_keibler_observation_qc.csv"), row.names = FALSE)
  features <- exercise6_model_features(preflight, observations, settings)
  feature_columns <- c("fba_flux", "fva_minimum", "fva_maximum", "fva_width", "zero_feasible", "can_uptake", "can_secrete")
  feature_columns <- feature_columns[feature_columns %in% names(features)]
  full_mcf7_features <- merge(observations$group_ids, full_mcf7_fva[, c("reaction_id", "minimum", "maximum", "fba_flux", "min_status", "max_status")], by.x = "x", by.y = "reaction_id", all.x = TRUE)
  names(full_mcf7_features)[names(full_mcf7_features) == "x"] <- "exchange_reaction_id"
  full_mcf7_features$observed_group_direction <- unname(observations$group_label[full_mcf7_features$exchange_reaction_id])
  utils::write.csv(full_mcf7_features, file.path(output_dir, "all_jain_exchange_assessment.csv"), row.names = FALSE)
  gtex_assessment <- merge(keibler$observations, full_gtex_fva[, c("reaction_id", "minimum", "maximum", "fba_flux", "min_status", "max_status")], by.x = "normalized_exchange_id", by.y = "reaction_id", all.x = TRUE)
  gtex_assessment$fva_minimum <- gtex_assessment$minimum
  gtex_assessment$fva_maximum <- gtex_assessment$maximum
  gtex_assessment$fba_status <- ifelse(gtex_assessment$min_status == "optimal" & gtex_assessment$max_status == "optimal", "optimal", "unsupported")
  gtex_assessment$fva_min_status <- gtex_assessment$min_status
  gtex_assessment$fva_max_status <- gtex_assessment$max_status
  gtex_assessment$fixed_model_category <- vapply(seq_len(nrow(gtex_assessment)), function(i) exercise6_fixed_category(gtex_assessment$observed_group_direction[[i]], gtex_assessment[i, ], settings$model_zero_tolerance), character(1))
  utils::write.csv(gtex_assessment, file.path(output_dir, "all_keibler_exchange_assessment.csv"), row.names = FALSE)

  gtex_labels <- vapply(keibler$group_ids, function(id) {
    values <- keibler$observations$replicate_direction[keibler$observations$normalized_exchange_id == id]
    if (any(values == "near_zero_or_ambiguous") || length(unique(values)) != 1L) "near_zero_or_ambiguous" else values[[1L]]
  }, character(1))
  gtex_features <- full_gtex_fva[, c("reaction_id", "fba_flux", "minimum", "maximum", "min_status", "max_status"), drop = FALSE]
  names(gtex_features) <- c("exchange_reaction_id", "fba_flux", "fva_minimum", "fva_maximum", "fva_min_status", "fva_max_status")
  gtex_features$fva_width <- gtex_features$fva_maximum - gtex_features$fva_minimum
  gtex_features$zero_feasible <- gtex_features$fva_minimum <= settings$model_zero_tolerance & gtex_features$fva_maximum >= -settings$model_zero_tolerance
  gtex_features$can_uptake <- gtex_features$fva_minimum < -settings$model_zero_tolerance
  gtex_features$can_secrete <- gtex_features$fva_maximum > settings$model_zero_tolerance
  gtex_features$observed_group_direction <- unname(stats::setNames(gtex_labels, keibler$group_ids)[gtex_features$exchange_reaction_id])
  gtex_folds <- exercise6_run_loocv(gtex_features, feature_columns)
  gtex_metrics <- rbind(exercise6_metrics(gtex_folds$classifier_prediction, gtex_folds$heldout_observed_direction, "GTEx/Keibler fold-trained GEM-feature direction classifier"), exercise6_metrics(gtex_folds$majority_baseline_prediction, gtex_folds$heldout_observed_direction, "GTEx/Keibler training-fold majority baseline"))
  utils::write.csv(gtex_features, file.path(output_dir, "gtex_keibler_classifier_features.csv"), row.names = FALSE)
  utils::write.csv(gtex_folds, file.path(output_dir, "gtex_keibler_loocv_predictions.csv"), row.names = FALSE)
  utils::write.csv(gtex_metrics, file.path(output_dir, "gtex_keibler_loocv_metrics.csv"), row.names = FALSE)
  source_qc <- observations$observations
  source_qc$source_row_id <- as.character(source_qc$source_row_id)
  utils::write.csv(source_qc, file.path(output_dir, "eligible_jain_observation_qc.csv"), row.names = FALSE)
  utils::write.csv(preflight$preflight, file.path(output_dir, "canonical_preflight_hashes.csv"), row.names = FALSE)

  classification <- features
  classification$fixed_model_category <- vapply(seq_len(nrow(classification)), function(i) exercise6_fixed_category(classification$observed_group_direction[[i]], classification[i, ], settings$model_zero_tolerance), character(1))
  classification$observation_unit <- settings$observation_unit
  classification$model_unit <- settings$model_unit
  classification$validation_settings_version <- settings$version
  utils::write.csv(classification, file.path(output_dir, "observation_vs_fva_classification.csv"), row.names = FALSE)
  coverage <- data.frame(category = c("eligible exchange groups", "supported FVA groups", "unsupported groups", "ambiguous observed groups", "fixed categories resolved"), count = c(nrow(features), sum(features$model_allowed_direction != "unsupported"), sum(features$model_allowed_direction == "unsupported"), sum(features$observed_group_direction == "near_zero_or_ambiguous"), sum(classification$fixed_model_category != "unsupported")), stringsAsFactors = FALSE)
  utils::write.csv(coverage, file.path(output_dir, "coverage_category_summary.csv"), row.names = FALSE)

  fold_rows <- list(); assignment_rows <- list()
  for (i in seq_len(nrow(features))) {
    heldout <- features$exchange_reaction_id[[i]]
    train <- features[-i, , drop = FALSE]
    heldout_rows <- source_qc$source_row_id[source_qc$normalized_exchange_id == heldout]
    train_rows <- source_qc$source_row_id[source_qc$normalized_exchange_id != heldout]
    leakage <- any(heldout_rows %in% train_rows) || any(train$exchange_reaction_id == heldout)
    exercise6_assert(!leakage, paste("Grouped fold leakage detected for", heldout))
    fit <- exercise6_fit_classifier(train, feature_columns)
    majority <- if (fit$n_training > 0L) {
      labels <- train$observed_group_direction[train$observed_group_direction %in% c("uptake", "secretion")]
      if (length(labels)) names(sort(table(labels), decreasing = TRUE))[1L] else NA_character_
    } else NA_character_
    heldout_observed <- features$observed_group_direction[[i]]
    scored <- heldout_observed %in% c("uptake", "secretion") && features$model_allowed_direction[[i]] != "unsupported" && fit$status == "fitted"
    reason <- if (scored) "" else if (heldout_observed == "near_zero_or_ambiguous") "ambiguous_heldout_label" else if (features$model_allowed_direction[[i]] == "unsupported") "unsupported_fva" else fit$reason
    prediction <- if (scored) exercise6_predict_classifier(fit, features[i, , drop = FALSE]) else NA_character_
    majority_prediction <- if (scored) majority else NA_character_
    fold_rows[[i]] <- data.frame(fold_id = i, heldout_exchange_id = heldout, training_group_count = nrow(train), fitted_training_count = fit$n_training, heldout_observed_direction = heldout_observed, classifier_prediction = prediction, majority_baseline_prediction = majority_prediction, scoreable = scored, unscored_reason = reason, fitted_state_hash = if (fit$status == "fitted") digest::digest(fit$state, algo = "sha256") else NA_character_, leakage_assertion = !leakage, stringsAsFactors = FALSE)
    assignment_rows[[i]] <- data.frame(fold_id = i, heldout_exchange_id = heldout, heldout_source_row_ids = paste(heldout_rows, collapse = ";"), training_exchange_ids = paste(train$exchange_reaction_id, collapse = ";"), training_source_row_count = length(train_rows), heldout_source_row_count = length(heldout_rows), no_heldout_source_row_in_training = !any(heldout_rows %in% train_rows), no_heldout_exchange_in_training = !any(train$exchange_reaction_id == heldout), stringsAsFactors = FALSE)
  }
  folds <- do.call(rbind, fold_rows); assignments <- do.call(rbind, assignment_rows)
  utils::write.csv(assignments, file.path(output_dir, "grouped_loocv_fold_assignments.csv"), row.names = FALSE)
  utils::write.csv(folds, file.path(output_dir, "grouped_loocv_predictions.csv"), row.names = FALSE)
  metrics <- rbind(exercise6_metrics(folds$classifier_prediction, folds$heldout_observed_direction, "fold-trained GEM-feature direction classifier"), exercise6_metrics(folds$majority_baseline_prediction, folds$heldout_observed_direction, "training-fold majority baseline"))
  confusion <- as.data.frame.matrix(table(factor(folds$heldout_observed_direction[folds$scoreable], levels = c("uptake", "secretion")), factor(folds$classifier_prediction[folds$scoreable], levels = c("uptake", "secretion"))))
  confusion$observed <- rownames(confusion); rownames(confusion) <- NULL
  utils::write.csv(confusion, file.path(output_dir, "loocv_confusion_matrix.csv"), row.names = FALSE)
  mcf7_confusion <- exercise6_make_confusion(folds$classifier_prediction, folds$heldout_observed_direction)
  gtex_confusion <- exercise6_make_confusion(gtex_folds$classifier_prediction, gtex_folds$heldout_observed_direction)
  utils::write.csv(mcf7_confusion, file.path(output_dir, "mcf7_jain_confusion_matrix.csv"), row.names = FALSE)
  utils::write.csv(gtex_confusion, file.path(output_dir, "gtex_keibler_confusion_matrix.csv"), row.names = FALSE)
  accuracy_summary <- rbind(metrics[1, ], gtex_metrics[1, ])
  names(accuracy_summary)[names(accuracy_summary) == "method"] <- "model_method"
  utils::write.csv(accuracy_summary, file.path(output_dir, "model_accuracy_summary.csv"), row.names = FALSE)
  exercise6_plot_confusion(list(MCF7_Jain = mcf7_confusion, GTEx_Keibler = gtex_confusion), output_dir)
  utils::write.csv(metrics, file.path(output_dir, "loocv_metrics.csv"), row.names = FALSE)
  settings_record <- data.frame(key = names(settings), value = vapply(settings, function(x) paste(x, collapse = ";"), character(1)), stringsAsFactors = FALSE)
  utils::write.csv(settings_record, file.path(output_dir, "validation_provenance_settings.csv"), row.names = FALSE)

  png(file.path(output_dir, "measured_direction_vs_model_allowed.png"), width = 1100, height = 700, res = 120)
  par(mar = c(11, 4, 4, 2) + 0.1)
  direction_table <- table(factor(classification$observed_group_direction, levels = c("uptake", "secretion", "near_zero_or_ambiguous")), factor(classification$model_allowed_direction, levels = c("uptake_only", "secretion_only", "uptake_or_secretion", "near_zero_or_ambiguous", "unsupported")))
  barplot(direction_table, beside = TRUE, legend.text = TRUE, las = 2, cex.names = 0.8, col = c("#2b6f9f", "#d47b3f", "#7f8c8d"), main = "Held-out Jain direction vs canonical model allowance", xlab = "Model-derived FVA direction category", ylab = "Exchange groups")
  dev.off()
  png(file.path(output_dir, "coverage_status.png"), width = 1000, height = 650, res = 120)
  par(mar = c(11, 4, 4, 2) + 0.1)
  barplot(coverage$count, names.arg = coverage$category, las = 2, cex.names = 0.8, col = c("#2b6f9f", "#4c956c", "#b56576", "#d47b3f", "#6c757d"), main = "Exercise 6 coverage and unresolved states", ylab = "Exchange groups")
  dev.off()
  png(file.path(output_dir, "loocv_predictions_observed.png"), width = 1200, height = 700, res = 120)
  par(mar = c(15, 5, 4, 2) + 0.1)
  plot(seq_len(nrow(folds)), rep(1, nrow(folds)), type = "n", xaxt = "n", yaxt = "n", xlab = "Held-out exchange fold", ylab = "Direction state", main = "Grouped LOOCV: every exchange group retained")
  axis(2, at = 1:3, labels = c("uptake", "secretion", "unscored/ambiguous"), las = 1)
  axis(1, at = seq_len(nrow(folds)), labels = folds$heldout_exchange_id, las = 2, cex.axis = 0.55)
  observed_y <- match(folds$heldout_observed_direction, c("uptake", "secretion", "near_zero_or_ambiguous")); observed_y[is.na(observed_y)] <- 3
  prediction_y <- match(folds$classifier_prediction, c("uptake", "secretion")); prediction_y[is.na(prediction_y)] <- 3
  points(seq_len(nrow(folds)) - 0.12, observed_y, pch = 16, col = "#2b6f9f")
  points(seq_len(nrow(folds)) + 0.12, prediction_y, pch = 17, col = "#d47b3f")
  legend("topright", legend = c("held-out observed direction", "classifier prediction/unscored"), pch = c(16, 17), col = c("#2b6f9f", "#d47b3f"), bty = "n")
  dev.off()
  list(settings = settings, preflight = preflight, observations = observations, features = features, classification = classification, folds = folds, metrics = metrics, coverage = coverage)
}