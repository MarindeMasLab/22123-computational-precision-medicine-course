# Run from hands-on_session/: Rscript R/run_full_screen.R
source("R/model_io.R")
source("R/ko_workflow.R")
out <- "Exercise4_disease_selective_knockouts/outputs"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
verified <- ko_context_preflight()
utils::write.csv(verified$qc, file.path(out, "native_wt_qc.csv"), row.names = FALSE)
rows <- ko_full_census(verified, file.path(out, "native_full_gene_ko.csv"))
paired <- ko_compare_census(rows)
symbols <- utils::read.csv("datasets_and_models/treated/cancer/transcriptomics/RECON1_DepMap_expression.csv",
  stringsAsFactors = FALSE, check.names = FALSE)[, c("recon1_gene_id", "gene_symbol")]
if (anyDuplicated(symbols$recon1_gene_id)) stop("Source symbol mapping is not unique.")
symbol_index <- match(sub("^G_", "", paired$gene_id), symbols$recon1_gene_id)
if (anyNA(symbol_index)) stop("Full census contains a gene missing from source symbol mapping.")
paired$gene_symbol <- symbols$gene_symbol[symbol_index]
utils::write.csv(paired, file.path(out, "native_full_gene_comparison.csv"), row.names = FALSE)
mapped <- ko_pathway_analysis(rows, paired)
utils::write.csv(mapped, file.path(out, "native_differential_ko_reaction_subsystems.csv"), row.names = FALSE)
if (nrow(mapped)) {
  groups <- split(mapped, list(mapped$context, mapped$subsystem), drop = TRUE)
  summary <- do.call(rbind, lapply(groups, function(x) data.frame(context = x$context[[1L]],
    subsystem = x$subsystem[[1L]], genes = length(unique(x$gene_id)),
    affected_reactions = length(unique(x$reaction_id)),
    gene_ids = paste(sort(unique(x$gene_id)), collapse = ";"), stringsAsFactors = FALSE)))
  rownames(summary) <- NULL
} else summary <- data.frame(context = character(), subsystem = character(), genes = integer(),
  affected_reactions = integer(), gene_ids = character())
utils::write.csv(summary, file.path(out, "native_differential_ko_pathways.csv"), row.names = FALSE)

png(file.path(out, "native_full_ko_ratios.png"), width = 1400, height = 850, res = 145)
valid <- paired$paired_eligible
plot(paired$healthy_ratio[valid], paired$cancer_ratio[valid], pch = 16, cex = 0.52,
  col = ifelse(paired$classification[valid] == "MCF7_lower_teaching_cutoffs", "firebrick",
    ifelse(paired$classification[valid] == "GTEx_lower_teaching_cutoffs", "steelblue", "grey45")),
  xlim = c(-0.04, max(1.05, paired$healthy_ratio[valid], na.rm = TRUE)),
  ylim = c(-0.04, max(1.05, paired$cancer_ratio[valid], na.rm = TRUE)),
  xlab = "GTEx/Keibler KO / own WT", ylab = "MCF7/RPMI KO / own WT",
  main = "Full Recon1 single-KO census (different native media; not selectivity validation)")
abline(0, 1, lty = 3); abline(h = c(1, 0.10), v = c(1, 0.50), lty = 2)
legend("bottomright", c("model ratios", "MCF7 lower (teaching)", "GTEx lower (teaching)"),
  col = c("grey45", "firebrick", "steelblue"), pch = 16, bty = "o", bg = "white", cex = 0.8)
mtext(sprintf("%d paired numeric / %d genes; %d excluded or failed (see complete CSV)",
  sum(valid), nrow(paired), sum(!valid)), side = 1, line = 4, cex = 0.8)
dev.off()

# One column per exact model gene product; unlike the ratio scatter this
# explicitly renders all screened IDs, including failures and no-GPR rows.
status_levels <- c("HiGHS_KO_LP", "WT_reused_identical_bounds", "not_tested", "failed")
status_colors <- c("#2d698c", "#64a573", "#b8b8b8", "#ce655b")
gene_order <- paired$gene_id
status_matrix <- sapply(names(verified$contexts), function(id) {
  x <- rows[rows$context == id, ]
  x <- x[match(gene_order, x$gene_id), ]
  match(x$method, status_levels)
})
stopifnot(!anyNA(status_matrix))
png(file.path(out, "native_full_ko_status_map.png"), width = 2400, height = 460, res = 150)
par(mar = c(4, 16, 3, 2))
image(seq_along(gene_order), seq_len(ncol(status_matrix)), status_matrix,
  col = status_colors, breaks = seq(0.5, length(status_levels) + 0.5, 1),
  axes = FALSE, xlab = "All Recon1 gene-product IDs in source-model order (see CSV for exact IDs)",
  ylab = "", main = "Complete native single-gene screen: execution status")
axis(2, at = seq_len(ncol(status_matrix)), labels = colnames(status_matrix), las = 1, cex.axis = 0.7)
legend("bottom", legend = status_levels, fill = status_colors, horiz = TRUE, bty = "n", inset = -0.1, xpd = TRUE, cex = 0.75)
dev.off()

if (nrow(summary)) {
  group_counts <- aggregate(genes ~ subsystem, data = summary, FUN = max)
  top <- head(group_counts[order(-group_counts$genes, group_counts$subsystem), ], 12L)
  png(file.path(out, "native_differential_ko_pathways.png"), width = 1500, height = 900, res = 145)
  par(mar = c(5, 18, 4, 2))
  labels <- ifelse(nchar(top$subsystem) > 50, paste0(substr(top$subsystem, 1, 47), "..."), top$subsystem)
  barplot(rev(top$genes), names.arg = rev(labels),
    horiz = TRUE, las = 1, col = "#5a85a5", cex.names = 0.68,
    xlab = "Distinct model gene IDs with changed bounds (not measured pathway activity)",
    main = "Differential-KO annotations (same GPR coverage in both contexts)")
  dev.off()
}
cat("Census:", nrow(rows), "rows; paired", sum(valid), "/", nrow(paired),
  "; categories:", paste(names(table(paired$classification)), table(paired$classification), collapse = ", "), "\n")