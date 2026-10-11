suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript rerun_gse163416_gene_level_20261009.R <input_dir> <prior_candidate_csv> <output_dir>")
input_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
prior_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  early = file.path(input_dir, "GSE163416_N-CAG+IM.all.mRNA.txt.gz"),
  dysplasia = file.path(input_dir, "GSE163416_N-Dys.all.mRNA.txt.gz"),
  cancer = file.path(input_dir, "GSE163416_N-GC.all.mRNA.txt.gz")
)
stopifnot(all(file.exists(files)))
tables <- lapply(files, function(path) read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE))
required <- c("mRNA", "RNA_type", "Gene", "Gene_type", "N1", "N2", "N3", "log2(Fold_change)", "p-value", "q-value")
stopifnot(all(vapply(tables, function(x) all(required %in% names(x)), logical(1))))

keys <- lapply(tables, function(x) paste(x$mRNA, x$Gene, x$Gene_type, sep = "\034"))
stopifnot(all(vapply(keys, function(x) !anyDuplicated(x), logical(1))))
sample_columns <- list(
  early = c("N1", "N2", "N3", "CAG+IM1", "CAG+IM2", "CAG+IM3"),
  dysplasia = c("N1", "N2", "N3", "Dys1", "Dys2", "Dys3"),
  cancer = c("N1", "N2", "N3", "GC1", "GC2", "GC3")
)
for (nm in names(tables)) {
  x <- tables[[nm]]
  stopifnot(all(sample_columns[[nm]] %in% names(x)))
  stopifnot(all(is.finite(as.matrix(x[, sample_columns[[nm]], drop = FALSE]))))
  stopifnot(all(as.matrix(x[, sample_columns[[nm]], drop = FALSE]) >= 0))
  stopifnot(all(as.matrix(x[, sample_columns[[nm]], drop = FALSE]) %% 1 == 0))
}

common_keys <- Reduce(intersect, keys)
for (nm in names(tables)[-1]) {
  i0 <- match(common_keys, keys[[1]])
  ik <- match(common_keys, keys[[nm]])
  stopifnot(all(as.matrix(tables[[1]][i0, c("N1", "N2", "N3")]) == as.matrix(tables[[nm]][ik, c("N1", "N2", "N3")])) )
}
all_keys <- sort(common_keys)
meta_index <- match(all_keys, keys[[1]])
meta <- tables[[1]][meta_index, c("mRNA", "RNA_type", "Gene", "Gene_type"), drop = FALSE]
missing_meta <- is.na(meta$mRNA)
if (any(missing_meta)) {
  for (nm in names(tables)[-1]) {
    j <- match(all_keys[missing_meta], keys[[nm]])
    use <- !is.na(j)
    if (any(use)) {
      at <- which(missing_meta)[use]
      meta[at, c("mRNA", "RNA_type", "Gene", "Gene_type")] <- tables[[nm]][j[use], c("mRNA", "RNA_type", "Gene", "Gene_type"), drop = FALSE]
      missing_meta[at] <- FALSE
    }
  }
}
stopifnot(!any(missing_meta), !anyNA(meta$Gene), all(nzchar(meta$Gene)))
stopifnot(!anyNA(meta$Gene), all(nzchar(meta$Gene)))
group_key <- paste(meta$Gene, meta$Gene_type, sep = "\034")
gene_meta <- meta[!duplicated(group_key), , drop = FALSE]
gene_meta$Gene_key <- group_key[!duplicated(group_key)]
idx <- match(group_key, group_key[!duplicated(group_key)])
sample_map <- c("N1", "N2", "N3", "CAG+IM1", "CAG+IM2", "CAG+IM3", "Dys1", "Dys2", "Dys3", "GC1", "GC2", "GC3")
sample_owner <- c(rep("early", 6), rep("dysplasia", 3), rep("cancer", 3))
tx_count <- vapply(seq_along(sample_map), function(i) {
  v <- numeric(length(all_keys))
  j <- match(all_keys, keys[[sample_owner[i]]])
  present <- !is.na(j)
  v[present] <- tables[[sample_owner[i]]][[sample_map[i]]][j[present]]
  v
}, numeric(length(all_keys)))
tx_count <- as.data.frame(tx_count)
colnames(tx_count) <- c("CSG1", "CSG2", "CSG3", "CAGIM1", "CAGIM2", "CAGIM3", "Dys1", "Dys2", "Dys3", "GC1", "GC2", "GC3")
gene_counts <- rowsum(as.matrix(tx_count), group = idx, reorder = FALSE)
rownames(gene_counts) <- gene_meta$Gene_key
storage.mode(gene_counts) <- "integer"

stage <- factor(rep(c("CSG", "CAGIM", "Dys", "GC"), each = 3), levels = c("CSG", "CAGIM", "Dys", "GC"))
design <- model.matrix(~ 0 + stage)
colnames(design) <- levels(stage)
y_all <- DGEList(counts = gene_counts, genes = gene_meta[, c("Gene", "Gene_type")])
keep <- filterByExpr(y_all, design = design)
y <- y_all[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
contrasts <- makeContrasts(
  CAGIM_vs_CSG = CAGIM - CSG,
  Dys_vs_CSG = Dys - CSG,
  GC_vs_CSG = GC - CSG,
  levels = design
)

tables_out <- lapply(seq_len(ncol(contrasts)), function(i) {
  test <- glmQLFTest(fit, contrast = contrasts[, i])
  tab <- topTags(test, n = Inf, sort.by = "none")$table
  tab$Gene_key <- rownames(tab)
  names(tab)[names(tab) == "logFC"] <- paste0("logFC_", colnames(contrasts)[i])
  names(tab)[names(tab) == "PValue"] <- paste0("p_", colnames(contrasts)[i])
  names(tab)[names(tab) == "FDR"] <- paste0("q_", colnames(contrasts)[i])
  tab[, c("Gene_key", grep("^(logFC|p_|q_)", names(tab), value = TRUE)), drop = FALSE]
})
names(tables_out) <- colnames(contrasts)

logcpm <- cpm(y, log = TRUE, prior.count = 2)
mean_by_stage <- sapply(levels(stage), function(s) rowMeans(logcpm[, stage == s, drop = FALSE]))
rho12 <- apply(logcpm, 1, function(v) suppressWarnings(cor(as.numeric(stage), v, method = "spearman")))
monotonic_up <- apply(mean_by_stage, 1, function(v) all(diff(v) >= 0) && any(diff(v) > 0))
monotonic_down <- apply(mean_by_stage, 1, function(v) all(diff(v) <= 0) && any(diff(v) < 0))
monotonic <- ifelse(monotonic_up, "monotonic_up", ifelse(monotonic_down, "monotonic_down", "non_monotonic"))
sample_logcpm_out <- data.frame(Gene_key = rownames(logcpm), logcpm, check.names = FALSE)
sample_logcpm_out$Gene <- gene_meta$Gene[match(sample_logcpm_out$Gene_key, gene_meta$Gene_key)]
sample_logcpm_out$Gene_type <- gene_meta$Gene_type[match(sample_logcpm_out$Gene_key, gene_meta$Gene_key)]
write.csv(sample_logcpm_out, file.path(out_dir, "gse163416_gene_level_sample_logCPM.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")

all_results <- Reduce(function(a, b) merge(a, b, by = "Gene_key", all = TRUE, sort = FALSE), tables_out)
all_results$Gene <- gene_meta$Gene[match(all_results$Gene_key, gene_meta$Gene_key)]
all_results$Gene_type <- gene_meta$Gene_type[match(all_results$Gene_key, gene_meta$Gene_key)]
all_results$monotonic_pattern <- monotonic[match(all_results$Gene_key, rownames(logcpm))]
all_results$spearman_rho_12_samples <- rho12[match(all_results$Gene_key, rownames(logcpm))]
for (j in seq_along(levels(stage))) all_results[[paste0("mean_logCPM_", levels(stage)[j])]] <- mean_by_stage[match(all_results$Gene_key, rownames(mean_by_stage)), j]
all_results$expressed_filter_pass <- TRUE

prior_tbl_for_loo <- read.csv(prior_path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
if (all(c("stable_gene_key", "evidence_tier") %in% names(prior_tbl_for_loo))) {
  prior_core_symbols <- unique(prior_tbl_for_loo$stable_gene_key[grepl("T1_REPLICATED_CORE", prior_tbl_for_loo$evidence_tier, fixed = TRUE)])
} else if (all(c("gene", "included_in_122_gene_core") %in% names(prior_tbl_for_loo))) {
  included <- as.logical(prior_tbl_for_loo$included_in_122_gene_core)
  stopifnot(!anyNA(included))
  prior_core_symbols <- unique(prior_tbl_for_loo$gene[included])
} else {
  stop("Prior candidate CSV must contain either stable_gene_key/evidence_tier or gene/included_in_122_gene_core columns.")
}
stopifnot(length(prior_core_symbols) > 0L, all(nzchar(prior_core_symbols)))
loo_gene_keys <- rownames(logcpm)[gene_meta$Gene[match(rownames(logcpm), gene_meta$Gene_key)] %in% prior_core_symbols & gene_meta$Gene_type[match(rownames(logcpm), gene_meta$Gene_key)] == "protein_coding"]
loo_rows <- lapply(loo_gene_keys, function(gk) {
  full_direction <- monotonic[match(gk, rownames(logcpm))]
  vapply(seq_along(stage), function(drop_i) {
    keep_i <- setdiff(seq_along(stage), drop_i)
    sub_means <- vapply(levels(stage), function(s) mean(logcpm[gk, keep_i[stage[keep_i] == s]]), numeric(1))
    sub_up <- all(diff(sub_means) >= 0) && any(diff(sub_means) > 0)
    sub_down <- all(diff(sub_means) <= 0) && any(diff(sub_means) < 0)
    sub_pattern <- if (sub_up) "monotonic_up" else if (sub_down) "monotonic_down" else "non_monotonic"
    sub_rho <- suppressWarnings(cor(as.numeric(stage[keep_i]), as.numeric(logcpm[gk, keep_i]), method = "spearman"))
    as.numeric(sub_pattern == full_direction && full_direction != "non_monotonic" && abs(sub_rho) >= 0.8)
  }, numeric(1))
})
loo_mat <- do.call(rbind, loo_rows)
rownames(loo_mat) <- loo_gene_keys
colnames(loo_mat) <- paste0("drop_", colnames(logcpm))
loo_out <- data.frame(Gene_key = loo_gene_keys, loo_rule_pass_n = rowSums(loo_mat), loo_rule_pass_fraction = rowMeans(loo_mat), loo_mat, check.names = FALSE)
loo_out$Gene <- gene_meta$Gene[match(loo_out$Gene_key, gene_meta$Gene_key)]
loo_out$Gene_type <- gene_meta$Gene_type[match(loo_out$Gene_key, gene_meta$Gene_key)]
loo_out$full_data_pattern <- monotonic[match(loo_out$Gene_key, rownames(logcpm))]
loo_out$spearman_rho_12_samples <- rho12[match(loo_out$Gene_key, rownames(logcpm))]
loo_out$previously_reported_core_122 <- loo_out$Gene %in% prior_core_symbols
write.csv(loo_out, file.path(out_dir, "gse163416_leave_one_sample_out_trend_stability.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")

all_results$discovery_candidate <- all_results$Gene_type == "protein_coding" &
  abs(all_results$logFC_CAGIM_vs_CSG) >= 1 & all_results$q_CAGIM_vs_CSG < 0.05 &
  abs(all_results$logFC_GC_vs_CSG) >= 1 & all_results$q_GC_vs_CSG < 0.05 &
  all_results$monotonic_pattern != "non_monotonic" & abs(all_results$spearman_rho_12_samples) >= 0.8

prior <- read.csv(prior_path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
stopifnot("stable_gene_key" %in% names(prior), "evidence_tier" %in% names(prior))
prior$previously_reported_130_union <- TRUE
prior$previously_reported_core_122 <- grepl("T1_REPLICATED_CORE", prior$evidence_tier, fixed = TRUE)
prior$previously_reported_early_only_8 <- prior$evidence_tier == "T2_EARLY_NOMINAL_ONLY"
all_results$previously_reported_130_union <- all_results$Gene %in% prior$stable_gene_key
all_results$previously_reported_core_122 <- all_results$Gene %in% prior$stable_gene_key[prior$previously_reported_core_122]
all_results$previously_reported_early_only_8 <- all_results$Gene %in% prior$stable_gene_key[prior$previously_reported_early_only_8]

prior_core <- unique(prior$stable_gene_key[prior$previously_reported_core_122])
new <- unique(all_results$Gene[all_results$discovery_candidate])
intersect_core <- intersect(new, prior_core)
old_tier <- setNames(prior$evidence_tier, prior$stable_gene_key)
all_results$prior_evidence_tier <- unname(old_tier[all_results$Gene])
all_results <- all_results[order(!all_results$discovery_candidate, all_results$q_GC_vs_CSG, all_results$q_CAGIM_vs_CSG), ]

write.csv(all_results, file.path(out_dir, "gse163416_gene_level_edger_results.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(all_results[all_results$discovery_candidate, ], file.path(out_dir, "gse163416_gene_level_discovery_candidates.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(data.frame(
  previous_core_gene = prior_core,
  present_in_new_discovery_set = prior_core %in% new,
  stringsAsFactors = FALSE
), file.path(out_dir, "previous_core_recovery.csv"), row.names = FALSE, fileEncoding = "UTF-8")

summary <- list(
  method = "Transcript rows with supplied integer sample counts were summed by Gene and Gene_type; filterByExpr; TMM normalization; edgeR robust dispersion and quasi-likelihood F tests; BH FDR across the expressed-gene universe separately per contrast. Candidate Spearman rho is computed across all 12 sample-level logCPM values, not four stage means.",
  dataset = "GSE163416",
  input_rows_transcript = nrow(tables[[1]]),
  transcript_ids_unique = length(unique(meta$mRNA)),
  feature_alignment = "Only transcript/Gene/Gene_type keys shared by all three contrast tables were retained; no absent table row was imputed as a zero count.",
  gene_gene_type_groups = nrow(gene_meta),
  expressed_gene_groups_tested = nrow(y),
  contrasts = list(
    CAGIM_vs_CSG = sum(all_results$q_CAGIM_vs_CSG < 0.05, na.rm = TRUE),
    Dys_vs_CSG = sum(all_results$q_Dys_vs_CSG < 0.05, na.rm = TRUE),
    GC_vs_CSG = sum(all_results$q_GC_vs_CSG < 0.05, na.rm = TRUE)
  ),
  prior_core_count = length(prior_core),
  refit_discovery_candidate_count = length(new),
  prior_core_recovered_count = length(intersect_core),
  prior_core_not_recovered_count = length(setdiff(prior_core, new)),
  prior_core_recovery_fraction = length(intersect_core) / length(prior_core),
  new_candidates_prior_core_overlap = sum(all_results$discovery_candidate & all_results$previously_reported_core_122, na.rm = TRUE),
  new_candidates_early_only_overlap = sum(all_results$discovery_candidate & all_results$previously_reported_early_only_8, na.rm = TRUE),
  new_candidates_not_in_prior_union = sum(all_results$discovery_candidate & !all_results$previously_reported_130_union, na.rm = TRUE),
  prior_core_early_contrast_fdr05 = sum(all_results$previously_reported_core_122 & all_results$q_CAGIM_vs_CSG < 0.05, na.rm = TRUE),
  prior_core_cancer_contrast_fdr05 = sum(all_results$previously_reported_core_122 & all_results$q_GC_vs_CSG < 0.05, na.rm = TRUE),
  candidate_rule = "protein-coding; |edgeR log2FC|>=1 and per-contrast BH FDR<0.05 for CAG+IM vs CSG and GC vs CSG; monotonic nonzero stage means in TMM logCPM; |Spearman rho|>=0.8 across 12 samples.",
  source_counts_are_raw_integer_values = TRUE,
  raw_fastq_reprocessed = FALSE,
  software = list(R = R.version.string, edgeR = as.character(packageVersion("edgeR"))),
  session = capture.output(sessionInfo())
)
jsonlite::write_json(summary, file.path(out_dir, "gse163416_gene_level_edger_summary.json"), pretty = TRUE, auto_unbox = TRUE, dataframe = "rows")
write.csv(data.frame(file = basename(files), transcript_rows = vapply(tables, nrow, integer(1)), distinct_symbols = vapply(tables, function(x) length(unique(x$Gene)), integer(1))), file.path(out_dir, "input_reconciliation.csv"), row.names = FALSE)
print(jsonlite::toJSON(summary[setdiff(names(summary), "session")], pretty = TRUE, auto_unbox = TRUE))
