suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: Rscript gene_level_candidate_leave_one_sample_out_20261010.R <input_dir> <primary_results.csv> <candidate_csv> <output_dir>")
input_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
results_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
candidates_path <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
out_dir <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  early = file.path(input_dir, "GSE163416_N-CAG+IM.all.mRNA.txt.gz"),
  dysplasia = file.path(input_dir, "GSE163416_N-Dys.all.mRNA.txt.gz"),
  cancer = file.path(input_dir, "GSE163416_N-GC.all.mRNA.txt.gz")
)
stopifnot(all(file.exists(files)))
tables <- lapply(files, function(path) read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE))
keys <- lapply(tables, function(x) paste(x$mRNA, x$Gene, x$Gene_type, sep = "\034"))
common_keys <- sort(Reduce(intersect, keys))

sample_map <- c("N1", "N2", "N3", "CAG+IM1", "CAG+IM2", "CAG+IM3", "Dys1", "Dys2", "Dys3", "GC1", "GC2", "GC3")
sample_owner <- c(rep("early", 6), rep("dysplasia", 3), rep("cancer", 3))
sample_names <- c("CSG1", "CSG2", "CSG3", "CAGIM1", "CAGIM2", "CAGIM3", "Dys1", "Dys2", "Dys3", "GC1", "GC2", "GC3")
stage_levels <- c("CSG", "CAGIM", "Dys", "GC")
stages <- factor(rep(stage_levels, each = 3), levels = stage_levels)

meta <- tables[[1]][match(common_keys, keys[[1]]), c("mRNA", "Gene", "Gene_type"), drop = FALSE]
gene_key_tx <- paste(meta$Gene, meta$Gene_type, sep = "\034")
unique_gene <- !duplicated(gene_key_tx)
gene_meta <- meta[unique_gene, c("Gene", "Gene_type"), drop = FALSE]
gene_meta$Gene_key <- gene_key_tx[unique_gene]
gene_index <- match(gene_key_tx, gene_meta$Gene_key)

tx_counts <- vapply(seq_along(sample_map), function(i) {
  j <- match(common_keys, keys[[sample_owner[i]]])
  stopifnot(!anyNA(j))
  as.numeric(tables[[sample_owner[i]]][[sample_map[i]]][j])
}, numeric(length(common_keys)))
colnames(tx_counts) <- sample_names
gene_counts <- rowsum(tx_counts, group = gene_index, reorder = FALSE)
rownames(gene_counts) <- gene_meta$Gene_key
storage.mode(gene_counts) <- "integer"

primary <- read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
candidates <- read.csv(candidates_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("Gene", "Gene_type", "discovery_candidate") %in% names(primary)))
stopifnot(all(c("Gene", "Gene_type") %in% names(candidates)))
candidate_keys <- unique(paste(candidates$Gene, candidates$Gene_type, sep = "\034"))
primary_keys <- unique(paste(primary$Gene[primary$discovery_candidate], primary$Gene_type[primary$discovery_candidate], sep = "\034"))
stopifnot(setequal(candidate_keys, primary_keys))

fit_one <- function(keep_idx, omit_label) {
  group <- droplevels(stages[keep_idx])
  design <- model.matrix(~0 + group)
  colnames(design) <- sub("^group", "", colnames(design))
  y0 <- DGEList(counts = gene_counts[, keep_idx, drop = FALSE], genes = gene_meta[, c("Gene", "Gene_type"), drop = FALSE])
  keep_gene <- filterByExpr(y0, design = design)
  y <- calcNormFactors(y0[keep_gene, , keep.lib.sizes = FALSE], method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  contrast_early <- makeContrasts(CAGIM - CSG, levels = design)
  contrast_cancer <- makeContrasts(GC - CSG, levels = design)
  early <- topTags(glmQLFTest(fit, contrast = contrast_early), n = Inf, sort.by = "none")$table
  cancer <- topTags(glmQLFTest(fit, contrast = contrast_cancer), n = Inf, sort.by = "none")$table
  eligible <- rownames(y)
  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  rows <- lapply(candidate_keys, function(gk) {
    g_parts <- strsplit(gk, "\034", fixed = TRUE)[[1]]
    gene <- g_parts[[1]]
    gene_type <- g_parts[[2]]
    if (!gk %in% eligible) {
      return(data.frame(omitted_sample = omit_label, Gene = gene, Gene_type = gene_type,
        expressed_filter_pass = FALSE, logFC_CAGIM_vs_CSG = NA_real_, q_CAGIM_vs_CSG = NA_real_,
        logFC_GC_vs_CSG = NA_real_, q_GC_vs_CSG = NA_real_, monotonic_pattern = "not_tested",
        spearman_rho = NA_real_, effect_FDR_pass = FALSE, trend_pass = FALSE, full_rule_pass = FALSE))
    }
    means <- sapply(levels(group), function(level) mean(logcpm[gk, group == level, drop = TRUE]))
    up <- all(diff(means) >= 0) && any(diff(means) > 0)
    down <- all(diff(means) <= 0) && any(diff(means) < 0)
    pattern <- if (up) "monotonic_up" else if (down) "monotonic_down" else "non_monotonic"
    rho <- suppressWarnings(cor(as.numeric(group), as.numeric(logcpm[gk, ]), method = "spearman"))
    e <- early[gk, , drop = FALSE]
    c <- cancer[gk, , drop = FALSE]
    effect_fdr <- abs(e$logFC[[1]]) >= 1 && e$FDR[[1]] < 0.05 && abs(c$logFC[[1]]) >= 1 && c$FDR[[1]] < 0.05
    trend <- pattern != "non_monotonic" && abs(rho) >= 0.8
    data.frame(omitted_sample = omit_label, Gene = gene, Gene_type = gene_type,
      expressed_filter_pass = TRUE, logFC_CAGIM_vs_CSG = e$logFC[[1]], q_CAGIM_vs_CSG = e$FDR[[1]],
      logFC_GC_vs_CSG = c$logFC[[1]], q_GC_vs_CSG = c$FDR[[1]], monotonic_pattern = pattern,
      spearman_rho = rho, effect_FDR_pass = effect_fdr, trend_pass = trend, full_rule_pass = effect_fdr && trend)
  })
  do.call(rbind, rows)
}

out <- do.call(rbind, lapply(seq_along(sample_names), function(i) fit_one(setdiff(seq_along(sample_names), i), sample_names[i])))
write.csv(out, file.path(out_dir, "gse163416_five_candidate_leave_one_sample_out.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
summary <- do.call(rbind, lapply(split(out, out$Gene), function(x) data.frame(
  Gene = x$Gene[1], omitted_samples = nrow(x), expressed_filter_pass_n = sum(x$expressed_filter_pass),
  effect_FDR_pass_n = sum(x$effect_FDR_pass), trend_pass_n = sum(x$trend_pass), full_rule_pass_n = sum(x$full_rule_pass),
  full_rule_pass_fraction = mean(x$full_rule_pass), q_CAGIM_min = min(x$q_CAGIM_vs_CSG, na.rm = TRUE),
  q_CAGIM_max = max(x$q_CAGIM_vs_CSG, na.rm = TRUE), q_GC_min = min(x$q_GC_vs_CSG, na.rm = TRUE),
  q_GC_max = max(x$q_GC_vs_CSG, na.rm = TRUE), stringsAsFactors = FALSE)))
write.csv(summary, file.path(out_dir, "gse163416_five_candidate_leave_one_sample_out_summary.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
writeLines(c(
  "Leave-one-sample-out refits use the same common transcript/Gene/Gene_type features, Gene+Gene_type aggregation, filterByExpr, TMM, robust edgeR quasi-likelihood model, two primary contrasts, and per-contrast BH correction as the full model.",
  "For each omission, the five primary full-model candidates are held fixed; the complete candidate rule is re-evaluated in the remaining 11 samples.",
  "These are influence diagnostics in a small discovery cohort, not independent validation or a substitute for prospective replication."
), file.path(out_dir, "README.txt"), useBytes = TRUE)
print(summary)
