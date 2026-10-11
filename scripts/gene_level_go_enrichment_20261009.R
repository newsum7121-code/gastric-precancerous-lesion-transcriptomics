suppressPackageStartupMessages({library(AnnotationDbi); library(org.Hs.eg.db)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript gene_level_go_enrichment_20261009.R <gene_level_results.csv> <candidate_csv> <output_dir>")
res <- read.csv(args[[1]], check.names = FALSE, stringsAsFactors = FALSE)
cand <- read.csv(args[[2]], check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
out <- args[[3]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
univ <- unique(res$Gene[res$Gene_type == "protein_coding" & res$expressed_filter_pass])
if (all(c("Gene", "Gene_type") %in% names(cand))) {
  gene_set <- unique(cand$Gene[cand$Gene_type == "protein_coding"])
} else if (all(c("stable_gene_key", "evidence_tier") %in% names(cand))) {
  gene_set <- unique(cand$stable_gene_key[grepl("T1_REPLICATED_CORE", cand$evidence_tier, fixed = TRUE)])
} else stop("Candidate file must provide Gene/Gene_type or stable_gene_key/evidence_tier")
map <- AnnotationDbi::select(org.Hs.eg.db, keys = univ, keytype = "SYMBOL", columns = c("ENTREZID", "GO", "ONTOLOGY"))
map <- map[!is.na(map$GO) & map$ONTOLOGY == "BP" & !is.na(map$ENTREZID), ]
map <- unique(map[, c("SYMBOL", "ENTREZID", "GO")])
eligible <- unique(map$SYMBOL)
selected <- intersect(gene_set, eligible)
go_terms <- split(map$SYMBOL, map$GO)
go_terms <- lapply(go_terms, unique)
N <- length(eligible); n <- length(selected)
rows <- lapply(names(go_terms), function(go) {
  members <- intersect(go_terms[[go]], eligible)
  k <- length(intersect(members, selected)); M <- length(members)
  if (M < 5L || M > 500L) return(NULL)
  data.frame(GO = go, overlap = k, selected_n = n, term_size = M, universe_n = N,
    pvalue = phyper(k - 1L, M, N - M, n, lower.tail = FALSE),
    genes = paste(sort(intersect(members, selected)), collapse = ";"), stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
tab$BH_q <- p.adjust(tab$pvalue, method = "BH")
tab <- tab[order(tab$pvalue), ]
write.csv(tab, file.path(out, "gse163416_gene_level_GO_BP_hypergeom.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(data.frame(gene_set = c("strict_gene_level_candidates", "protein_coding_expressed_universe", "GO_eligible_universe"), n = c(length(selected), length(univ), N)), file.path(out, "gse163416_gene_level_GO_BP_universe_summary.csv"), row.names = FALSE)
writeLines(c(paste0("candidate symbols: ", paste(selected, collapse = ", ")), paste0("candidate mapped symbols n=", length(selected)), paste0("protein-coding expressed symbols n=", length(univ)), paste0("GO BP eligible universe symbols n=", N), paste0("tested BP terms (size 5-500): ", nrow(tab)), paste0("BH q<0.05 terms: ", sum(tab$BH_q < 0.05, na.rm = TRUE))), file.path(out, "gse163416_gene_level_GO_BP_summary.txt"), useBytes = TRUE)
print(readLines(file.path(out, "gse163416_gene_level_GO_BP_summary.txt")))
