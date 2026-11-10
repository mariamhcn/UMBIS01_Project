# Helpers for Q3: functional enrichment (ORA then GSEA) of the Q2 pseudobulk DE results,
# with clusterProfiler.
#
# Gene sets are handled in two ways:
#  - GO Biological Process: clusterProfiler's own functions (enrichGO / gseGO) with the org.Hs.eg.db
#    annotation, as in the Datanovia tutorial; this also allows simplify() to remove redundant GO terms;
#  - other collections (Hallmark, Reactome, SynGO): generic functions (enricher / GSEA) with a
#    TERM2GENE table (columns gs_name, gene_symbol).

DE_REQUIRED_COLUMNS <- c("cell_type", "contrast", "gene", "log2FoldChange", "stat", "pvalue", "padj")

# Read the Q2 table and check it follows the shared contract (see README.md at the project root).
load_de_results <- function(path) {
  de <- utils::read.csv(path, stringsAsFactors = FALSE)
  missing <- setdiff(DE_REQUIRED_COLUMNS, colnames(de))
  if (length(missing) > 0) {
    stop("Q2 table is missing columns: ", paste(missing, collapse = ", "))
  }
  de
}

# ---------------------------------------------------------------------------
# ORA
# ---------------------------------------------------------------------------

# Lists of significant genes for one cell type: up and down separately, for each disease.
# Up and down are separated because an ORA ignores direction: mixing them would hide
# whether a function is induced or repressed.
make_deg_lists <- function(de_ct, alpha = 0.05, contrasts = c(SCZ = "SCZ_vs_CTRL", BD = "BD_vs_CTRL")) {
  lists <- list()
  for (dx in names(contrasts)) {
    d <- de_ct[de_ct$contrast == contrasts[[dx]] & !is.na(de_ct$padj), ]
    lists[[paste0(dx, "_up")]]   <- unique(d$gene[d$padj < alpha & d$log2FoldChange > 0])
    lists[[paste0(dx, "_down")]] <- unique(d$gene[d$padj < alpha & d$log2FoldChange < 0])
  }
  lists
}

# ORA of several gene lists at once with compareCluster (one column per list in the dotplot).
# term2gene = NULL means GO Biological Process through enrichGO.
run_ora_compare <- function(gene_lists, universe, term2gene = NULL,
                            alpha = 0.05, min_size = 15, max_size = 500, min_genes = 5) {
  gene_lists <- gene_lists[vapply(gene_lists, length, integer(1)) >= min_genes]
  if (length(gene_lists) == 0) return(NULL)
  res <- tryCatch({
    if (is.null(term2gene)) {
      clusterProfiler::compareCluster(
        geneClusters  = gene_lists, fun = "enrichGO",
        OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
        universe = universe, pAdjustMethod = "BH",
        pvalueCutoff = alpha, qvalueCutoff = 1,
        minGSSize = min_size, maxGSSize = max_size, readable = FALSE)
    } else {
      clusterProfiler::compareCluster(
        geneClusters  = gene_lists, fun = "enricher",
        TERM2GENE = term2gene, universe = universe, pAdjustMethod = "BH",
        pvalueCutoff = alpha, qvalueCutoff = 1,
        minGSSize = min_size, maxGSSize = max_size)
    }
  }, error = function(e) NULL)   # compareCluster stops when no list has any enriched term
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  res
}

# ORA of one gene list (used for the cnetplot example).
run_ora_one <- function(genes, universe, term2gene = NULL, alpha = 0.05, min_size = 15, max_size = 500) {
  if (is.null(term2gene)) {
    clusterProfiler::enrichGO(genes, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                              universe = universe, pAdjustMethod = "BH", pvalueCutoff = alpha,
                              qvalueCutoff = 1, minGSSize = min_size, maxGSSize = max_size)
  } else {
    clusterProfiler::enricher(genes, TERM2GENE = term2gene, universe = universe, pAdjustMethod = "BH",
                              pvalueCutoff = alpha, qvalueCutoff = 1, minGSSize = min_size, maxGSSize = max_size)
  }
}

# Remove redundant GO terms (keeps the most significant term of each group of similar terms).
simplify_go <- function(x, cutoff = 0.7) {
  clusterProfiler::simplify(x, cutoff = cutoff, by = "p.adjust", select_fun = min)
}

# ---------------------------------------------------------------------------
# GSEA
# ---------------------------------------------------------------------------

# Ranked gene list for GSEA: one value per gene, sorted decreasingly.
# The DESeq2 Wald statistic combines effect size and its uncertainty, and keeps the sign.
make_ranked_list <- function(de_one, metric = "stat") {
  de_one <- de_one[!is.na(de_one[[metric]]) & !is.na(de_one$gene), ]
  # if a symbol appears twice, keep the strongest value
  de_one <- de_one[order(-abs(de_one[[metric]])), ]
  de_one <- de_one[!duplicated(de_one$gene), ]
  ranks <- stats::setNames(de_one[[metric]], de_one$gene)
  sort(ranks, decreasing = TRUE)
}

# GSEA for one ranked list, returning the clusterProfiler object (needed for gseaplot2 and dotplot).
# term2gene = NULL means GO Biological Process through gseGO.
# pvalueCutoff = 1 keeps every tested set, so SCZ and BD can be compared on the same pathways afterwards.
run_gsea_object <- function(ranks, term2gene = NULL, min_size = 15, max_size = 500, seed = 42) {
  set.seed(seed)
  if (is.null(term2gene)) {
    clusterProfiler::gseGO(geneList = ranks, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL",
                           ont = "BP", minGSSize = min_size, maxGSSize = max_size,
                           pvalueCutoff = 1, eps = 1e-10, seed = TRUE, verbose = FALSE)
  } else {
    clusterProfiler::GSEA(geneList = ranks, TERM2GENE = term2gene,
                          minGSSize = min_size, maxGSSize = max_size,
                          pvalueCutoff = 1, eps = 1e-10, seed = TRUE, verbose = FALSE)
  }
}

# Run GSEA for every cell type x contrast of one collection and stack the results in a data frame.
run_gsea_all <- function(de, term2gene = NULL, metric = "stat", ...) {
  groups <- unique(de[, c("cell_type", "contrast")])
  out <- lapply(seq_len(nrow(groups)), function(i) {
    ct <- groups$cell_type[i]
    co <- groups$contrast[i]
    ranks <- make_ranked_list(de[de$cell_type == ct & de$contrast == co, ], metric)
    res <- as.data.frame(run_gsea_object(ranks, term2gene, ...))
    if (nrow(res) == 0) return(NULL)
    cbind(cell_type = ct, contrast = co, res)
  })
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# SCZ vs BD comparison
# ---------------------------------------------------------------------------

# Put SCZ and BD side by side for each cell type x pathway.
pair_contrasts <- function(gsea_all,
                           contrast_a = "SCZ_vs_CTRL",
                           contrast_b = "BD_vs_CTRL") {
  keep <- c("cell_type", "ID", "Description", "NES", "p.adjust", "core_enrichment")
  a <- gsea_all[gsea_all$contrast == contrast_a, keep]
  b <- gsea_all[gsea_all$contrast == contrast_b, keep]
  merge(a, b, by = c("cell_type", "ID", "Description"), suffixes = c("_SCZ", "_BD"))
}

# Label each pathway: shared (same direction), opposite, specific to one disease, or not significant.
classify_pathways <- function(paired, alpha = 0.05) {
  sig_scz <- paired$p.adjust_SCZ < alpha
  sig_bd  <- paired$p.adjust_BD  < alpha
  same    <- sign(paired$NES_SCZ) == sign(paired$NES_BD)
  paired$category <- ifelse(sig_scz & sig_bd & same,  "shared_concordant",
                     ifelse(sig_scz & sig_bd & !same, "shared_opposite",
                     ifelse(sig_scz & !sig_bd,         "SCZ_specific",
                     ifelse(!sig_scz & sig_bd,         "BD_specific", "not_significant"))))
  paired
}

# Jaccard index between the leading-edge genes of SCZ and BD for one pathway.
# Low overlap on a shared pathway = same function reached through different genes (hypothesis H1).
leading_edge_jaccard <- function(core_a, core_b) {
  a <- strsplit(core_a, "/", fixed = TRUE)[[1]]
  b <- strsplit(core_b, "/", fixed = TRUE)[[1]]
  length(intersect(a, b)) / length(union(a, b))
}

add_leading_edge_overlap <- function(paired) {
  paired$leading_edge_jaccard <- mapply(leading_edge_jaccard,
                                        paired$core_enrichment_SCZ,
                                        paired$core_enrichment_BD)
  paired
}
