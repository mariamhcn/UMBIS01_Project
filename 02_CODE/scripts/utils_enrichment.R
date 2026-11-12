# Helpers for Q3: functional enrichment (ORA then GSEA) of the Q2 pseudobulk DE results,
# with clusterProfiler.
#
# Gene sets are handled in two ways:
#  - GO Biological Process: clusterProfiler's own functions (enrichGO / gseGO) with the org.Hs.eg.db
#    annotation, as in the clusterProfiler book; this also allows simplify() to remove redundant GO terms;
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

# GSEA for one ranked list, returning the clusterProfiler object (used by run_gsea_all).
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
    # progress line in the render log (not in the HTML), to follow the long GSEA step
    cat(sprintf("[%s] GSEA %s - %s\n", format(Sys.time(), "%H:%M:%S"), ct, co), file = stderr())
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

# ---------------------------------------------------------------------------
# Figures of Part A (ORA)
# The order and the families of the cell types come from utils_figures.R (cell_family, order_cell_types).
# ---------------------------------------------------------------------------

ORA_LISTS <- c(SCZ_up = "SCZ up", SCZ_down = "SCZ down", BD_up = "BD up", BD_down = "BD down")

# Overview of Part A: (A) significant genes per cell type; (B) number of functions found by ORA in each list.
plot_ora_overview <- function(deg_counts, ora_table, cell_order) {
  centred <- ggplot2::theme(plot.title = ggplot2::element_text(colour = "black", face = "bold", size = 13, hjust = 0.5),
                            plot.title.position = "panel")
  a <- deg_counts |>
    dplyr::filter(disease %in% c("SCZ", "BD")) |>
    dplyr::select(cell_type, disease, up, down) |>
    tidyr::pivot_longer(c(up, down), names_to = "direction", values_to = "n") |>
    dplyr::mutate(signed = ifelse(direction == "up", n, -n),
                  cell_type = factor(cell_type, levels = rev(cell_order)),
                  disease = factor(disease, levels = c("SCZ", "BD"), labels = c("SCZ vs controls", "BD vs controls")))
  p_a <- ggplot2::ggplot(a, ggplot2::aes(signed, cell_type, fill = direction)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey30") +
    ggplot2::facet_wrap(~ disease) +
    ggplot2::scale_fill_manual(values = c(up = "#e66101", down = "#5e3c99"), breaks = c("up", "down"),
                               labels = c(up = "increased in patients", down = "decreased in patients"), name = NULL) +
    ggplot2::scale_x_continuous(labels = abs) +
    ggplot2::labs(x = "Number of significant genes", y = NULL, title = "Significant genes per cell type") +
    theme_project() + centred + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  b <- expand.grid(cell_type = cell_order, Cluster = names(ORA_LISTS), stringsAsFactors = FALSE) |>
    dplyr::left_join(dplyr::count(ora_table, cell_type, Cluster, name = "n_terms"), by = c("cell_type", "Cluster")) |>
    dplyr::mutate(n_terms = dplyr::coalesce(n_terms, 0L),
                  cell_type = factor(cell_type, levels = rev(cell_order)),
                  list = factor(ORA_LISTS[Cluster], levels = ORA_LISTS))
  p_b <- ggplot2::ggplot(b, ggplot2::aes(list, cell_type)) +
    ggplot2::geom_tile(ggplot2::aes(fill = n_terms > 0), colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = n_terms), size = 3.4) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "grey93", `TRUE` = "#e2a9f1"), guide = "none") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = NULL, y = NULL, title = "Functions found by ORA") +
    theme_project() + centred +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())
  (p_a + p_b + patchwork::plot_layout(widths = c(2.2, 1))) +
    patchwork::plot_annotation(tag_levels = "A", theme = theme_project(),
      title = "All cell types have significant genes, but ORA finds functions in only some of them",
      subtitle = paste0("(A) Number of genes significantly increased or decreased in patients (adjusted p < 0.05).\n",
                        "(B) Number of GO Biological Process terms over-represented in each of these lists (grey = none).")) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 15, hjust = 0, vjust = 1),
                   plot.tag.position = c(0, 1))
}

# One panel = one cell type: the n_top most significant functions of each of its four lists.
ora_panel <- function(d, ct, n_top, max_signif) {
  d <- d |>
    dplyr::group_by(Cluster) |>
    dplyr::slice_min(p.adjust, n = n_top, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(list = factor(ORA_LISTS[Cluster], levels = ORA_LISTS),
                  term = vapply(strwrap_list <- lapply(clean_set_name(Description, 80), strwrap, width = 34),
                                paste, character(1), collapse = "\n"),
                  gene_ratio = Count / DESize) |>
    dplyr::arrange(list, p.adjust)
  d$term <- factor(d$term, levels = rev(unique(d$term)))
  ggplot2::ggplot(d, ggplot2::aes(list, term)) +
    ggplot2::geom_point(ggplot2::aes(size = gene_ratio, colour = -log10(p.adjust))) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_colour_gradient(low = "#e2a9f1", high = "#5b1f75", limits = c(0, max_signif),
                                   name = "-log10(adjusted p)") +
    ggplot2::scale_size_continuous(range = c(2, 6), limits = c(0, 1), name = "Gene ratio") +
    ggplot2::labs(x = NULL, y = NULL, title = ct) +
    theme_project(10) +
    ggplot2::theme(plot.title = ggplot2::element_text(colour = "black", face = "bold", size = 13, hjust = 0.5),
                   plot.title.position = "panel",
                   axis.text.y = ggplot2::element_text(size = 9), axis.text.x = ggplot2::element_text(size = 9))
}

# One figure per family of cell types, one lettered panel (A, B, C...) per cell type with at least one function.
# Returns NULL when no cell type of the family has a result; otherwise the figure and a suitable height (inches).
ora_family_figure <- function(ora_table, family, collection_name, cell_order, n_top = 3, ncol = 2) {
  cts <- cell_order[cell_family(cell_order) == family & cell_order %in% ora_table$cell_type]
  if (length(cts) == 0) return(NULL)
  d_family   <- ora_table[ora_table$cell_type %in% cts, ]
  max_signif <- max(2, -log10(d_family$p.adjust))
  panels <- lapply(cts, function(ct) ora_panel(d_family[d_family$cell_type == ct, ], ct, n_top, max_signif))
  # height of a panel ~ number of functions shown; a row of panels is as tall as its tallest panel
  n_shown <- vapply(panels, function(p) max(2.5, 1.1 * nrow(p$data)), numeric(1))
  row_h   <- tapply(n_shown, ceiling(seq_along(n_shown) / ncol), max)
  p <- patchwork::wrap_plots(panels, ncol = ncol, heights = row_h) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A", theme = theme_project(),
      title = paste0(family, ": functions over-represented among the significant genes"),
      subtitle = paste0("ORA on ", collection_name, "; ", n_top,
                        " most significant functions per list (SCZ up, SCZ down, BD up, BD down); one panel per cell type"),
      caption = "Only cell types with at least one function are shown. Gene ratio = fraction of the list's genes that belong to the function. Darker = more significant.") &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 15, hjust = 0, vjust = 1),
                   plot.tag.position = c(0, 1), legend.position = "bottom")
  list(plot = p, height = 2.2 + 0.32 * sum(row_h) + 0.6 * length(row_h))
}

# Saves the family figures of one collection and returns the paths of the saved files, named by family
# (the notebook displays them in tabs).
save_ora_family_figures <- function(ora_table, collection_name, file_prefix, cell_order, fig_dir) {
  short <- c("Excitatory neurons" = "EN", "Inhibitory neurons" = "IN", "Glial cells" = "glia", "Vascular cells" = "vascular")
  paths <- character(0)
  for (family in FAMILIES) {
    fig <- ora_family_figure(ora_table, family, collection_name, cell_order)
    if (is.null(fig)) next
    path <- save_figure(fig$plot, paste0(file_prefix, "_", short[[family]]), fig_dir, width = 14, height = fig$height)
    paths <- c(paths, setNames(path, family))
  }
  paths
}

# Do SCZ and BD change the same genes? Share of each cell type's significant genes by group, with the Jaccard index.
plot_deg_overlap <- function(deg_overlap, cell_order) {
  groups <- c(SCZ_only = "SCZ only", shared = "Shared, same direction",
              opposite = "Shared, opposite directions", BD_only = "BD only")
  d <- deg_overlap |>
    tidyr::pivot_longer(dplyr::all_of(names(groups)), names_to = "group", values_to = "n") |>
    dplyr::mutate(fraction = n / n_DEG,
                  group = factor(groups[group], levels = groups),
                  cell_type = factor(cell_type, levels = rev(cell_order)))
  labels <- dplyr::mutate(deg_overlap, cell_type = factor(cell_type, levels = rev(cell_order)))
  ggplot2::ggplot(d, ggplot2::aes(fraction, cell_type, fill = group)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(data = labels, inherit.aes = FALSE, hjust = 0, size = 3,
                       ggplot2::aes(x = 1.02, y = cell_type,
                                    label = paste0("J = ", sprintf("%.2f", jaccard_DEG), "  (", n_DEG, " genes)"))) +
    ggplot2::scale_fill_manual(values = c("SCZ only" = "#2166ac", "Shared, same direction" = "#8e3fb0",
                                          "Shared, opposite directions" = "#f0a030", "BD only" = "#b2182b"), name = NULL) +
    ggplot2::scale_x_continuous(labels = function(x) paste0(round(100 * x), " %"), breaks = c(0, .25, .5, .75, 1),
                                expand = ggplot2::expansion(mult = c(0, 0.3))) +
    ggplot2::labs(x = "Share of the cell type's significant genes", y = NULL,
                  title = "SCZ and BD share few significant genes in every cell type",
                  subtitle = "Significant genes (adjusted p < 0.05) of each cell type, by disease",
                  caption = "J (Jaccard index) = genes shared in the same direction / all significant genes: 0 = no gene in common, 1 = same genes.") +
    theme_project() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

# ---------------------------------------------------------------------------
# Figures of Part B (GSEA)
# The order and the families of the cell types come from utils_figures.R (cell_family, order_cell_types).
# ---------------------------------------------------------------------------

GSEA_GROUP_COLOURS <- c("Shared, same direction" = "#8e3fb0", "Shared, opposite directions" = "#f0a030",
                        "SCZ only, confirmed" = "#2166ac", "SCZ only, not confirmed" = "#a6c8e6",
                        "BD only, confirmed" = "#b2182b", "BD only, not confirmed" = "#f0a8a8")

# Running enrichment score of one gene set along a ranked list (weighted Kolmogorov-Smirnov, as in GSEA, Subramanian et al. 2005).
running_score <- function(ranks, set_genes) {
  in_set <- names(ranks) %in% set_genes
  w <- abs(ranks) * in_set
  hit  <- cumsum(w) / sum(w)
  miss <- cumsum(!in_set) / sum(!in_set)
  hit - miss
}

# How to read a GSEA result: (A) running score, (B) position of the set's genes, (C) the ranked list.
plot_gsea_example <- function(ranks, set_genes, title, nes, padj) {
  rs <- running_score(ranks, set_genes)
  peak <- which.max(abs(rs))
  d <- data.frame(rank = seq_along(ranks), score = rs, stat = unname(ranks), in_set = names(ranks) %in% set_genes)
  up <- rs[peak] > 0
  le <- if (up) c(1, peak) else c(peak, length(ranks))
  shade <- ggplot2::annotate("rect", xmin = le[1], xmax = le[2], ymin = -Inf, ymax = Inf, fill = "#e2a9f1", alpha = 0.35)
  centred <- ggplot2::theme(plot.title = ggplot2::element_text(colour = "black", face = "bold", size = 12, hjust = 0.5),
                            plot.title.position = "panel")
  p_a <- ggplot2::ggplot(d, ggplot2::aes(rank, score)) + shade +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_line(colour = "#8e3fb0", linewidth = 1) +
    ggplot2::annotate("point", x = peak, y = rs[peak], size = 3, colour = "#5b1f75") +
    ggplot2::annotate("label", x = peak, y = rs[peak] * 0.85, hjust = -0.1, vjust = if (up) 1 else 0, size = 3.5,
                      label = paste0("peak = enrichment score\n-> NES = ", sprintf("%.2f", nes))) +
    ggplot2::labs(x = NULL, y = "Running score", title = "Running enrichment score") +
    theme_project() + centred + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  p_b <- ggplot2::ggplot(d[d$in_set, ], ggplot2::aes(x = rank)) + shade +
    ggplot2::geom_linerange(ggplot2::aes(ymin = 0, ymax = 1), linewidth = 0.3) +
    ggplot2::scale_x_continuous(limits = c(1, length(ranks))) +
    ggplot2::labs(x = NULL, y = NULL, title = "Position of the genes of the set in the ranking") +
    theme_project() + centred +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  p_c <- ggplot2::ggplot(d, ggplot2::aes(rank, stat)) +
    ggplot2::geom_area(fill = "grey70") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::annotate("text", x = length(ranks) * 0.02, y = max(d$stat) * 0.8, hjust = 0, size = 3.5,
                      label = "increased in patients") +
    ggplot2::annotate("text", x = length(ranks) * 0.98, y = min(d$stat) * 0.8, hjust = 1, size = 3.5,
                      label = "decreased in patients") +
    ggplot2::labs(x = "Rank of the gene (all tested genes)", y = "Wald statistic",
                  title = "All genes, ranked from the most increased to the most decreased") +
    theme_project() + centred
  (p_a / p_b / p_c + patchwork::plot_layout(heights = c(3, 1, 2))) +
    patchwork::plot_annotation(tag_levels = "A", theme = theme_project(), title = title,
      subtitle = paste0("NES = ", sprintf("%.2f", nes), ", adjusted p = ", signif(padj, 2),
                        ". Shaded: the leading edge (genes before the peak, which drive the result)")) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 15, hjust = 0, vjust = 1),
                   plot.tag.position = c(0, 1))
}

# Where does GSEA find something? Number of significant gene sets per cell type, up and down, SCZ and BD.
plot_gsea_counts <- function(gsea_all, collection_name, cell_order, alpha = 0.05) {
  d <- gsea_all |>
    dplyr::filter(collection == collection_name, disease %in% c("SCZ", "BD")) |>
    dplyr::group_by(cell_type, disease) |>
    dplyr::summarise(up = sum(p.adjust < alpha & NES > 0), down = sum(p.adjust < alpha & NES < 0),
                     n_tested = dplyr::n(), .groups = "drop") |>
    tidyr::pivot_longer(c(up, down), names_to = "direction", values_to = "n") |>
    dplyr::mutate(signed = ifelse(direction == "up", n, -n),
                  cell_type = factor(cell_type, levels = rev(cell_order)),
                  disease = factor(disease, levels = c("SCZ", "BD"), labels = c("SCZ vs controls", "BD vs controls")))
  ggplot2::ggplot(d, ggplot2::aes(signed, cell_type, fill = direction)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey30") +
    ggplot2::facet_wrap(~ disease) +
    ggplot2::scale_fill_manual(values = c(up = "#e66101", down = "#5e3c99"), breaks = c("up", "down"),
                               labels = c(up = "shifted up (NES > 0)", down = "shifted down (NES < 0)"), name = NULL) +
    ggplot2::scale_x_continuous(labels = abs) +
    ggplot2::labs(x = "Number of significant gene sets (adjusted p < 0.05)", y = NULL,
                  title = paste0("Gene sets shifted in patients, per cell type (", collection_name, ")"),
                  subtitle = paste0("GSEA on all genes; ", max(d$n_tested), " gene sets tested per cell type at most"),
                  caption = "Right of the line: sets whose genes are grouped among the genes increased in patients; left: among the decreased genes.") +
    theme_project() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

# Which gene sets? NES heatmap of the gene sets significant in the largest number of cell types.
# One lettered panel per comparison (SCZ, BD and, if asked, the direct SCZ vs BD test), stacked or side by side;
# the cell types are grouped by family and the rows are in the same order in every panel.
plot_gsea_heatmap <- function(gsea_all, collection_name, cell_order, n_sets = 30, alpha = 0.05,
                              diseases = c("SCZ", "BD"), side_by_side = TRUE, title = NULL, subtitle = NULL,
                              caption = NULL) {
  d_all <- dplyr::filter(gsea_all, collection == collection_name)
  ranking <- d_all |>
    dplyr::filter(disease %in% c("SCZ", "BD")) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(n_sig = sum(p.adjust < alpha), best = min(p.adjust), .groups = "drop") |>
    dplyr::filter(n_sig > 0) |>
    dplyr::arrange(dplyr::desc(n_sig), best)
  keep <- utils::head(ranking$ID, n_sets)
  # rows: sets with similar SCZ and BD profiles next to each other
  wide <- d_all |>
    dplyr::filter(ID %in% keep, disease %in% c("SCZ", "BD")) |>
    dplyr::mutate(column = paste(disease, cell_type)) |>
    dplyr::select(ID, column, NES) |>
    tidyr::pivot_wider(names_from = column, values_from = NES, values_fill = 0)
  row_order <- if (nrow(wide) > 2) wide$ID[stats::hclust(stats::dist(as.matrix(wide[, -1])))$order] else wide$ID
  d <- dplyr::filter(d_all, ID %in% keep, disease %in% diseases)
  labels <- stats::setNames(clean_set_name(d$Description, 50), d$ID)
  lim <- max(abs(d$NES), na.rm = TRUE)
  short_family <- c("Excitatory neurons" = "Excitatory neurons", "Inhibitory neurons" = "Inhibitory neurons",
                    "Glial cells" = "Glia", "Vascular cells" = "Vasc.")
  panel_names <- c(SCZ = "SCZ vs controls", BD = "BD vs controls", "SCZ vs BD" = "SCZ vs BD (direct comparison)")
  d <- d |> dplyr::mutate(ID = factor(ID, levels = row_order),
                          cell_type = factor(cell_type, levels = cell_order),
                          family = factor(short_family[cell_family(as.character(cell_type))], levels = short_family),
                          star = ifelse(p.adjust < alpha, "*", ""))
  one <- function(dis, show_y) {
    ggplot2::ggplot(dplyr::filter(d, disease == dis), ggplot2::aes(cell_type, ID, fill = NES)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = star), size = 4, vjust = 0.75) +
      ggplot2::facet_grid(. ~ family, scales = "free_x", space = "free_x") +
      ggplot2::scale_fill_gradient2(low = "#5e3c99", mid = "white", high = "#e66101", midpoint = 0,
                                    limits = c(-lim, lim), name = "NES") +
      ggplot2::scale_y_discrete(labels = labels, drop = FALSE) +
      ggplot2::labs(x = NULL, y = NULL, title = panel_names[[dis]]) +
      theme_project() +
      ggplot2::theme(plot.title = ggplot2::element_text(colour = "black", face = "bold", size = 13, hjust = 0.5),
                     plot.title.position = "panel",
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                     axis.text.y = if (show_y) ggplot2::element_text() else ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(), panel.spacing.x = ggplot2::unit(2, "pt"),
                     strip.text = ggplot2::element_text(size = 8))
  }
  panels <- lapply(seq_along(diseases), function(k) one(diseases[k], show_y = !side_by_side || k == 1))
  layout <- if (side_by_side) patchwork::plot_layout(nrow = 1, guides = "collect") else patchwork::plot_layout(ncol = 1, guides = "collect")
  tags <- if (length(diseases) > 1) "A" else NULL
  (patchwork::wrap_plots(panels) + layout) +
    patchwork::plot_annotation(tag_levels = tags, theme = theme_project(),
      title = if (is.null(title)) paste0("Gene sets shifted in SCZ and BD (", collection_name, ")") else title,
      subtitle = if (is.null(subtitle)) paste0("The ", length(keep), " gene sets significant in the largest number of cell types; * adjusted p < ", alpha) else subtitle,
      caption = if (is.null(caption)) paste0("Orange: genes of the set grouped among the genes increased in patients; purple: among the decreased genes.",
                                             if ("SCZ vs BD" %in% diseases) " Direct comparison: orange = higher in SCZ than in BD, purple = higher in BD." else "") else caption) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 15, hjust = 0, vjust = 1),
                   plot.tag.position = c(0, 1), legend.position = "right")
}

# Together or separately? (A) category of the significant gene sets, per cell type;
# (B) similarity of the whole profiles (Spearman rho of the NES of all sets in SCZ and BD).
plot_shared_specific <- function(paired, collection_name, cell_order) {
  centred <- ggplot2::theme(plot.title = ggplot2::element_text(colour = "black", face = "bold", size = 13, hjust = 0.5),
                            plot.title.position = "panel")
  a <- paired |>
    dplyr::filter(collection == collection_name, group != "Not significant") |>
    dplyr::count(cell_type, group) |>
    dplyr::mutate(group = factor(group, levels = names(GSEA_GROUP_COLOURS)),
                  cell_type = factor(cell_type, levels = rev(cell_order)))
  p_a <- ggplot2::ggplot(a, ggplot2::aes(n, cell_type, fill = group)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_fill_manual(values = GSEA_GROUP_COLOURS, name = NULL, drop = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::labs(x = "Number of significant gene sets", y = NULL,
                  title = "Shared or specific? (significant sets)") +
    theme_project() + centred + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  b <- paired |>
    dplyr::group_by(collection, cell_type) |>
    dplyr::summarise(rho = stats::cor(NES_SCZ, NES_BD, method = "spearman"), .groups = "drop") |>
    dplyr::filter(collection %in% c("GO_BP", "Hallmark", "Reactome")) |>
    dplyr::mutate(cell_type = factor(cell_type, levels = rev(cell_order)),
                  collection = factor(collection, levels = c("Hallmark", "GO_BP", "Reactome")))
  p_b <- ggplot2::ggplot(b, ggplot2::aes(rho, cell_type, colour = collection, shape = collection)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey50") +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::scale_colour_manual(values = c(Hallmark = "#8e3fb0", GO_BP = "#1b8a78", Reactome = "#d46a00"), name = NULL) +
    ggplot2::scale_shape_manual(values = c(Hallmark = 16, GO_BP = 17, Reactome = 15), name = NULL) +
    ggplot2::scale_x_continuous(limits = c(-0.2, 1), breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(x = "Spearman rho of the NES (SCZ vs BD)", y = NULL,
                  title = "Similar profiles? (all sets)") +
    theme_project() + centred +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(colour = "grey93"))
  (p_a + p_b + patchwork::plot_layout(widths = c(1.4, 1))) +
    patchwork::plot_annotation(tag_levels = "A", theme = theme_project(),
      title = "Shared and disease-specific gene sets, and similarity of the SCZ and BD profiles",
      subtitle = paste0("(A) ", collection_name, " gene sets significant in at least one disease, by category. ",
                        "(B) Correlation of the NES of all tested sets between SCZ and BD (1 = same profile, 0 = unrelated)."),
      caption = "Confirmed = also significant in the direct SCZ vs BD test; light colours = seen in one disease only, without evidence of a difference.") &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 15, hjust = 0, vjust = 1),
                   plot.tag.position = c(0, 1), legend.position = "bottom")
}

# Same pathway, same genes? Overlap (Jaccard) of the genes driving each shared pathway in SCZ and BD,
# compared with the overlap of the significant genes (diamond, Step A4).
# Overlap expected by chance: if each disease picked its driving genes at random inside the set
# (a genes in SCZ, b in BD, among n genes of the set), they would share a * b / n genes on average.
add_expected_jaccard <- function(shared, gsea_all) {
  sizes <- gsea_all |>
    dplyr::filter(disease == "SCZ") |>
    dplyr::select(collection, cell_type, ID, setSize)
  shared |>
    dplyr::left_join(sizes, by = c("collection", "cell_type", "ID")) |>
    dplyr::mutate(a = lengths(strsplit(core_enrichment_SCZ, "/")),
                  b = lengths(strsplit(core_enrichment_BD, "/")),
                  expected_jaccard = (a * b / setSize) / (a + b - a * b / setSize))
}

# Paired test per cell type: in its shared pathways, does the observed overlap differ from the overlap expected
# by chance? Wilcoxon signed-rank test on the pairs (observed, expected) of each pathway, two-sided,
# Benjamini-Hochberg correction over the cell types.
test_leading_edge <- function(shared) {
  shared |>
    dplyr::filter(collection %in% c("GO_BP", "Hallmark", "Reactome")) |>
    dplyr::group_by(cell_type) |>
    dplyr::summarise(n_pathways = dplyr::n(),
                     median_observed = stats::median(leading_edge_jaccard),
                     median_expected = stats::median(expected_jaccard),
                     share_above = mean(leading_edge_jaccard > expected_jaccard),
                     p_value = suppressWarnings(stats::wilcox.test(leading_edge_jaccard, expected_jaccard,
                                                                   paired = TRUE)$p.value),
                     .groups = "drop") |>
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH"))
}

plot_leading_edge <- function(shared, le_test, cell_order, alpha = 0.05, seed = 42) {
  d <- dplyr::filter(shared, collection %in% c("GO_BP", "Hallmark", "Reactome"))
  cols <- c("Expected by chance" = "grey55", "Observed" = "#8e3fb0")
  long <- d |>
    tidyr::pivot_longer(c(expected_jaccard, leading_edge_jaccard), names_to = "what", values_to = "J") |>
    dplyr::mutate(what = factor(what, levels = c("expected_jaccard", "leading_edge_jaccard"), labels = names(cols)),
                  family = cell_family(cell_type))
  # one bracket per cell type, just above its highest dot, with the result of the test
  tops <- long |> dplyr::group_by(cell_type) |> dplyr::summarise(top = max(J, na.rm = TRUE), .groups = "drop")
  brackets <- le_test |>
    dplyr::left_join(tops, by = "cell_type") |>
    dplyr::mutate(label = dplyr::case_when(p_adj < 0.0001 ~ "****", p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                                           p_adj < alpha ~ "*", TRUE ~ "ns"),
                  y = top + 0.04, family = cell_family(cell_type))
  x_labels <- stats::setNames(paste0(le_test$cell_type, "\n(", le_test$n_pathways, ")"), le_test$cell_type)
  y_max <- max(brackets$y) + 0.06

  one <- function(fam) {
    cts <- cell_order[cell_order %in% le_test$cell_type & cell_family(cell_order) == fam]
    l <- dplyr::filter(long, family == fam) |> dplyr::mutate(cell_type = factor(cell_type, levels = cts))
    b <- dplyr::filter(brackets, family == fam) |>
      dplyr::mutate(pos = match(cell_type, cts))   # position of the cell type on the axis
    ggplot2::ggplot(l, ggplot2::aes(cell_type, J, colour = what)) +
      ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8, seed = seed),
                          size = 0.7, alpha = 0.3) +
      ggplot2::geom_boxplot(fill = NA, outlier.shape = NA, width = 0.7, linewidth = 0.5,
                            position = ggplot2::position_dodge(width = 0.8)) +
      ggplot2::geom_segment(data = b, inherit.aes = FALSE, linewidth = 0.4,
                            ggplot2::aes(x = pos - 0.2, xend = pos + 0.2, y = y, yend = y)) +
      ggplot2::geom_segment(data = b, inherit.aes = FALSE, linewidth = 0.4,
                            ggplot2::aes(x = pos - 0.2, xend = pos - 0.2, y = y, yend = y - 0.015)) +
      ggplot2::geom_segment(data = b, inherit.aes = FALSE, linewidth = 0.4,
                            ggplot2::aes(x = pos + 0.2, xend = pos + 0.2, y = y, yend = y - 0.015)) +
      ggplot2::geom_text(data = b, inherit.aes = FALSE, size = 3.3,
                         ggplot2::aes(x = pos, y = y + 0.012, label = label, vjust = ifelse(label == "ns", 0, 0.4))) +
      ggplot2::scale_colour_manual(values = cols, name = NULL) +
      ggplot2::scale_x_discrete(labels = x_labels) +
      ggplot2::scale_y_continuous(limits = c(0, y_max), breaks = seq(0, 1, 0.25)) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1, shape = 15, linetype = 0))) +
      ggplot2::labs(x = NULL, y = "Jaccard index", title = fam) +
      theme_project() +
      ggplot2::theme(plot.title = ggplot2::element_text(colour = "black", face = "bold", size = 13, hjust = 0.5),
                     plot.title.position = "panel",
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  n_ct <- function(fam) sum(cell_order %in% le_test$cell_type & cell_family(cell_order) == fam)
  p_a <- one("Excitatory neurons"); p_b <- one("Inhibitory neurons")
  p_c <- one("Glial cells");        p_d <- one("Vascular cells")
  # design of the figure: excitatory and inhibitory neurons on the first row, glia and vascular cells on the second;
  # each panel is as wide as its number of cell types
  n <- vapply(FAMILIES, n_ct, numeric(1))
  # excitatory neurons alone on the first row; inhibitory neurons, glia and vascular cells below
  design <- paste0(strrep("A", sum(n[2:4])), "\n", strrep("B", n[2]), strrep("C", n[3]), strrep("D", n[4]))
  (p_a + p_b + p_c + p_d + patchwork::plot_layout(design = design, guides = "collect")) +
    patchwork::plot_annotation(tag_levels = "A", theme = theme_project(),
      title = "The pathways shared by SCZ and BD are driven by partly the same genes",
      subtitle = "One dot per shared pathway: overlap of the genes that drive it in SCZ and in BD (Jaccard index), observed and expected by chance",
      caption = paste0("Jaccard index: 0 = no gene in common, 1 = the same genes. In brackets under each cell type: number of pathways shared by SCZ and BD ",
                       "(gene sets shifted in the same direction in both diseases; GO, Hallmark, Reactome).\n",
                       "Observed vs expected: paired two-sided Wilcoxon signed-rank test per cell type, Benjamini-Hochberg adjusted (ns: not significant; * < 0.05; ** < 0.01; *** < 0.001; **** < 0.0001).")) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 15, hjust = 0, vjust = 1),
                   plot.tag.position = c(0, 1), legend.position = "bottom")
}

# Writes the markdown of a tabset (one tab per element of a named list of PNG paths).
# To use in a chunk with results = "asis", inside a ::: {.panel-tabset} block.
show_figures_in_tabs <- function(files_by_tab) {
  for (tab in names(files_by_tab)) {
    cat("\n## ", tab, " {.unnumbered .unlisted}\n\n", sep = "")
    if (length(files_by_tab[[tab]]) == 0) cat("No result.\n\n")
    for (f in files_by_tab[[tab]]) cat("![](", f, ")\n\n", sep = "")
  }
}
