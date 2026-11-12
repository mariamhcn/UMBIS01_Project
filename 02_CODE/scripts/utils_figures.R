# Shared figure helpers: one look for every figure of the project (notebooks, slides, report).
#
#  - theme_project(): same fonts, sizes and layout everywhere;
#  - save_figure(): saves each figure as a PNG (300 dpi), for the slides and the report;
#  - clean_set_name(): turns gene set identifiers (HALLMARK_..., GOBP_...) into readable labels;
#  - cell_family(), order_cell_types(): the same order of the cell types in every figure of the project;
#  - plot_reasoning_diagram(): the diagram that shows how a notebook answers its question (Overview).
#
# This file holds only the common style, shared by all notebooks. The figures specific to one question are drawn by
# that question's script (utils_enrichment.R for Q3, utils_aging_scores.R for Q5).

# Project colours (pastel = backgrounds, dark = accents), as in the slides and the report
project_colours <- list(
  preprocessing = c(light = "#cfd6f8", dark = "#3849a2"),
  Q1 = c(light = "#a8e6dc", dark = "#1b8a78"),
  Q2 = c(light = "#ffcf9e", dark = "#d46a00"),
  Q3 = c(light = "#e2a9f1", dark = "#8e3fb0"),
  Q4 = c(light = "#7ed957", dark = "#2f7d12"),
  Q5 = c(light = "#fbe38c", dark = "#a88600")
)

# Disease colours, identical in every notebook
col_disease <- c(SCZ = "#2166ac", BD = "#b2182b", "SCZ vs BD" = "grey40")

theme_project <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title          = ggplot2::element_text(face = "bold", size = base_size + 2, colour = "#3849a2"),
      plot.subtitle       = ggplot2::element_text(size = base_size, colour = "grey25"),
      plot.caption        = ggplot2::element_text(size = base_size - 2, colour = "grey40", hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      strip.background    = ggplot2::element_rect(fill = "grey95", colour = "grey70"),
      strip.text          = ggplot2::element_text(face = "bold"),
      legend.position     = "bottom",
      panel.grid.minor    = ggplot2::element_blank()
    )
}

# Save one figure as <name>.png (300 dpi, white background).
# width / height in inches: 10 x 6 fits a 16:9 slide and the width of a report page.
save_figure <- function(p, name, dir, width = 10, height = 6, dpi = 300) {
  path <- file.path(dir, paste0(name, ".png"))
  ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi, bg = "white")
  invisible(path)
}

# Same cleaning inside a clusterProfiler result, before plotting it with enrichplot (dotplot, cnetplot...).
clean_result_names <- function(x) {
  if (methods::.hasSlot(x, "compareClusterResult")) {
    x@compareClusterResult$Description <- clean_set_name(x@compareClusterResult$Description)
  } else if (methods::.hasSlot(x, "result")) {
    x@result$Description <- clean_set_name(x@result$Description)
  }
  x
}

# "HALLMARK_INFLAMMATORY_RESPONSE" -> "Inflammatory response"; GO descriptions are kept as they are.
clean_set_name <- function(x, max_char = 60) {
  out <- ifelse(grepl("^[A-Z0-9_]+$", x),
                gsub("_", " ", sub("^(HALLMARK|GOBP|GOCC|GOMF|REACTOME|KEGG|SYNGO|WP)_", "", x)),
                x)
  was_upper <- grepl("^[A-Z0-9 ]+$", out)
  out <- ifelse(was_upper, tolower(out), out)
  # put common acronyms back in capitals
  acronyms <- c("TNFA", "NFKB", "MTORC1", "MTOR", "KRAS", "IL6", "IL2", "JAK", "STAT3", "STAT5", "DNA", "RNA",
                "ATP", "TCA", "G2M", "E2F", "MYC", "UV", "P53", "TGF", "PI3K", "AKT", "WNT", "NMDA", "AMPA",
                "GABA", "ROS", "UPR", "MHC", "ECM", "GPCR", "ERK", "MAPK", "CREB",
                "FGFR1", "FGFR2", "FGFR3", "FGFR4", "TP53", "RHOA", "BMP", "NTRK", "NTRKS", "SREBF", "SREBP",
                "PDH", "EIF2AK4", "GCN2", "SRP", "HSF1", "TNFR1", "NOTCH3", "MET", "RET")
  for (a in acronyms) out[was_upper] <- gsub(paste0("\\b", tolower(a), "\\b"), a, out[was_upper])
  out[was_upper] <- gsub("\\bsars cov 2\\b", "SARS-CoV-2", out[was_upper])
  out[was_upper] <- gsub("\\bnfkappab\\b", "NF-kappaB", out[was_upper])
  out[was_upper] <- gsub("\\bgtpase\\b", "GTPase", out[was_upper])
  out <- paste0(toupper(substr(out, 1, 1)), substr(out, 2, nchar(out)))
  ifelse(nchar(out) > max_char, paste0(substr(out, 1, max_char - 3), "..."), out)
}

# ---------------------------------------------------------------------------------------
# Cell types: families and order, the same in every figure of the project
# ---------------------------------------------------------------------------------------

# The four families of cell types (subclasses), in the order used in the figures.
FAMILIES <- c("Excitatory neurons", "Inhibitory neurons", "Glial cells", "Vascular cells")

# Family of each cell type: EN_* = excitatory neurons, IN_* = inhibitory neurons,
# Astro / Micro / Oligo / OPC = glial cells, the others (Endo, PC, SMC, VLMC...) = vascular cells.
cell_family <- function(x) {
  ifelse(grepl("^EN_", x), "Excitatory neurons",
  ifelse(grepl("^IN_", x), "Inhibitory neurons",
  ifelse(x %in% c("Astro", "Micro", "Oligo", "OPC"), "Glial cells", "Vascular cells")))
}

# Excitatory neurons, then inhibitory neurons, glial cells and vascular cells; alphabetical inside a family.
# The alphabetical order is the same on every computer (as sort_subclasses() in utils_h5ad.R): case is ignored and
# "_" comes before letters and digits, so EN_L6_CT comes before EN_L6B whatever the language settings.
order_cell_types <- function(x) {
  u <- unique(as.character(x))
  u[order(match(cell_family(u), FAMILIES), gsub("_", " ", toupper(u)), method = "radix")]
}

# ---------------------------------------------------------------------------
# Reasoning diagram of a notebook (Overview)
# ---------------------------------------------------------------------------

# Columns of boxes linked by arrows, a conclusion box at the bottom.
# cols: data.frame(head, sub, fill, border); boxes: data.frame(col, step, text, out)
plot_reasoning_diagram <- function(cols, boxes, conclusion, title, subtitle, w = 0.42, h = 0.41) {
  cols$x <- seq_len(nrow(cols))
  boxes <- merge(boxes, cols[, c("x", "fill", "border")], by.x = "col", by.y = "x")
  boxes$x <- boxes$col
  boxes <- do.call(rbind, lapply(split(boxes, boxes$col), function(d) { d$y <- rev(seq_len(nrow(d))) + (4 - nrow(d)); d }))
  ymax <- 4
  arrows_in <- do.call(rbind, lapply(split(boxes, boxes$col), function(d) if (nrow(d) > 1) data.frame(x = d$x[-nrow(d)], y = d$y[-nrow(d)]) else NULL))
  n <- nrow(cols)
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = arrows_in, ggplot2::aes(x = x, xend = x, y = y - h, yend = y - 1 + h),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"), type = "closed"), colour = "grey45") +
    (if (n > 1) ggplot2::annotate("segment", x = seq_len(n - 1) + w + 0.01, xend = seq_len(n - 1) + 1 - w - 0.01,
                                  y = ymax + 0.64, yend = ymax + 0.64,
                                  arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm"), type = "closed"),
                                  colour = "grey30", linewidth = 1.2)) +
    ggplot2::geom_rect(data = boxes, ggplot2::aes(xmin = x - w, xmax = x + w, ymin = y - h, ymax = y + h,
                                                  fill = I(fill), colour = I(border)), linewidth = 0.6) +
    ggplot2::geom_text(data = boxes, ggplot2::aes(x = x - w + 0.03, y = y + h - 0.06, label = step),
                       hjust = 0, vjust = 1, fontface = "bold", size = 4.4) +
    ggplot2::geom_text(data = boxes, ggplot2::aes(x = x, y = y + 0.04, label = text), size = 4.6, lineheight = 0.95) +
    ggplot2::geom_label(data = subset(boxes, out != ""), ggplot2::aes(x = x + w - 0.03, y = y - h + 0.06, label = out),
                        hjust = 1, vjust = 0, size = 4, fontface = "bold", fill = "white", label.size = 0.3) +
    ggplot2::geom_text(data = cols, ggplot2::aes(x = x, y = ymax + 0.75, label = head), fontface = "bold", size = 5.6) +
    ggplot2::geom_text(data = cols, ggplot2::aes(x = x, y = ymax + 0.53, label = sub), size = 4.3, colour = "grey35") +
    ggplot2::annotate("rect", xmin = 1 - w, xmax = n + w, ymin = -0.05, ymax = 0.4, fill = "grey95", colour = "grey40", linewidth = 0.5) +
    ggplot2::annotate("text", x = (n + 1) / 2, y = 0.175, size = 4.6, label = conclusion) +
    ggplot2::coord_cartesian(xlim = c(1 - w - 0.03, n + w + 0.03), ylim = c(-0.1, ymax + 0.9), expand = FALSE) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", colour = "#3849a2", size = 18),
                   plot.subtitle = ggplot2::element_text(size = 13.5, colour = "grey30"),
                   plot.margin = ggplot2::margin(10, 10, 10, 10), plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}
