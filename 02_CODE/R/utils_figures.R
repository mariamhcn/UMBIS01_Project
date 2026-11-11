# Shared figure helpers: one look for every figure of the project (notebooks, slides, report).
#
#  - theme_project(): same fonts, sizes and layout everywhere;
#  - save_figure(): saves each figure as a PNG (300 dpi), for the slides and the report;
#  - clean_set_name(): turns gene set identifiers (HALLMARK_..., GOBP_...) into readable labels.

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
                "GABA", "ROS", "UPR", "MHC", "ECM", "GPCR", "ERK", "MAPK", "CREB")
  for (a in acronyms) out[was_upper] <- gsub(paste0("\\b", tolower(a), "\\b"), a, out[was_upper])
  out <- paste0(toupper(substr(out, 1, 1)), substr(out, 2, nchar(out)))
  ifelse(nchar(out) > max_char, paste0(substr(out, 1, max_char - 3), "..."), out)
}
