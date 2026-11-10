# Gene set helpers shared by Q3 (enrichment) and Q5 (ageing scores).
# Every gene set is stored as a TERM2GENE data frame: columns gs_name, gene_symbol.

# Read a gene set file of 01_DATA/genesets/ (CSV with gs_name, gene_symbol).
# Returns NULL with a warning if the file is not there yet, so notebooks keep running.
read_geneset_file <- function(path) {
  if (!file.exists(path)) {
    warning("Gene set file not found, skipped: ", path)
    return(NULL)
  }
  sets <- utils::read.csv(path, stringsAsFactors = FALSE)
  stopifnot(all(c("gs_name", "gene_symbol") %in% colnames(sets)))
  unique(sets[, c("gs_name", "gene_symbol")])
}

# Convert a TERM2GENE data frame into a named list of gene vectors.
term2gene_to_list <- function(term2gene) {
  split(term2gene$gene_symbol, term2gene$gs_name)
}

# Fraction of each gene set found in the data: a set with poor coverage gives unreliable results.
geneset_coverage <- function(set_list, genes_in_data) {
  data.frame(
    gs_name   = names(set_list),
    n_genes   = vapply(set_list, length, integer(1)),
    n_present = vapply(set_list, function(g) sum(g %in% genes_in_data), integer(1)),
    row.names = NULL
  ) |>
    transform(fraction_present = round(n_present / n_genes, 2))
}
