# Functions to read the data file of the HBCC cohort (.h5ad).
# This file exists because the data set is too large to be loaded at once:
# the functions read it one donor at a time.
# Used by notebook 01_preprocessing.qmd.
#
# Technical notes, for those who want to modify these functions.
# An .h5ad file is an HDF5 file with a fixed structure (AnnData):
#   - /obs : one row per nucleus (donor, diagnosis, age, subclass...);
#   - /var : one row per gene (Ensembl identifier, gene symbol);
#   - /X   : the count matrix, stored in CSR format (one row per nucleus) as three vectors:
#            data (non-zero values), indices (gene number of each value, from 0) and
#            indptr (where each nucleus starts in data and indices, from 0).
# Text columns are stored as "categorical": integer codes (from 0) + a dictionary of categories.
# The full matrix (1,486,324 nuclei, 6.1 billion values) needs about 50 GB of memory,
# so it is always read one block of nuclei (one donor) at a time.
#
# Author of the reading method: Myriam (scripts 01, 02 and 04 of the preprocessing).

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
})

# ---------------------------------------------------------------------------------------
# Metadata (obs) and gene names (var)
# ---------------------------------------------------------------------------------------

# Reads one column of a data frame stored in the file (obs or var).
# A categorical column becomes a factor, with the categories in the order of the file.
read_h5ad_column <- function(file, group, column) {
  path <- paste0("/", group, "/", column)
  fid <- H5Fopen(file, flags = "H5F_ACC_RDONLY")
  # categorical = a group holding codes + categories (a plain column is a dataset)
  is_categorical <- tryCatch(H5Lexists(fid, paste0(path, "/codes")), error = function(e) FALSE)
  H5Fclose(fid)
  if (is_categorical) {
    codes      <- as.integer(h5read(file, paste0(path, "/codes")))
    categories <- as.character(h5read(file, paste0(path, "/categories")))
    codes[codes < 0] <- NA                             # -1 means missing value
    factor(categories[codes + 1], levels = categories)
  } else {
    as.vector(h5read(file, path))
  }
}

# Reads the nucleus metadata: one row per nucleus, in the order of the file.
read_h5ad_obs <- function(file, columns) {
  content <- h5ls(file, recursive = 2, datasetinfo = FALSE)
  obs_columns <- content$name[content$group == "/obs"]
  missing_columns <- setdiff(columns, obs_columns)
  if (length(missing_columns) > 0) stop("Columns not found in /obs: ", paste(missing_columns, collapse = ", "))
  meta <- lapply(columns, function(col) read_h5ad_column(file, "obs", col))
  names(meta) <- columns
  h5closeAll()
  as.data.frame(meta, stringsAsFactors = FALSE)
}

# Gene symbols (made unique, as Seurat requires) and Ensembl identifiers.
read_gene_names <- function(file) {
  symbols <- as.character(read_h5ad_column(file, "var", "feature_name"))
  ensembl <- as.character(h5read(file, "/var/_index"))
  h5closeAll()
  data.frame(ensembl_id = ensembl, symbol = symbols, name = make.unique(symbols))
}

# ---------------------------------------------------------------------------------------
# Count matrix (X), one block of nuclei at a time
# ---------------------------------------------------------------------------------------

# indptr exceeds 2 billion: it is read as doubles so that it is not truncated.
read_indptr <- function(file) {
  as.numeric(h5read(file, "/X/indptr", bit64conversion = "double"))
}

# Reads the nuclei first..last (row numbers in the file, from 1) as a sparse matrix nuclei x genes.
read_nuclei_block <- function(file, indptr, first, last, n_genes) {
  start    <- indptr[first]                            # position (from 0) of the first count
  end      <- indptr[last + 1]                         # end position (excluded)
  n_values <- end - start
  values <- h5read(file, "/X/data",    start = start + 1, count = n_values)
  genes  <- h5read(file, "/X/indices", start = start + 1, count = n_values)
  new("dgRMatrix",
      p   = as.integer(indptr[first:(last + 1)] - start),
      j   = as.integer(genes),
      x   = as.numeric(values),
      Dim = c(as.integer(last - first + 1), as.integer(n_genes)))
}

# Joins sparse matrices side by side (same rows), without cbind().
# cbind() joins the pieces two by two and copies the whole matrix at each step (9.8 GB of 10 at
# the first attempt). A sparse matrix is only three vectors (x: values, i: row of each value,
# p: where each column starts), so the vectors are put end to end and p is shifted.
join_columns <- function(pieces) {
  n_values <- vapply(pieces, function(m) length(m@x), numeric(1))
  if (sum(n_values) >= 2^31) stop("Too many non-zero values for one sparse matrix")
  offset <- c(0, cumsum(n_values))
  p_all <- c(0L, as.integer(unlist(lapply(seq_along(pieces), function(k) pieces[[k]]@p[-1] + offset[k]))))
  out <- new("dgCMatrix",
             i   = unlist(lapply(pieces, function(m) m@i)),
             p   = p_all,
             x   = unlist(lapply(pieces, function(m) m@x)),
             Dim = c(nrow(pieces[[1]]), length(p_all) - 1L))
  validObject(out)
  out
}

# Row numbers of the nuclei of each donor. The method needs the nuclei of a donor to follow
# each other in the file: this is checked here.
rows_per_donor <- function(meta, donor_ids) {
  rows <- split(seq_len(nrow(meta)), as.character(meta$donor_id))[donor_ids]
  contiguous <- vapply(rows, function(r) all(diff(r) == 1), logical(1))
  if (!all(contiguous)) stop("The nuclei of some donors do not follow each other in the file")
  rows
}

# Prints the progress of a loop over donors every 20 donors.
report_progress <- function(k, n, start_time) {
  if (k %% 20 == 0 || k == n) {
    elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
    cat(k, "donors of", n, "-", elapsed, "min\n")
  }
}

# ---------------------------------------------------------------------------------------
# Random subsample of nuclei (Seurat object)
# ---------------------------------------------------------------------------------------

# Draws at most n_per_donor nuclei per donor and returns their row numbers (sorted, per donor).
draw_nuclei <- function(rows_by_donor, n_per_donor, seed) {
  set.seed(seed)
  lapply(rows_by_donor, function(r) sort(sample(r, min(length(r), n_per_donor))))
}

# Reads the chosen nuclei, donor by donor: genes x nuclei, raw counts.
read_chosen_nuclei <- function(file, chosen_by_donor, gene_names) {
  indptr  <- read_indptr(file)
  n_genes <- length(gene_names)
  pieces  <- vector("list", length(chosen_by_donor))
  start_time <- Sys.time()
  for (k in seq_along(chosen_by_donor)) {
    chosen <- chosen_by_donor[[k]]
    block  <- read_nuclei_block(file, indptr, min(chosen), max(chosen), n_genes)
    kept   <- block[chosen - min(chosen) + 1, , drop = FALSE]
    pieces[[k]] <- t(as(kept, "CsparseMatrix"))          # genes x nuclei
    rm(block, kept); invisible(gc())                     # R does not free the discarded block on its own
    report_progress(k, length(chosen_by_donor), start_time)
  }
  h5closeAll()
  counts <- join_columns(pieces)
  rownames(counts) <- gene_names
  colnames(counts) <- paste0("nucleus_", unlist(chosen_by_donor))
  counts
}

# ---------------------------------------------------------------------------------------
# Pseudobulk: sum of the counts of all nuclei of each subclass, per donor
# ---------------------------------------------------------------------------------------

# For one donor: an indicator matrix (subclasses x nuclei, 1 when the nucleus belongs to the
# subclass) multiplied by the block (nuclei x genes) adds up the nuclei of each subclass.
# All subclasses are kept, even those absent in this donor (empty columns).
pseudobulk_block <- function(block, subclass_of_nuclei, subclasses) {
  indicator <- fac2sparse(factor(subclass_of_nuclei, levels = subclasses), drop.unused.levels = FALSE)
  t(as(indicator %*% block, "CsparseMatrix"))           # genes x subclasses
}

# Alphabetical order of the subclasses, the same on every computer: sort() depends on the language
# settings of the session (EN_L6B before or after EN_L6_CT). Here, case is ignored and "_" comes
# before letters and digits, as on the computers where the first pseudobulk matrix was made.
sort_subclasses <- function(x) {
  x[order(gsub("_", " ", toupper(x)), method = "radix")]
}

# Pseudobulk matrix of the selected donors: genes x (donor x subclass), raw counts,
# with a table describing each column.
build_pseudobulk <- function(file, meta, donors, gene_names) {
  subclasses <- sort_subclasses(unique(as.character(meta$subclass)))
  rows       <- rows_per_donor(meta, donors$donor_id)
  indptr     <- read_indptr(file)
  n_genes    <- length(gene_names)
  pieces     <- vector("list", length(rows))
  start_time <- Sys.time()
  for (k in seq_along(rows)) {
    r <- rows[[k]]
    block <- read_nuclei_block(file, indptr, min(r), max(r), n_genes)
    pieces[[k]] <- pseudobulk_block(block, as.character(meta$subclass[r]), subclasses)
    rm(block); invisible(gc())
    report_progress(k, length(rows), start_time)
  }
  h5closeAll()
  counts <- join_columns(pieces)

  description <- data.frame(donor    = rep(names(rows), each = length(subclasses)),
                            subclass = rep(subclasses, times = length(rows)))
  description$column   <- paste(description$donor, description$subclass, sep = "__")
  description$n_nuclei <- as.vector(vapply(rows, function(r)
    as.integer(table(factor(as.character(meta$subclass[r]), levels = subclasses))),
    integer(length(subclasses))))
  position <- match(description$donor, donors$donor_id)
  description$group    <- as.character(donors$group[position])
  description$age      <- as.numeric(donors$age[position])
  description$sex      <- as.character(donors$sex[position])
  description$ancestry <- as.character(donors$genetic_ancestry[position])

  rownames(counts) <- gene_names
  colnames(counts) <- description$column
  list(counts = counts, description = description)
}
