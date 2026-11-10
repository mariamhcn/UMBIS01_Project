# Helpers for Q5: transcriptomic ageing scores per cell type and donor.

# Input: the pseudobulk matrix of the preprocessing (counts summed per donor x subclass, all nuclei).
# Column names written by the preprocessing (scripts/04_pseudobulk.R) and the names used in this notebook:
PB_COLUMNS <- c(column = "column", donor_id = "donor", cell_type = "subclass", n_nuclei = "n_nuclei",
                diagnosis = "group", age = "age", sex = "sex", ancestry = "ancestry")

# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

# Reads the pseudobulk file (a list: `counts` = genes x columns, raw counts; `description` = one row per column),
# checks the columns and renames them (donor -> donor_id, subclass -> cell_type, group -> diagnosis, Control -> CTRL).
load_pseudobulk <- function(path, control_label = "Control") {
  pb <- readRDS(path)
  desc <- pb$description
  missing_cols <- setdiff(PB_COLUMNS, colnames(desc))
  if (length(missing_cols) > 0) stop("Missing columns in the pseudobulk description: ", paste(missing_cols, collapse = ", "))
  desc <- desc[, PB_COLUMNS]
  colnames(desc) <- names(PB_COLUMNS)
  desc$diagnosis <- factor(ifelse(desc$diagnosis == control_label, "CTRL", desc$diagnosis), levels = c("CTRL", "SCZ", "BD"))
  if (anyNA(desc$diagnosis)) stop("group must be Control, SCZ or BD")
  stopifnot(identical(colnames(pb$counts), desc$column))
  list(counts = pb$counts, desc = desc)
}

# Cell types analysed: a cell type is kept if at least `min_donors` donors of EACH group have a profile
# built from at least `min_nuclei` nuclei (shared rule E7). Returns one row per cell type, for the report.
select_celltypes <- function(desc, cell_types, min_nuclei = 20, min_donors = 10) {
  desc |>
    dplyr::filter(cell_type %in% cell_types, n_nuclei >= min_nuclei) |>
    dplyr::count(cell_type, diagnosis, name = "n_donors") |>
    tidyr::pivot_wider(names_from = diagnosis, values_from = n_donors, values_fill = 0) |>
    dplyr::right_join(data.frame(cell_type = cell_types), by = "cell_type") |>
    dplyr::mutate(dplyr::across(-cell_type, ~ dplyr::coalesce(.x, 0L)),
                  kept = CTRL >= min_donors & SCZ >= min_donors & BD >= min_donors)
}

# ---------------------------------------------------------------------------
# Scores
# ---------------------------------------------------------------------------

# Scores the pseudobulk profiles of one cell type (one column per donor) with AddModuleScore (method seen in the TP):
# the profiles are log-normalised (NormalizeData), then each signature is compared with control genes of similar
# expression level, so the score does not depend on sequencing depth. Signatures with too few genes are skipped.
# Returns one row per donor and signature.
score_pseudobulk_celltype <- function(counts_ct, donor_ids, signatures, seed = 42, min_genes = 5) {
  colnames(counts_ct) <- donor_ids
  counts_ct <- counts_ct[Matrix::rowSums(counts_ct) > 0, , drop = FALSE]   # genes not detected in this cell type
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts_ct))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  signatures <- lapply(signatures, intersect, rownames(obj))
  signatures <- signatures[vapply(signatures, length, integer(1)) >= min_genes]
  if (length(signatures) == 0) return(NULL)
  obj <- Seurat::AddModuleScore(obj, features = signatures, name = "aging_sig_", seed = seed)
  scores <- obj@meta.data[, paste0("aging_sig_", seq_along(signatures)), drop = FALSE]
  colnames(scores) <- names(signatures)
  data.frame(donor_id = rownames(scores), scores, check.names = FALSE) |>
    tidyr::pivot_longer(-donor_id, names_to = "signature", values_to = "score")
}

# log2-CPM of pseudobulk profiles (genes x donors), keeping genes with a mean CPM >= min_cpm.
pseudobulk_logcpm <- function(counts, min_cpm = 1) {
  pb <- as.matrix(counts)
  cpm <- t(t(pb) / colSums(pb)) * 1e6
  keep <- rowMeans(cpm) >= min_cpm
  log2(cpm[keep, , drop = FALSE] + 1)
}

# Express each donor score as a z-score relative to controls of the same cell type,
# so effects are comparable across signatures and cell types (1 unit = 1 control SD).
standardize_by_controls <- function(donor_scores, control_label = "CTRL") {
  donor_scores |>
    dplyr::group_by(cell_type, signature) |>
    dplyr::mutate(
      ctrl_mean = mean(score[diagnosis == control_label], na.rm = TRUE),
      ctrl_sd   = stats::sd(score[diagnosis == control_label], na.rm = TRUE),
      score_z   = (score - ctrl_mean) / ctrl_sd
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-ctrl_mean, -ctrl_sd)
}

# Does each score change with age in controls? Spearman correlation per cell type x signature.
validate_in_controls <- function(donor_scores, control_label = "CTRL", min_donors = 10) {
  donor_scores |>
    dplyr::filter(diagnosis == control_label, !is.na(score_z)) |>
    dplyr::group_by(cell_type, signature) |>
    dplyr::filter(dplyr::n() >= min_donors) |>
    dplyr::summarise(
      n_donors = dplyr::n(),
      rho      = suppressWarnings(stats::cor(age, score_z, method = "spearman")),
      p_value  = suppressWarnings(stats::cor.test(age, score_z, method = "spearman", exact = FALSE)$p.value),
      .groups  = "drop"
    ) |>
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH"))
}

# Ageing index of a cell type: mean of the signatures validated in controls in that cell type,
# each one signed so that a higher value always means "older-looking" (age_direction = +1 or -1),
# then re-standardised against controls. `exclude` (optional) lists donor x signature scores to leave out.
add_aging_index <- function(donor_scores, valid_pairs, exclude = NULL, name = "Ageing_index") {
  d <- donor_scores |>
    dplyr::inner_join(valid_pairs, by = c("cell_type", "signature"))
  if (!is.null(exclude)) d <- dplyr::anti_join(d, exclude, by = c("donor_id", "signature"))
  idx <- d |>
    dplyr::group_by(cell_type, donor_id) |>
    dplyr::summarise(score        = mean(score_z * age_direction, na.rm = TRUE),
                     n_nuclei     = min(n_nuclei),
                     n_signatures = dplyr::n(),
                     .groups = "drop") |>
    dplyr::mutate(signature = name)
  meta_cols <- setdiff(colnames(donor_scores), c(colnames(idx), "score_z"))
  meta <- dplyr::distinct(donor_scores[, c("donor_id", meta_cols)])
  idx <- dplyr::left_join(idx, meta, by = "donor_id") |> standardize_by_controls()
  dplyr::bind_rows(donor_scores, idx)
}

# ---------------------------------------------------------------------------
# Age-disease concordance (all genes)
# ---------------------------------------------------------------------------

# Effect of age on every gene in controls, for one cell type:
# log-CPM of the pseudobulk profiles of the control donors, then a linear model gene ~ age + covariates, fitted for all
# genes at once. Returns the t statistic of age per gene (positive = rises with age), and the data needed for permutations.
age_effects_in_controls <- function(counts_ct, donor_ids, donor_meta, covariates = c("sex", "ancestry"),
                                    control_label = "CTRL", min_donors = 10) {
  ctrl <- donor_ids %in% donor_meta$donor_id[donor_meta$diagnosis == control_label]
  if (sum(ctrl) < min_donors) return(NULL)
  counts <- counts_ct[, ctrl, drop = FALSE]
  colnames(counts) <- donor_ids[ctrl]
  logcpm <- pseudobulk_logcpm(counts)
  meta <- donor_meta[match(colnames(logcpm), donor_meta$donor_id), ]
  keep <- stats::complete.cases(meta[, c("age", covariates)])
  meta <- meta[keep, ]; logcpm <- logcpm[, keep, drop = FALSE]
  X <- stats::model.matrix(stats::as.formula(paste("~ age", if (length(covariates)) paste("+", paste(covariates, collapse = " + ")) else "")), meta)
  list(t_age = age_t_statistic(t(logcpm), X), Y = t(logcpm), X = X)
}

# t statistic of the age coefficient for every column (gene) of Y, with design matrix X.
age_t_statistic <- function(Y, X) {
  fit <- stats::lm.fit(X, Y)
  df <- nrow(X) - ncol(X)
  sigma2 <- colSums(fit$residuals^2) / df
  xtx_inv <- chol2inv(qr.R(fit$qr))[which(colnames(X) == "age"), which(colnames(X) == "age")]
  stats::setNames(fit$coefficients["age", ] / sqrt(sigma2 * xtx_inv), colnames(Y))
}

# Spearman correlation between the age effect (controls) and the disease effect (Q2 Wald statistic) over all genes,
# with a permutation test: ages are shuffled among controls to see which correlations arise by chance.
age_disease_concordance <- function(age_fit, disease_stat, n_perm = 200, seed = 42) {
  genes <- intersect(names(age_fit$t_age), names(disease_stat))
  if (length(genes) < 100) return(NULL)
  rho <- stats::cor(age_fit$t_age[genes], disease_stat[genes], method = "spearman", use = "complete.obs")
  set.seed(seed)
  Y <- age_fit$Y[, genes, drop = FALSE]
  null <- replicate(n_perm, {
    Xp <- age_fit$X
    Xp[, "age"] <- sample(Xp[, "age"])
    stats::cor(age_t_statistic(Y, Xp), disease_stat[genes], method = "spearman", use = "complete.obs")
  })
  data.frame(n_genes = length(genes), rho = rho,
             null_low = stats::quantile(null, 0.025), null_high = stats::quantile(null, 0.975),
             p_perm = (sum(null >= rho) + 1) / (n_perm + 1), row.names = NULL)
}

# Overlap between the genes significantly changed with age (controls) and with the disease (Q2), with:
#  - a one-sided Fisher exact test on the 2 x 2 table (age gene yes/no x disease gene yes/no): more overlap than chance?
#  - a binomial test on the overlapping genes: do more than half change in the same direction with age and disease?
# Universe: genes tested in both analyses (Q2 padj not NA). Age p-values come from the t statistic of age_t_statistic().
age_disease_overlap <- function(age_fit, disease_stat, disease_padj, alpha = 0.05, min_overlap = 5) {
  df_age <- nrow(age_fit$X) - ncol(age_fit$X)
  p_age  <- 2 * stats::pt(abs(age_fit$t_age), df = df_age, lower.tail = FALSE)
  genes  <- intersect(names(p_age)[!is.na(p_age)],
                      intersect(names(disease_padj)[!is.na(disease_padj)], names(disease_stat)[!is.na(disease_stat)]))
  if (length(genes) < 100) return(NULL)
  age_sig <- stats::p.adjust(p_age[genes], method = "BH") < alpha
  dis_sig <- disease_padj[genes] < alpha
  tab <- table(age = factor(age_sig, levels = c(TRUE, FALSE)), disease = factor(dis_sig, levels = c(TRUE, FALSE)))
  ft  <- stats::fisher.test(tab, alternative = "greater")
  both   <- genes[age_sig & dis_sig]
  n_same <- sum(sign(age_fit$t_age[both]) == sign(disease_stat[both]))
  p_dir  <- if (length(both) >= min_overlap) {
    stats::binom.test(n_same, length(both), p = 0.5, alternative = "greater")$p.value
  } else NA_real_
  data.frame(n_genes = length(genes), n_age = sum(age_sig), n_disease = sum(dis_sig), n_both = length(both),
             expected_both = round(sum(age_sig) * sum(dis_sig) / length(genes), 1),
             odds_ratio = unname(ft$estimate), p_fisher = ft$p.value,
             pct_same_direction = if (length(both) > 0) round(100 * n_same / length(both), 1) else NA_real_,
             p_direction = p_dir, row.names = NULL)
}

# ---------------------------------------------------------------------------
# Disease models
# ---------------------------------------------------------------------------

# Difference between two coefficients of an lm fit (e.g. SCZ effect - BD effect), with its standard error.
contrast_coefs <- function(fit, term_a, term_b) {
  b <- stats::coef(fit); v <- stats::vcov(fit)
  if (!all(c(term_a, term_b) %in% names(b))) return(NULL)
  est <- unname(b[term_a] - b[term_b])
  se  <- sqrt(v[term_a, term_a] + v[term_b, term_b] - 2 * v[term_a, term_b])
  p   <- 2 * stats::pt(abs(est / se), df = fit$df.residual, lower.tail = FALSE)
  data.frame(term = paste(term_a, "-", term_b), estimate = est, std_error = se, p_value = p)
}

tidy_coefs <- function(fit, terms) {
  co <- summary(fit)$coefficients
  terms <- intersect(terms, rownames(co))
  data.frame(term = terms, estimate = co[terms, "Estimate"], std_error = co[terms, "Std. Error"],
             p_value = co[terms, "Pr(>|t|)"], row.names = NULL)
}

# Disease effect at equal age, per cell type x signature:
#  - offset model:      score_z ~ age + diagnosis + covariates  -> patients look "older" than controls of the same age
#  - interaction model: score_z ~ age * diagnosis + covariates  -> the score changes faster with age in patients
# For each model: SCZ vs CTRL, BD vs CTRL, and the direct comparison SCZ vs BD (difference of the two effects).
fit_disease_models <- function(donor_scores, covariates = c("sex", "ancestry"), min_donors_per_group = 10) {
  groups <- unique(donor_scores[, c("cell_type", "signature")])
  out <- lapply(seq_len(nrow(groups)), function(i) {
    d <- donor_scores[donor_scores$cell_type == groups$cell_type[i] &
                        donor_scores$signature == groups$signature[i], ]
    d <- d[stats::complete.cases(d[, c("score_z", "age", "diagnosis", covariates)]), ]
    d$diagnosis <- factor(d$diagnosis, levels = c("CTRL", "SCZ", "BD"))
    if (any(table(d$diagnosis) < min_donors_per_group)) return(NULL)
    cov_part <- if (length(covariates) > 0) paste("+", paste(covariates, collapse = " + ")) else ""
    fit_off <- stats::lm(stats::as.formula(paste("score_z ~ age + diagnosis", cov_part)), data = d)
    fit_int <- stats::lm(stats::as.formula(paste("score_z ~ age * diagnosis", cov_part)), data = d)
    res <- rbind(
      cbind(model = "offset", comparison = c("SCZ", "BD"),
            tidy_coefs(fit_off, c("diagnosisSCZ", "diagnosisBD"))),
      cbind(model = "offset", comparison = "SCZ_vs_BD",
            contrast_coefs(fit_off, "diagnosisSCZ", "diagnosisBD")),
      cbind(model = "slope", comparison = c("SCZ", "BD"),
            tidy_coefs(fit_int, c("age:diagnosisSCZ", "age:diagnosisBD"))),
      cbind(model = "slope", comparison = "SCZ_vs_BD",
            contrast_coefs(fit_int, "age:diagnosisSCZ", "age:diagnosisBD"))
    )
    # apparent age gap (offset model): disease effect / age effect = how many years older-looking,
    # with a delta-method standard error
    b <- stats::coef(fit_off); v <- stats::vcov(fit_off)
    terms <- c("diagnosisSCZ", "diagnosisBD", "age")
    a <- b["age"]
    # numerator weights: SCZ, BD, SCZ - BD
    weights <- list(SCZ = c(1, 0), BD = c(0, 1), SCZ_vs_BD = c(1, -1))
    res$age_slope <- unname(a)
    res$years <- NA_real_; res$years_se <- NA_real_
    for (cmp in names(weights)) {
      w <- weights[[cmp]]
      num <- sum(w * b[terms[1:2]])
      grad <- c(w / a, -num / a^2)                       # derivative of num / a
      k <- res$model == "offset" & res$comparison == cmp
      res$years[k]    <- unname(num / a)
      res$years_se[k] <- sqrt(as.numeric(t(grad) %*% v[terms, terms] %*% grad))
    }
    cbind(cell_type = groups$cell_type[i], signature = groups$signature[i],
          n_CTRL = sum(d$diagnosis == "CTRL"), n_SCZ = sum(d$diagnosis == "SCZ"),
          n_BD = sum(d$diagnosis == "BD"), res)
  })
  do.call(rbind, out)
}

# Benjamini-Hochberg correction within each family of tests:
# model (offset / slope) x level (ageing index = main test / single signatures = detail) x comparison type
# (vs controls / direct SCZ vs BD).
adjust_by_family <- function(models, index_name = "Ageing_index") {
  models |>
    dplyr::mutate(level = ifelse(signature == index_name, "index", "signatures"),
                  comparison_type = ifelse(comparison == "SCZ_vs_BD", "direct", "vs_CTRL")) |>
    dplyr::group_by(model, level, comparison_type) |>
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup()
}
