# Helpers for Q5: transcriptomic ageing scores per cell type and donor.

META_REQUIRED_COLUMNS <- c("donor_id", "diagnosis", "age", "sex", "pmi")

# ---------------------------------------------------------------------------
# Analysis units (shared rule of the project, see README.md at the project root)
# ---------------------------------------------------------------------------

# A subclass is analysed on its own if at least `min_donors` donors of EACH group have at least
# `min_cells` nuclei of it. Otherwise its nuclei are merged into their class.
# Returns a list: `unit` (one value per nucleus) and `decision` (one row per subclass, for the report).
assign_analysis_units <- function(meta, subclass_col = "cell_subclass", class_col = "cell_class",
                                  min_cells = 20, min_donors = 10) {
  counts <- meta |>
    dplyr::count(subclass = .data[[subclass_col]], class = .data[[class_col]],
                 diagnosis, donor_id, name = "n_nuclei") |>
    dplyr::filter(n_nuclei >= min_cells) |>
    dplyr::count(subclass, class, diagnosis, name = "n_donors") |>
    tidyr::pivot_wider(names_from = diagnosis, values_from = n_donors, values_fill = 0)
  groups <- intersect(c("CTRL", "SCZ", "BD"), colnames(counts))
  all_sub <- dplyr::distinct(meta, subclass = .data[[subclass_col]], class = .data[[class_col]])
  decision <- dplyr::left_join(all_sub, counts, by = c("subclass", "class")) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(groups), ~ dplyr::coalesce(.x, 0L)))
  decision$enough_donors <- apply(decision[, groups, drop = FALSE] >= min_donors, 1, all)
  # merged subclasses are pooled into "<class>_other" when other subclasses of the class are kept,
  # and into "<class>" when none is kept
  class_has_kept <- tapply(decision$enough_donors, decision$class, any)
  decision$unit <- ifelse(decision$enough_donors, decision$subclass,
                          ifelse(class_has_kept[decision$class], paste0(decision$class, "_other"), decision$class))
  unit <- decision$unit[match(meta[[subclass_col]], decision$subclass)]
  list(unit = unit, decision = decision)
}

# ---------------------------------------------------------------------------
# Scores
# ---------------------------------------------------------------------------

# Score every nucleus of one cell type with AddModuleScore (method seen in the TP).
# Each signature is compared with control genes of similar expression level,
# so the score is not driven by sequencing depth. Signatures with too few genes in the data are skipped.
score_one_celltype <- function(sobj_ct, signatures, seed = 42, min_genes = 5) {
  genes_in_data <- rownames(sobj_ct)
  signatures <- lapply(signatures, intersect, genes_in_data)
  signatures <- signatures[vapply(signatures, length, integer(1)) >= min_genes]
  if (length(signatures) == 0) return(NULL)
  sobj_ct <- Seurat::AddModuleScore(sobj_ct, features = signatures, name = "aging_sig_", seed = seed)
  scores <- sobj_ct@meta.data[, paste0("aging_sig_", seq_along(signatures)), drop = FALSE]
  colnames(scores) <- names(signatures)
  scores
}

# Pseudobulk log-CPM per donor (sum of raw counts of all nuclei of a donor).
pseudobulk_logcpm <- function(counts, donor_ids, min_cpm = 1) {
  design <- Matrix::sparse.model.matrix(~ 0 + donor, data.frame(donor = factor(donor_ids)))
  colnames(design) <- sub("^donor", "", colnames(design))
  pb <- as.matrix(counts %*% design)
  cpm <- t(t(pb) / colSums(pb)) * 1e6
  keep <- rowMeans(cpm) >= min_cpm
  log2(cpm[keep, , drop = FALSE] + 1)
}

# Mean score per donor, keeping only donors with enough nuclei of that cell type
# (a mean over a handful of nuclei is too noisy).
aggregate_by_donor <- function(cell_scores, donor_ids, cell_type, min_cells = 20) {
  long <- data.frame(donor_id = donor_ids, cell_scores, check.names = FALSE) |>
    tidyr::pivot_longer(-donor_id, names_to = "signature", values_to = "score")
  long |>
    dplyr::group_by(donor_id, signature) |>
    dplyr::summarise(score = mean(score), n_cells = dplyr::n(), .groups = "drop") |>
    dplyr::filter(n_cells >= min_cells) |>
    dplyr::mutate(cell_type = cell_type, .before = 1)
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
                     n_cells      = min(n_cells),
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
# pseudobulk log-CPM per control donor, then a linear model gene ~ age + covariates, fitted for all genes at once.
# Returns the t statistic of age per gene (positive = rises with age), and the data needed for permutations.
age_effects_in_controls <- function(sobj_ct, donor_meta, covariates = c("sex", "pmi"),
                                    control_label = "CTRL", min_donors = 10) {
  ctrl <- donor_meta$donor_id[donor_meta$diagnosis == control_label]
  cells <- colnames(sobj_ct)[sobj_ct$donor_id %in% ctrl]
  if (length(unique(sobj_ct$donor_id[cells])) < min_donors) return(NULL)
  counts <- SeuratObject::LayerData(sobj_ct, assay = "RNA", layer = "counts")[, cells, drop = FALSE]
  logcpm <- pseudobulk_logcpm(counts, sobj_ct$donor_id[cells])
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
fit_disease_models <- function(donor_scores, covariates = c("sex", "pmi"), min_donors_per_group = 10) {
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
