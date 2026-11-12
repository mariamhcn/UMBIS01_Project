# Helpers for Q5: transcriptomic ageing scores per cell type and donor.
#
# Contents (all the functions of Q5 are in this file):
#  PART 1 – CALCULATIONS: input, scores, age-disease concordance (Part B), disease models
#  PART 2 – FIGURES of Q5 (F03, F04 and the reasoning diagram); the common style (theme, colours, saving) is in
#           utils_figures.R

# Input: the pseudobulk matrix of the preprocessing (counts summed per donor x subclass, all nuclei).
# Column names written by the preprocessing (scripts/04_pseudobulk.R) and the names used in this notebook:
PB_COLUMNS <- c(column = "column", donor_id = "donor", cell_type = "subclass", n_nuclei = "n_nuclei",
                diagnosis = "group", age = "age", sex = "sex", ancestry = "ancestry")

# ===========================================================================
# PART 1 – CALCULATIONS
# ===========================================================================

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

# Composite ageing score of a cell type: mean of the signatures validated in controls in that cell type,
# each one signed so that a higher value always means "older-looking" (age_direction = +1 or -1),
# then re-standardised against controls. `exclude` (optional) lists donor x signature scores to leave out.
add_composite_score <- function(donor_scores, valid_pairs, exclude = NULL, name = "Composite_ageing_score") {
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
# with a two-sided permutation test: ages are shuffled among controls to see which correlations arise by chance
# (a negative correlation, i.e. a younger-looking profile, would also be reported).
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
             p_perm = (sum(abs(null) >= abs(rho)) + 1) / (n_perm + 1), row.names = NULL)
}

# Overlap between the genes significantly changed with age (controls) and with the disease (Q2), with:
#  - a one-sided Fisher exact test on the 2 x 2 table (age gene yes/no x disease gene yes/no): more overlap than chance?
#    (an enrichment test: less overlap than chance has no biological meaning here, as for the ORA of Q3)
#  - a two-sided binomial test on the overlapping genes: does the share changing in the same direction with age and
#    disease differ from one half?
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
    stats::binom.test(n_same, length(both), p = 0.5, alternative = "two.sided")$p.value
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

# 95% confidence interval of a ratio num / den of two estimates (Fieller's method), from their variances and
# covariance. If den (here the effect of age) is not significantly different from 0, the interval has no upper
# limit: NA is returned ("not estimable").
fieller_ci <- function(num, den, v_num, v_den, cov_num_den, df, level = 0.95) {
  t2 <- stats::qt(1 - (1 - level) / 2, df)^2
  A <- den^2 - t2 * v_den
  B <- num * den - t2 * cov_num_den
  C <- num^2 - t2 * v_num
  if (A <= 0) return(c(NA_real_, NA_real_))
  sort(unname((B + c(-1, 1) * sqrt(B^2 - A * C)) / A))
}

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
    # with a 95% confidence interval by Fieller's method (ratio of two regression coefficients)
    b <- stats::coef(fit_off); v <- stats::vcov(fit_off)
    terms <- c("diagnosisSCZ", "diagnosisBD")
    a <- b["age"]
    # numerator weights: SCZ, BD, SCZ - BD
    weights <- list(SCZ = c(1, 0), BD = c(0, 1), SCZ_vs_BD = c(1, -1))
    res$age_slope <- unname(a)
    res$years <- NA_real_; res$years_lower <- NA_real_; res$years_upper <- NA_real_
    for (cmp in names(weights)) {
      w   <- weights[[cmp]]
      num <- sum(w * b[terms])
      ci  <- fieller_ci(num, a, v_num = as.numeric(t(w) %*% v[terms, terms] %*% w), v_den = v["age", "age"],
                        cov_num_den = sum(w * v[terms, "age"]), df = fit_off$df.residual)
      k <- res$model == "offset" & res$comparison == cmp
      res$years[k]       <- unname(num / a)
      res$years_lower[k] <- ci[1]
      res$years_upper[k] <- ci[2]
    }
    cbind(cell_type = groups$cell_type[i], signature = groups$signature[i],
          n_CTRL = sum(d$diagnosis == "CTRL"), n_SCZ = sum(d$diagnosis == "SCZ"),
          n_BD = sum(d$diagnosis == "BD"), res)
  })
  do.call(rbind, out)
}

# Benjamini-Hochberg correction within each family of tests:
# model (offset / slope) x level (composite ageing score = main test / single signatures = detail) x comparison type
# (vs controls / direct SCZ vs BD).
adjust_by_family <- function(models, index_name = "Composite_ageing_score") {
  models |>
    dplyr::mutate(level = ifelse(signature == index_name, "index", "signatures"),
                  comparison_type = ifelse(comparison == "SCZ_vs_BD", "direct", "vs_CTRL")) |>
    dplyr::group_by(model, level, comparison_type) |>
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup()
}

# ===========================================================================
# PART 2 – FIGURES
# ===========================================================================

# ---------------------------------------------------------------------------
# Figure: composite ageing score of patients and controls at equal age
# ---------------------------------------------------------------------------

# Boxplots of the composite ageing score adjusted for age, sex and ancestry, one panel per family (lettered A, B, C)
# and one facet per cell type. adjusted = score - score expected for a control of the same age, sex and ancestry
# (linear model fitted on controls only). Brackets: p_adj of the offset model (models), one per pair of groups.
plot_composite_boxplots <- function(model_data, models, cell_order, cols, composite = "Composite_ageing_score",
                                    covariates = c("sex", "ancestry"), seed = 42) {
  form <- stats::as.formula(paste("score_z ~ age +", paste(covariates, collapse = " + ")))
  grp  <- c(CTRL = "Controls", SCZ = "SCZ", BD = "BD")
  adj <- model_data |>
    dplyr::filter(signature == composite) |>
    dplyr::group_by(cell_type) |>
    dplyr::group_modify(function(d, k) {
      fit <- stats::lm(form, data = d[d$diagnosis == "CTRL", ])
      d$adjusted <- d$score_z - stats::predict(fit, newdata = d)
      d
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(group = factor(grp[as.character(diagnosis)], levels = grp), family = cell_family(cell_type))
  stars <- function(p) ifelse(p < 1e-4, "****", ifelse(p < 1e-3, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
  top <- adj |> dplyr::group_by(cell_type) |> dplyr::summarise(top = max(adjusted, na.rm = TRUE))
  brackets <- models |>
    dplyr::filter(signature == composite, model == "offset") |>
    dplyr::mutate(comparison = as.character(comparison), label = stars(p_adj),
                  g1 = c(SCZ = "Controls", BD = "Controls", SCZ_vs_BD = "SCZ")[comparison],
                  g2 = c(SCZ = "SCZ", BD = "BD", SCZ_vs_BD = "BD")[comparison],
                  level = c(SCZ = 1, SCZ_vs_BD = 2, BD = 3)[comparison]) |>
    dplyr::inner_join(top, by = "cell_type") |>
    dplyr::mutate(y = top + 0.4 + 0.6 * (level - 1), x1 = match(g1, grp), x2 = match(g2, grp),
                  family = cell_family(cell_type))
  cell_order <- cell_order[cell_order %in% adj$cell_type]
  one_family <- function(fam) {
    cts <- cell_order[cell_family(cell_order) == fam]
    d <- dplyr::filter(adj, family == fam) |> dplyr::mutate(cell_type = factor(cell_type, levels = cts))
    b <- dplyr::filter(brackets, family == fam) |> dplyr::mutate(cell_type = factor(cell_type, levels = cts))
    ggplot2::ggplot(d, ggplot2::aes(group, adjusted, colour = group)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15, height = 0, seed = seed),
                          size = 1, alpha = 0.5, show.legend = FALSE) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = group), alpha = 0.15, outlier.shape = NA, width = 0.6,
                            linewidth = 0.6, key_glyph = "rect") +
      ggplot2::geom_segment(data = b, inherit.aes = FALSE, ggplot2::aes(x = x1, xend = x2, y = y, yend = y), linewidth = 0.35) +
      ggplot2::geom_segment(data = b, inherit.aes = FALSE, ggplot2::aes(x = x1, xend = x1, y = y, yend = y - 0.15), linewidth = 0.35) +
      ggplot2::geom_segment(data = b, inherit.aes = FALSE, ggplot2::aes(x = x2, xend = x2, y = y, yend = y - 0.15), linewidth = 0.35) +
      ggplot2::geom_text(data = b, inherit.aes = FALSE, size = 3.2,
                         ggplot2::aes(x = (x1 + x2) / 2, y = y, label = label, vjust = ifelse(label == "ns", -0.25, 0.3))) +
      ggplot2::facet_wrap(~ cell_type, nrow = 1, scales = "free_y") +
      ggplot2::scale_colour_manual(values = cols, name = NULL, guide = "none") +
      ggplot2::scale_fill_manual(values = cols, name = NULL) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1))) +
      ggplot2::labs(x = NULL, y = "Adjusted composite\nageing score (SD)", title = fam) +
      theme_project(11) +
      ggplot2::theme(plot.title = ggplot2::element_text(colour = "black", face = "bold", size = 14, hjust = 0.5),
                     plot.title.position = "panel",
                     strip.text = ggplot2::element_text(size = 12, face = "bold", colour = "black"),
                     strip.background = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank())
  }
  families <- intersect(c("Excitatory neurons", "Inhibitory neurons", "Glial cells", "Vascular cells"),
                        cell_family(cell_order))
  patchwork::wrap_plots(lapply(families, one_family), ncol = 1) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      tag_levels = "A", theme = theme_project(),
      title = "At equal age, patients have a higher composite ageing score than controls in most cell types",
      subtitle = "Composite ageing score adjusted for age, sex and ancestry: 0 = the score expected for a control of the same age",
      caption = paste0("One dot per donor; above 0 = older-looking than expected for the age. Only the cell types with a composite score are shown (at least one signature following age in controls).\n",
                       "Brackets: difference between two groups at equal age (linear model adjusted for age, sex and ancestry), Benjamini-Hochberg adjusted (ns: not significant; * < 0.05; ** < 0.01; *** < 0.001; **** < 0.0001).")) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16, hjust = 0, vjust = 1),
                   plot.tag.position = c(0, 1), legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 12))
}

# ---------------------------------------------------------------------------
# Figure: shift of each score of the composite ageing score
# ---------------------------------------------------------------------------

# Heatmap with the same layout as Figure F01 (cell types grouped by family, a bracket naming each family), restricted
# to the cell types with a composite score and to the scores kept in it (validated in controls). One panel per
# comparison with controls (SCZ, BD); colour = shift towards an older profile, stars = adjusted p of that score.
plot_score_shifts <- function(models, composite = "Composite_ageing_score", alpha = 0.05) {
  stars <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 0.01, "**", ifelse(p < alpha, "*", "")))
  d <- models |>
    dplyr::filter(signature != composite, model == "offset", comparison %in% c("SCZ", "BD"), validated)
  lev <- rev(order_cell_types(d$cell_type))
  d <- d |>
    dplyr::mutate(cell_type = factor(cell_type, levels = lev), signature = factor(signature),
                  comparison = factor(as.character(comparison), levels = c("SCZ", "BD"),
                                      labels = c("SCZ vs controls", "BD vs controls")))
  fam <- data.frame(y = seq_along(lev), family = cell_family(lev)) |>
    dplyr::group_by(family) |>
    dplyr::summarise(ymin = min(y) - 0.4, ymax = max(y) + 0.4, ymid = mean(y), .groups = "drop") |>
    dplyr::mutate(label = sub(" ", "\n", family))
  lim <- max(abs(d$aging_shift), na.rm = TRUE)
  heat <- ggplot2::ggplot(d, ggplot2::aes(signature, cell_type, fill = aging_shift)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = stars(p_adj)), size = 5, vjust = 0.75) +
    ggplot2::facet_wrap(~ comparison) +
    ggplot2::scale_fill_gradient2(low = "#5e3c99", mid = "white", high = "#e66101", midpoint = 0,
                                  limits = c(-lim, lim), name = "Shift towards\nolder profile\n(control SD)") +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(x = "Signature", y = NULL,
                  title = "The older-looking profile is carried by the lifespan and microglial ageing genes",
                  subtitle = paste0("Shift of each score kept in the composite ageing score, in patients compared with controls of the same age; * adjusted p < ", alpha),
                  caption = "Orange: older-looking than controls of the same age. Grey: not a validated ageing score in that cell type. No score differs between SCZ and BD (direct test).") +
    theme_project() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1), panel.grid = ggplot2::element_blank(),
                   legend.position = "right", panel.background = ggplot2::element_rect(fill = "grey90"),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold", size = 12))
  # brackets naming the families, drawn in a narrow plot on the left with the same vertical scale
  brackets <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0.5, length(lev) + 0.5), expand = c(0, 0)) +
    ggplot2::annotate("segment", x = 0.95, xend = 0.95, y = fam$ymin, yend = fam$ymax) +
    ggplot2::annotate("segment", x = 0.95, xend = 1, y = fam$ymin, yend = fam$ymin) +
    ggplot2::annotate("segment", x = 0.95, xend = 1, y = fam$ymax, yend = fam$ymax) +
    ggplot2::annotate("text", x = 0.9, y = fam$ymid, label = fam$label, hjust = 1, fontface = "bold", size = 3.8,
                      lineheight = 0.9) +
    ggplot2::theme_void()
  brackets + heat + patchwork::plot_layout(widths = c(0.12, 1))
}

# ---------------------------------------------------------------------------
# Figure: how the notebook answers Q5 (reasoning diagram, no data)
# ---------------------------------------------------------------------------

# Content of the diagram drawn by plot_reasoning_diagram() (utils_figures.R): one column per stage, one box per step.
plot_reasoning <- function() {
  cols <- data.frame(head = c("1. Build the ruler", "2. Measure the patients", "3. Check with all genes"),
                     sub = c("Part A, steps A1 to A4", "Part A, steps A5 to A7", "Part B"),
                     fill = c("#fdf3c4", "#fbe38c", "#f5d35a"), border = "#a88600")
  boxes <- data.frame(
    col  = c(1, 1, 1, 1, 2, 2, 2, 3, 3),
    step = c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "B1", "B2"),
    text = c("Published gene lists,\none per ageing hallmark",
             "One score per donor\nand cell type",
             "Which scores really\nfollow age in\nhealthy controls?",
             "Combine them: one\ncomposite ageing score\nper cell type",
             "Do patients look older\nthan controls of the\nsame age? By how much?",
             "Which hallmark\ndrives the shift?",
             "A head start,\nor faster ageing?",
             "With all genes: does the\ndisease move genes\nlike age does?",
             "Are they the same\nsignificant genes?"),
    out = c("", "", "F01", "", "F02 + table", "F03", "table", "F04", "table"))
  plot_reasoning_diagram(cols, boxes,
    conclusion = "Conclusion: answer to Q5, H3, and comparison with the literature",
    title = "How the notebook answers Q5",
    subtitle = "We first build an ageing ruler in healthy controls, then measure the patients with it, then check the result without chosen gene lists")
}
