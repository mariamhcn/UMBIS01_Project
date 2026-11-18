# Functions for the cell-type composition test (Q1), with propeller (speckle package).
# Used by notebook 02_Q1_composition.qmd, for the main analysis and the sensitivity analyses.
#
# Author of the method: Myriam (script 03 of Q1).

suppressPackageStartupMessages({
  library(speckle)
  library(limma)
})

# Contrast weights: -1 for the reference group, +1 for the tested group, 0 for every other column.
make_contrast <- function(design, reference_group, tested_group) {
  weights <- rep(0, ncol(design))
  weights[colnames(design) == paste0("group", reference_group)] <- -1
  weights[colnames(design) == paste0("group", tested_group)]    <-  1
  weights
}

# propeller test on counts (donors x subclasses) for a list of comparisons c(reference, tested).
# Model: logit(proportion) ~ 0 + group + age + sex + ancestry (ancestry only if it varies).
# Returns one row per subclass and comparison: t, p-value, FDR and the raw ratio of mean proportions.
#
# Reading the output of propeller.ttest: Tstatistic, P.Value and FDR come from the ADJUSTED model;
# PropRatio is NOT adjusted (ratio of the raw mean proportions of the two groups).
test_composition <- function(counts, info, comparisons) {
  info <- droplevels(info)
  props <- convertDataToList(t(counts), data.type = "counts", transform = "logit")
  if (nlevels(info$ancestry) > 1) {
    design <- model.matrix(~ 0 + group + age + sex + ancestry, data = info)
  } else {
    design <- model.matrix(~ 0 + group + age + sex, data = info)
  }
  pieces <- lapply(comparisons, function(comparison) {
    res <- propeller.ttest(props, design = design,
                           contrasts = make_contrast(design, comparison[1], comparison[2]),
                           robust = TRUE, trend = FALSE, sort = FALSE)
    data.frame(subclass = rownames(res),
               contrast = paste0(comparison[2], "_vs_", comparison[1]),
               t        = res$Tstatistic,
               p_value  = res$P.Value,
               FDR      = res$FDR,
               ratio    = res$PropRatio)
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}
