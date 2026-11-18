# Functions to select the donors of the HBCC cohort, as patients are screened for a clinical trial.
# Used by notebook 01_preprocessing.qmd. They work on the nucleus metadata only (no expression).
#
# Author of the selection method: Myriam (script 01 of the preprocessing).

# ---------------------------------------------------------------------------------------
# From nuclei to donors
# ---------------------------------------------------------------------------------------

# One row per donor (people are counted, not nuclei), with the number of nuclei of each donor.
donors_from_nuclei <- function(meta) {
  donors <- meta[!duplicated(meta$donor_id),
                 c("donor_id", "disease", "development_stage", "sex",
                   "genetic_ancestry", "ASCVD_status", "Tardive_dyskinesia")]
  donors$n_nuclei <- as.integer(table(meta$donor_id)[as.character(donors$donor_id)])
  donors$donor_id <- as.character(donors$donor_id)
  donors$disease  <- as.character(donors$disease)
  donors$sex      <- as.character(donors$sex)
  rownames(donors) <- NULL
  donors
}

# Age in years from the text of the file: "45-year-old stage" -> 45, "8-month-old stage" -> 0.67.
# A format not recognised gives NA (checked in the notebook).
parse_age <- function(stage) {
  stage  <- as.character(stage)
  number <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", stage)))
  age <- rep(NA_real_, length(stage))
  in_years  <- grepl("year-old", stage)
  in_months <- grepl("month-old", stage)
  age[in_years]  <- number[in_years]
  age[in_months] <- number[in_months] / 12
  age
}

# ---------------------------------------------------------------------------------------
# One TRUE/FALSE column per diagnosis
# ---------------------------------------------------------------------------------------

# Diagnoses known to the selection; whatever remains in the text is "another diagnosis".
# The order matters: long expressions are removed before the short ones they contain.
KNOWN_DIAGNOSES <- c("schizoaffective bipolar disorder", "schizoaffective depressive disorder",
                     "schizophrenia", "bipolar II disorder", "bipolar I disorder", "bipolar disorder",
                     "type 1 diabetes mellitus", "type 2 diabetes mellitus", "diabetes mellitus",
                     "atherosclerosis", "anorexia nervosa", "bulimia nervosa", "tardive dyskinesia",
                     "normal")

# The "disease" column lists all the diagnoses of a donor, separated by " || ".
flag_diagnoses <- function(donors) {
  dx <- donors$disease
  # "Schizoaffective bipolar disorder" contains the word "bipolar": it is flagged first, then
  # removed from the text so that it is not confused with bipolar disorder.
  donors$has_schizoaff_bipolar <- grepl("schizoaffective bipolar disorder", dx)
  dx_without_schizoaff_bip     <- gsub("schizoaffective bipolar disorder", "", dx)

  donors$has_SCZ <- grepl("schizophrenia", dx) | grepl("schizoaffective depressive disorder", dx)   # I5
  donors$has_BD  <- grepl("bipolar", dx_without_schizoaff_bip)                                      # I6
  donors$has_tardive_dyskinesia <- grepl("tardive dyskinesia", dx)
  donors$has_metabolic          <- grepl("diabetes", dx) | grepl("atherosclerosis", dx)
  donors$has_eating_disorder    <- grepl("anorexia nervosa", dx) | grepl("bulimia nervosa", dx)
  donors$is_normal              <- dx == "normal"

  remaining <- dx
  for (known in KNOWN_DIAGNOSES) remaining <- gsub(known, "", remaining, fixed = TRUE)
  remaining <- trimws(gsub("\\|\\|", "", remaining))
  donors$other_diagnosis     <- remaining
  donors$has_other_diagnosis <- remaining != ""
  donors
}

# ---------------------------------------------------------------------------------------
# Inclusion and exclusion criteria
# ---------------------------------------------------------------------------------------

# Returns the group of each donor (Control, SCZ, BD or Excluded) and the reason of exclusion.
# Each exclude(...) line is one criterion; a donor keeps the FIRST reason that excludes it.
# Written as a function so that it can be repeated with other parameters (sensitivity analyses).
#   schizoaff_bip: "exclude" (main analysis), or "SCZ" / "BD" to assign these donors to a group.
select_donors <- function(donors, age_min, keep_td, keep_metabolic_controls,
                          schizoaff_bip = "exclude") {
  has_SCZ <- donors$has_SCZ
  has_BD  <- donors$has_BD
  if (schizoaff_bip == "SCZ") { has_SCZ <- has_SCZ | donors$has_schizoaff_bipolar;  has_BD <- has_BD & !donors$has_schizoaff_bipolar }
  if (schizoaff_bip == "BD")  { has_BD  <- has_BD  | donors$has_schizoaff_bipolar;  has_SCZ <- has_SCZ & !donors$has_schizoaff_bipolar }

  is_patient <- has_SCZ | has_BD
  is_control <- !has_SCZ & !has_BD & !donors$has_schizoaff_bipolar

  reason <- rep(NA_character_, nrow(donors))
  exclude <- function(reason, condition, text) {
    reason[condition & is.na(reason)] <- text            # a reason already given is never overwritten
    reason
  }
  reason <- exclude(reason, is.na(donors$age) | is.na(donors$sex),      "E6 - missing data")
  if (schizoaff_bip == "exclude") {
    reason <- exclude(reason, donors$has_schizoaff_bipolar,             "E2 - schizoaffective disorder, bipolar type")
  }
  reason <- exclude(reason, has_SCZ & has_BD,                           "E1 - SCZ and BD")
  reason <- exclude(reason, is_patient & donors$has_other_diagnosis,    "E4 - other diagnosis (patient)")
  reason <- exclude(reason, is_control & (donors$has_other_diagnosis |
                                          donors$has_eating_disorder |
                                          donors$has_tardive_dyskinesia), "E5 - other diagnosis (control)")
  if (!keep_td) {
    reason <- exclude(reason, is_patient & donors$has_tardive_dyskinesia, "Tardive dyskinesia")
  }
  if (!keep_metabolic_controls) {
    reason <- exclude(reason, is_control & donors$has_metabolic,       "E5 - control with diabetes or atherosclerosis")
  }
  reason <- exclude(reason, donors$age < age_min,                      "I2 - age below threshold")

  group <- rep(NA_character_, nrow(donors))
  group[has_SCZ]        <- "SCZ"
  group[has_BD]         <- "BD"
  group[is_control]     <- "Control"
  group[!is.na(reason)] <- "Excluded"
  data.frame(group = group, reason = reason)
}

# ---------------------------------------------------------------------------------------
# Description of the groups
# ---------------------------------------------------------------------------------------

# Table 1: one column per group. Tardive dyskinesia is read in the text of the diagnoses, because
# the Tardive_dyskinesia column of the file is inconsistent (see metadata_agreement()).
make_table1 <- function(sel) {
  percent <- function(count, total) paste0(count, " (", round(100 * count / total), " %)")
  n <- table(sel$group)
  table1 <- data.frame(
    N                  = as.vector(n),
    Age_median         = tapply(sel$age, sel$group, median),
    Age_min            = tapply(sel$age, sel$group, min),
    Age_max            = tapply(sel$age, sel$group, max),
    Women              = percent(tapply(sel$sex == "female", sel$group, sum), n),
    Atherosclerosis    = percent(tapply(sel$ASCVD_status == "Yes", sel$group, sum), n),
    Diabetes           = percent(tapply(grepl("diabetes", sel$disease), sel$group, sum), n),
    Tardive_dyskinesia = percent(tapply(grepl("tardive dyskinesia", sel$disease), sel$group, sum), n),
    Nuclei_median      = tapply(sel$n_nuclei, sel$group, median))
  t(table1)
}

# Agreement between a yes/no column of the file and the text of the diagnoses.
# A column with discordant donors is not used.
metadata_agreement <- function(donors, column_name, searched_word) {
  column_yes <- as.character(donors[[column_name]]) == "Yes"
  text_yes   <- grepl(searched_word, donors$disease)
  data.frame(column             = column_name,
             positive_in_column = sum(column_yes),
             positive_in_text   = sum(text_yes),
             positive_in_both   = sum(column_yes & text_yes),
             discordant         = sum(column_yes != text_yes))
}

# Number of donors with at least nuclei_min nuclei of each subclass, per group (criterion E7).
donors_per_subclass <- function(nuclei_donor_subclass, group_of_rows, nuclei_min) {
  enough <- nuclei_donor_subclass >= nuclei_min
  sapply(c("Control", "SCZ", "BD"), function(g) colSums(enough[group_of_rows == g, , drop = FALSE]))
}
