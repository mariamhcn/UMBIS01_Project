# Marker gene files

Used by notebook `02_CODE/01_preprocessing.qmd` (Step 8) to check the annotation of the cell types.

| File | Content | Source |
|------|---------|--------|
| `PsychAD_SupplementaryTable10.csv` | Gene panel designed by the authors of the PsychAD atlas for their spatial transcriptomics experiment (Xenium, 100 genes), with the class, subclass and subtype that each gene marks. Columns: `tier`, `class`, `subclass`, `subtype`, `marker`, `Ensembl_ID`. The notebook keeps the subclass-level markers (`tier` = 1, empty `subtype`) of the non-neuronal subclasses. | Lee et al., PsychAD atlas, Supplementary Table 10 (file not modified) |

The other marker genes used in the notebook are written directly in its code, with their source:
Figure 2d of Lee et al. (neuronal subclasses) and Ma et al. 2022, *Science* 377, eabo7257 (taxonomy used by the authors to
name the cell types).
