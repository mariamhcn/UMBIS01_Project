# UMBIS01 - Single-Nucleus Transcriptomic Analysis of Schizophrenia and Bipolar Disorder using PsychAD cohort
**Cell-type-specific molecular alterations, functional convergence, intercellular communication and transcriptomic ageing in the human prefrontal cortex**

## 1. Project overview
Schizophrenia (SCZ) and bipolar disorder (BD) are complex psychiatric disorders with partially overlapping genetic and molecular features. However, the extent to which they share transcriptional alterations across neuronal and non-neuronal cell populations remains unclear.
This project investigates these alterations using publicly available single-nucleus RNA sequencing (snRNA-seq) data from the PsychAD atlas, focusing on the Human Brain Collection Core (HBCC) cohort.
Main research question : Do schizophrenia and bipolar disorder share cell-type-specific molecular alterations in the human prefrontal cortex, and could these changes be associated with ageing-like transcriptomic profiles and differences in inferred neuronal-glial communication ?
Five complementary analyses address cellular composition (Q1), differential expression (Q2), functional enrichment (Q3), cell-cell communication (Q4) and transcriptomic ageing (Q5).


## 2. Dataset
- **Source :** PsychAD Human Prefrontal Cortex Atlas
- **Cohort :** Human Brain Collection Core (HBCC)
- **Region :** Dorsolateral prefrontal cortex
- **Technology :** snRNA-seq
- **Data access :** CZ CELLxGENE - PsychAD data collection
  
The initial HBCC dataset comprises approximately 1.49 million nuclei from 300 donors. **Following predefined clinical, diagnostic and age-related selection criteria, 224 donors aged 20 years or older were retained.**

| Diagnosis | Donors |
|---|---|
| Control | 121 |
| Schizophrenia (SCZ) | 49 |
| Bipolar disorder (BD) | 54 |
| Total | 224 |

The original PsychAD annotations were preserved, covering 27 cell subclasses.
Preprocessing generated two complementary files :
- **A Seurat object**, subsampled to a maximum of 200 nuclei per donor, used for visualization and cell-cell communication inference.
- **A pseudobulk count matrix**, generated from all available nuclei by aggregating raw counts per donor and cell subclass, used for donor-level statistical analyses.
The full pseudobulk matrix contains 34 176 genes and 6 048 potential donor-subclass profiles, including empty combinations.


## 3. Analyses and main findings

| Analysis | Methods | Main findings |
|---|---|---|
| **Q1 - Cell composition** | Donor-level proportions, `propeller` / `limma` | Lower **relative microglial proportion** in SCZ and BD, with stronger sensitivity-analysis support in SCZ; a less robust oligodendrocyte signal in BD. No significant direct SCZ-BD composition difference. |
| **Q2 - Differential expression** | Donor-level pseudobulk, `DESeq2`, three contrasts per subclass | **23 subclasses**, **4,042 eligible pseudobulk profiles**, **69 subclass-contrast comparisons**. `EN_L2_3_IT` in SCZ: 2 455 DEGs; astrocytes in BD: 1 103 DEGs. Shared but heterogeneous signatures; microglia show the strongest global SCZ–BD effect-size concordance (Spearman ρ = 0.535). |
| **Q3 - Functional enrichment** | ORA and GSEA (`clusterProfiler`); GO, Hallmark, Reactome, SynGO | Convergence on inflammatory and stress-related pathways is strongest in **glial and vascular** populations. **Synaptic programs differ** by disorder and neuronal subclass. Shared pathways involve a partially shared core of contributing genes. |
| **Q4 - Cell–cell communication** | `CellChat` and `NeuronChat` on 5 100 sampled nuclei from 17 subclasses | Broadly preserved CellChat connectivity, with differences in inferred neuronal-oligodendroglial interactions and NRXN/ADGRL/PTPR/NCAM signaling. NeuronChat predicts **60, 16, and 48 active GABAergic sender–receiver pairs** in Control, SCZ, and BD, respectively. These are descriptive predictions, **not functional measurements or donor-level significance tests**. |
| **Q5 - Transcriptomic ageing** | Signature scores, donor-level models, age-disease gene concordance | Multiple neuronal and glial populations show **older-looking expression profiles** in SCZ and BD (e.g. `EN_L2_3_IT`: approximately +14 and +13 apparent years). **No evidence of a faster age-related rate of change** was detected; an ageing-like offset is not proof of accelerated biological ageing. |

### Q1 - Cell-type composition
Cell-type proportions were compared across diagnoses at the donor level, adjusting for age, sex and genetic ancestry.
Microglia showed a lower relative abundance in both SCZ and BD, with stronger sensitivity-analysis support in SCZ. A lower oligodendrocyte proportion was also detected in BD, although this finding was less robust. No significant composition difference was detected in the direct SCZ-BD comparison.
These results refer to relative proportions, not absolute cell loss.

### Q2 - Cell-type-specific differential expression
Two notebooks address complementary aspects of differential expression.
	**Q2_1 - Pseudobulk DESeq2**
After filtering, 4 042 donor-subclass profiles across 23 eligible subclasses were analyzed. DESeq2 models included age, sex, genetic ancestry and diagnosis, with three contrasts per subclass : SCZ versus Control, BD versus Control and SCZ versus BD.
The analysis identified heterogeneous transcriptional alterations across populations. For example, 2 455 DEGs were detected in EN_L2_3_IT neurons in SCZ versus controls, while 1,103 DEGs were detected in astrocytes in BD versus controls (FDR < 0.05).
	**Q2_2 - Shared and distinct signatures**
Cross-disease comparisons identified concordant shared DEGs in excitatory neurons, astrocytes, pericytes and microglia. Microglia showed particularly strong genome-wide effect-size concordance between SCZ and BD (Spearman ρ = 0.535).
Overall, SCZ and BD exhibit partially overlapping but cell-subclass-dependent transcriptional alterations. Significance in only one disease comparison does not establish disease specificity.

### Q3 - Functional enrichment
ORA and GSEA were performed using GO Biological Process, Hallmark, Reactome and SynGO gene sets.
Results indicate that SCZ and BD converge on inflammatory and stress-related pathways, particularly in glial and vascular populations. Shared processes include interferon responses, IL6/STAT3 and TNFA signaling, hypoxia and mTORC1-related pathways.
Synaptic alterations are more heterogeneous and differ according to diagnostic group and neuronal subclass.
The results support partial functional convergence, even when the significantly altered genes are not identical between disorders.

### Q4 - Cell-cell communication
CellChat and NeuronChat were applied to a balanced dataset comprising 17 sufficiently represented subclasses and 5 100 nuclei.
CellChat identified broadly preserved inferred network connectivity, alongside differences in communication strength involving neurons, OPCs and oligodendrocytes. Prominent pathways included NRXN, ADGRL, PTPR and NCAM.
NeuronChat highlighted reduced inferred GABAergic connectivity, particularly in SCZ (16 active sender-receiver pairs, compared with 60 in controls and 48 in BD). Glutamatergic networks retained 90 active pairs in each diagnosis, with greater inferred communication weight in SCZ.
These analyses are exploratory and reflect predicted molecular interactions, not directly measured functional signaling or donor-level statistical differences.

### Q5 - Transcriptomic ageing
Ageing associated gene signatures and gene level age-disease concordance were evaluated across cellular populations.
Results suggest an older-looking transcriptomic profile in several neuronal and glial populations in SCZ and BD, including upper-layer excitatory neurons, interneurons and oligodendrocyte-lineage populations.
For example, the estimated apparent-age shifts included approximately +14 years in SCZ and +13 years in BD for EN_L2_3_IT, and approximately +21 and +22 years, respectively, for OPCs.
However, no significant increase in the rate of age-related change was detected. The findings therefore support an ageing-like expression shift rather than demonstrated accelerated biological ageing.
These observations remain subject to the limitations of cross-sectional post-mortem transcriptomics and potential clinical confounding.


## 4. Repository structure

The repository is organized as follows :
```text
UMBIS01_Project/                     # GitHub repository (local folder may be renamed)
├── 01_DATA/
│   ├── genesets/                  # Functional and ageing gene sets
│   ├── markers/                   # Marker reference lists
│   └── processed/                 # Processed data; large objects may be absent
├── 02_CODE/
│   ├── scripts/                   # Shared R functions, themes and color palettes
│   ├── 01_preprocessing.qmd
│   ├── 02_Q1_composition.qmd
│   ├── 03_Q2_1_pseudobulk_DE.qmd
│   ├── 03_Q2_2_shared_specific_signatures.qmd
│   ├── 04_Q3_pathway_enrichment.qmd
│   ├── 05_Q4_cell_communication.qmd
│   └── 06_Q5_transcriptomic_aging.qmd
├── 03_RESULTS/
│   ├── README.md                  # Missing large .rds files and regeneration instructions
│   ├── preprocessing/
│   ├── pseudobulk/
│   ├── Q1_composition/
│   ├── Q2_differential_expression/
│   ├── Q3_enrichment/
│   ├── Q4_cell_communication/
│   └── Q5_aging/
├── 04_FIGURES/
│   ├── preprocessing/
│   ├── Q1_composition/
│   ├── Q2_differential_expression/
│   ├── Q3_enrichment/
│   ├── Q4_cell_communication/
│   └── Q5_aging/
├── 05_DOCS/                       # Final report and presentation slides
├── renv/
│   ├── activate.R                 # renv activation script
│   ├── settings.json              # renv project configuration
│   └── .gitignore                 # Excludes renv/library/ and local cache files
├── .Rprofile                      # Automatically activates renv in the R project
├── renv.lock                      # Exact R package versions and sources
├── .gitignore
├── PsychAD_Project.Rproj
└── README.md
```

The locally installed package library (`renv/library/`) and some large input/output files are not committed to Git. The Quarto notebooks have corresponding HTML reports. See [`03_RESULTS/README.md`](03_RESULTS/README.md) for information on omitted intermediate results.

### Directory descriptions
	•	01_DATA/: input metadata, gene sets, marker lists and processed data.
	•	02_CODE/: executable Quarto notebooks, corresponding HTML reports and reusable R scripts.
	•	02_CODE/scripts/: shared helper functions for data processing, statistical analyses and consistent figures, including graphical themes and color palettes.
	•	03_RESULTS/: statistical tables, analysis summaries and saved intermediate results.
	•	04_FIGURES/: exported plots, grouped by analysis.
	•	05_DOCS/: final scientific report and oral presentation slides.


## 5. Data availability and large files
Not all input data and intermediate results could be uploaded to GitHub because of their large file sizes.
In particular, the original HBCC .h5ad dataset, large processed objects (.rds) and some computationally intensive analysis outputs may be absent from the published repository.
Their absence does not necessarily indicate missing analyses. The corresponding procedures are documented in the Quarto notebooks and supporting R scripts.
These files can be regenerated by executing the relevant code, provided that the original datasets, required packages, dependencies and sufficient computing resources are available.
The published HTML notebooks contain recorded analyses, tables, figures and interpretations, allowing results to be inspected without repeating every computational step.
See the proposed 01_DATA/README.md for additional information on input files and regeneration.
For further details on omitted intermediate results, including large CellChat and NeuronChat `.rds` files and instructions for their regeneration, see [03_RESULTS/README.md](03_RESULTS/README.md).

## 6. Reproducibility

### Software and package environment

The analyses were performed in **R**, using RStudio for project development and **Quarto** for notebook rendering. The committed [`renv.lock`](renv.lock) records the project's R version (**4.6.0**), Bioconductor version (**3.23**), and the versions and sources of the R packages needed to reproduce the analyses, including the GitHub versions of **CellChat** and **NeuronChat**.

The repository includes `renv.lock`, `.Rprofile`, `renv/activate.R` and `renv/settings.json`. The installed package library (`renv/library/`) is intentionally excluded from Git.

### Restoring the R environment

After cloning the repository, open `PsychAD_Project.Rproj` in RStudio. From the **project root**, run the following in the R console:

```r
install.packages("renv")  # Only if renv is not already available
renv::restore()
renv::status()
```

`renv::restore()` installs the recorded package versions; **do not run `renv::init()` again**. Quarto, the original scientific datasets, operating-system libraries and sufficient memory are still required separately. Package compatibility is easiest to reproduce with the R version recorded in the lockfile.

### Reproducing the workflow

The source notebooks are stored in `02_CODE/`. Render them in the following order after restoring the dependencies and obtaining the original HBCC input data:

```bash
cd 02_CODE
quarto render 01_preprocessing.qmd --to html
quarto render 02_Q1_composition.qmd --to html
quarto render 03_Q2_1_pseudobulk_DE.qmd --to html
quarto render 03_Q2_2_shared_specific_signatures.qmd --to html
quarto render 04_Q3_pathway_enrichment.qmd --to html
quarto render 05_Q4_cell_communication.qmd --to html
quarto render 06_Q5_transcriptomic_aging.qmd --to html
```

Some notebooks load helper scripts via relative paths such as `scripts/utils_figures.R`, so their expected execution directory is `02_CODE/`. Other inputs and outputs are resolved relative to the repository root using `here::here()`.

**Large files : ** some `.rds` inputs and intermediate results are deliberately excluded from GitHub. See [`03_RESULTS/README.md`](03_RESULTS/README.md) for regeneration instructions. Q4's CellChat and NeuronChat inference chunks are disabled by default (`#| eval: false`); change them to `#| eval: true` when regenerating the network `.rds` files. Rendering Q4 without those files and without re-enabling inference will fail. The original/processed HBCC inputs must also be available or regenerated first.

The analyses use documented filtering criteria, covariates and fixed random seeds where appropriate. Q4's pooled-network comparisons are exploratory predictions, not donor-level significance tests.

### Reproducibility considerations

- Donors are the biological units for statistical disease comparisons.
- Raw counts are preserved for count-based analyses.
- Filtering criteria, covariates and significance thresholds are documented.
- Random seeds are specified for stochastic procedures.
- Reusable plotting functions promote consistency across analyses.
- Model outputs and figures are saved in dedicated directories.
- Computationally inferred communication is distinguished from statistically tested gene-expression differences.
- The package environment is recorded in `renv.lock`, and notebook sessions are documented through `sessionInfo()`.

## 7. Project organization

- Myriam S. : Donor selection, preprocessing and Q1
- Mariam H. : Q2 differential expression and Q4 cell-cell communication
- Ann-Solène A. : Q3 functional enrichment and Q5 transcriptomic ageing

All analyses were discussed and reviewed collaboratively.

## 8. Main limitations

- The dataset is observational and post-mortem.
- Differences between diagnostic groups may be affected by genetic ancestry, age, treatment history, cause of death and other unmeasured factors.
- Cell-type proportions do not directly measure cell loss. Gene-expression differences and enriched pathways do not establish causal mechanisms or functional changes.
- CellChat and NeuronChat infer potential signaling networks without functional validation or donor-level statistical adjustment.
- Ageing-like transcriptional profiles do not demonstrate an increased biological ageing rate.
- Findings are interpreted as disease-associated molecular patterns rather than causal or experimentally validated mechanisms.

## 9. References
The following publications and methodological resources are cited in the project notebooks. The list is **deduplicated**, citation numbering within each notebook remains local to that notebook.

1. Fullard J. F. et al. (2025). Population-scale cross-disorder atlas of the human prefrontal cortex at single-cell resolution. *Scientific Data*, 12, 954. https://doi.org/10.1038/s41597-025-04687-5
2. Lee D. et al. (2026). Single-cell atlas of transcriptomic vulnerability across brain disorders. *Nature*, 657, 988–1002. https://doi.org/10.1038/s41586-025-09573-z
3. Lee D. et al. (2024). Single-cell atlas of transcriptomic vulnerability across multiple neurodegenerative and neuropsychiatric diseases. *medRxiv*. https://doi.org/10.1101/2024.10.31.24316513
4. Girdhar K. et al. PsychAD lifespan atlas of the human prefrontal cortex. *medRxiv*. https://doi.org/10.1101/2024.11.06.24316592
5. Yang H. et al. (2026). Lifespan single-cell transcriptomic atlas of the human prefrontal cortex. *Nature*, 657, 1003–1015. https://doi.org/10.1038/s41586-026-10271-7
6. Hao Y. et al. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology*, 42, 293–304.
7. Ma S. et al. (2022). Molecular and cellular evolution of the primate dorsolateral prefrontal cortex. *Science*, 377, eabo7257.
8. Squair J. W. et al. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications*, 12, 5692. https://doi.org/10.1038/s41467-021-25960-2
9. Phipson B. et al. (2022). propeller: testing for differences in cell type proportions in single cell data. *Bioinformatics*, 38, 4720–4726.
10. Ritchie M. E. et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43, e47.
11. Love M. I., Huber W. & Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. https://doi.org/10.1186/s13059-014-0550-8
12. Love M. I. et al. *Analyzing RNA-seq data with DESeq2*. Bioconductor vignette. https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
13. Benjamini Y. & Hochberg Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B*, 57, 289–300.
14. Gandal M. J. et al. (2018). Shared molecular neuropathology across major psychiatric disorders parallels polygenic overlap. *Science*, 359, 693–697. https://doi.org/10.1126/science.aad6469
15. Ruzicka W. B. et al. (2024). Single-cell multi-cohort dissection of the schizophrenia transcriptome. *Science*, 384, eadg5136. https://doi.org/10.1126/science.adg5136
16. Harvard Chan Bioinformatics Core. *Functional analysis of pseudobulk DE*. https://hbctraining.github.io/Pseudobulk-for-scRNAseq-Quarto/lessons/10_functional_analysis_pseudobulk.html
17. Wu T. et al. (2021). clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *The Innovation*, 2, 100141.
18. Subramanian A. et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 102, 15545–15550.
19. Zhou Y. et al. (2019). Metascape provides a biologist-oriented resource for the analysis of systems-level datasets. *Nature Communications*, 10, 1523.
20. Liberzon A. et al. (2015). The Molecular Signatures Database hallmark gene set collection. *Cell Systems*, 1, 417–425.
21. Koopmans F. et al. (2019). SynGO: an evidence-based, expert-curated knowledge base for the synapse. *Neuron*, 103, 217–234.
22. Trubetskoy V. et al. (2022). Mapping genomic loci implicates genes and synaptic biology in schizophrenia. *Nature*, 604, 502–508.
23. Mullins N. et al. (2021). Genome-wide association study of more than 40,000 bipolar disorder cases provides new insights into the underlying biology. *Nature Genetics*, 53, 817–829.
24. Lanz T. A. et al. (2019). Postmortem transcriptional profiling reveals widespread increase in inflammation in schizophrenia. *Translational Psychiatry*, 9, 151.
25. Ling E. et al. (2024). A concerted neuron–astrocyte program declines in ageing and schizophrenia. *Nature*, 627, 604–611.
26. Tomita H. et al. (2004). Effect of agonal and postmortem factors on gene expression profile. *Biological Psychiatry*, 55, 346–352.
27. Jin S. et al. (2021). Inference and analysis of cell-cell communication using CellChat. *Nature Communications*, 12, 1088. https://doi.org/10.1038/s41467-021-21246-9
28. Jin S., Plikus M. V. & Nie Q. (2025). CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics. *Nature Protocols*, 20, 180–219. https://doi.org/10.1038/s41596-024-01045-4
29. Zhao W. et al. (2023). Inferring neuron-neuron communications from single-cell transcriptomics through NeuronChat. *Nature Communications*, 14, 1128. https://doi.org/10.1038/s41467-023-36800-w
30. Okuda T. et al. (2024). Alterations in inhibitory neuron subtype-selective transcripts in the prefrontal cortex: comparisons across schizophrenia and mood disorders. *Psychological Medicine*, 54, 3896–3905.
31. Dienel S. J., Fish K. N. & Lewis D. A. (2023). The nature of prefrontal cortical GABA neuron alterations in schizophrenia: markedly lower somatostatin and parvalbumin gene expression without missing neurons. *American Journal of Psychiatry*, 180, 495–507. https://doi.org/10.1176/appi.ajp.20220676
32. Gao R. et al. (2026). Postmortem brain single-nucleus and bulk gene expression analyses identify shared and distinct abnormalities in bipolar disorder and major depressive disorder. *Translational Psychiatry*, 16, 442.
33. López-Otín C. et al. (2023). Hallmarks of aging: an expanding universe. *Cell*, 186, 243–278.
34. Mattson M. P. & Arumugam T. V. (2018). Hallmarks of brain aging: adaptive and pathological modification by metabolic states. *Cell Metabolism*, 27, 1176–1199.
35. Kaufmann T. et al. (2019). Common brain disorders are associated with heritable patterns of apparent aging of the brain. *Nature Neuroscience*, 22, 1617–1623.
36. McKinney B. C. et al. (2017). DNA methylation evidence against the accelerated aging hypothesis of schizophrenia. *npj Schizophrenia*, 3, 13.
37. Fröhlich A. S. et al. (2024). Single-nucleus transcriptomic profiling of human orbitofrontal cortex reveals convergent effects of aging and psychiatric disease. *Nature Neuroscience*, 27, 2021–2032.
38. Saul D. et al. (2022). A new gene set identifies senescent cells and predicts senescence-associated pathways across tissues. *Nature Communications*, 13, 4827.
39. Galatro T. F. et al. (2017). Transcriptomic analysis of purified human cortical microglia reveals age-associated changes. *Nature Neuroscience*, 20, 1162–1171.
40. The Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. *Genetics*, 224, iyad031.
41. Ciernik L. et al. (2026). Robust and efficient annotation of cell states through gene signature scoring. *Genome Research*, 36, 630–644.
42. Fieller E. C. (1954). Some problems in interval estimation. *Journal of the Royal Statistical Society B*, 16, 175–185.
