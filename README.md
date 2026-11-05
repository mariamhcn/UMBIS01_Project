# PsychAD — Single-Nucleus Transcriptomic Analysis of Schizophrenia and Bipolar Disorder

### Exploring cell-type-specific molecular alterations, intercellular communication and aging-associated signatures in the human prefrontal cortex

## 1. Project Overview

Schizophrenia (SCZ) and bipolar disorder (BD) are complex psychiatric disorders that share certain clinical and genetic features. However, the extent to which they exhibit common or distinct molecular alterations across neuronal and glial populations remains an important research question.

This project investigates transcriptomic alterations associated with SCZ and BD using publicly available human brain single-nucleus RNA sequencing (snRNA-seq) data.

Our analysis focuses on the **Human Brain Collection Core (HBCC) cohort** from the **PsychAD atlas**, available through CZ CELLxGENE.

**Main research question:**

*Do schizophrenia and bipolar disorder share cell-type-specific molecular alterations in the human prefrontal cortex, and could these changes reflect accelerated aging-like transcriptomic profiles and affect predicted neuronal–glial signaling networks?*

The project combines differential gene expression, functional enrichment, cell–cell communication inference and aging-related transcriptional signatures.

## 2. Dataset

**Source:** PsychAD — Human Prefrontal Cortex Single-Cell Atlas  
**Cohort:** HBCC  
**Brain region:** Dorsolateral prefrontal cortex (DLPFC)  
**Technology:** Single-nucleus RNA sequencing (snRNA-seq)

**Disease groups:**
- Control donors
- Schizophrenia (SCZ)
- Bipolar disorder (BD)

The original HBCC cohort comprises approximately 300 donors and 1.4 million nuclei before applying project-specific inclusion criteria.

### Donor Selection

Our main analysis focuses on adult donors aged **20 years or older**, with additional clinical and diagnostic exclusion criteria.

Donors with neurodegenerative conditions or relevant secondary psychiatric/neurological diagnoses are excluded according to the predefined selection strategy. Bipolar-type schizoaffective diagnoses are excluded from the main SCZ–BD comparison to maintain mutually exclusive groups.

Age-threshold sensitivity analyses (≥17 and ≥24 years) are planned to evaluate the robustness of the results.

The final donor counts and eligibility criteria will be documented after dataset preprocessing.

**Dataset access:**  
https://cellxgene.cziscience.com/collections/84ce6837-548d-4a1f-919f-0bc0d9a3952f

## 3. Scientific Objectives

The project is organized around five complementary research questions.

| Question | Objective | Main methods |
|---|---|---|
| Q1 — Cellular alterations | Identify neuronal and glial populations showing disease-associated changes | Cell-type composition, donor-level proportions |
| Q2 — Differential expression | Identify shared and disease-associated transcriptional signatures | Pseudobulk, DESeq2 |
| Q3 — Functional convergence | Determine whether distinct gene signatures converge on common biological pathways | GSEA, ORA, clusterProfiler |
| Q4 — Cell–cell communication | Investigate changes in predicted neuronal–glial signaling networks | CellChat v2 |
| Q5 — Transcriptomic aging | Assess whether disease-associated expression profiles resemble age-related molecular changes | Gene signature scoring, donor-level regression |

Together, these analyses explore disease-associated alterations at the cellular, molecular, functional and intercellular signaling levels.

## 4. Analysis Workflow

### Step 1 — Data Preprocessing and Quality Control

The preprocessing pipeline includes donor selection, quality assessment, metadata validation, normalization and examination of cell-type annotations.

The PsychAD reference annotations are used as the starting point, rather than performing de novo annotation.

Batch effects and integration requirements will be assessed before applying any correction.

A common Seurat object will be generated for downstream analyses while preserving the original raw expression counts.

**Main tools:** Seurat, SeuratObject, zellkonverter, ggplot2.

### Step 2 — Cell-Type Composition (Q1)

Cellular populations will be compared across Control, SCZ and BD groups.

Analyses will focus on donor-level cell-type proportions and the identification of potentially disease-associated populations.

### Step 3 — Differential Gene Expression (Q2)

Raw gene counts will be aggregated by donor and cell type using a **pseudobulk approach**.

Differential expression will be tested using DESeq2, with adjustments for available demographic and technical covariates, including age, sex and post-mortem interval (PMI).

Three comparisons will be considered:

- SCZ vs Control
- BD vs Control
- SCZ vs BD

Gene-level effects will be compared to identify concordant, discordant and disease-associated transcriptional signatures.

### Step 4 — Functional Enrichment (Q3)

Differential expression results will be used to investigate enriched biological processes.

**Gene Set Enrichment Analysis (GSEA)** will be the primary approach, complemented by Over-Representation Analysis (ORA).

Gene sets will include Gene Ontology Biological Process, MSigDB Hallmark, Reactome and synapse-related annotations, including SynGO.

The aim is to determine whether SCZ and BD affect similar biological functions, even when their differentially expressed genes differ.

### Step 5 — Cell–Cell Communication (Q4)

CellChat v2 will be used to infer and compare potential intercellular signaling interactions across Control, SCZ and BD groups.

The analysis will focus on predicted ligand–receptor interactions, communication strengths and signaling pathways among neuronal and glial populations.

NeuronChat may be considered as an additional analysis focusing on neurotransmitter-mediated signaling.

These results will be interpreted as computational predictions rather than experimentally validated communications.

### Step 6 — Aging-Associated Transcriptomic Signatures (Q5)

Published gene signatures associated with cellular senescence, inflammation, mitochondrial function, myelination and brain aging will be evaluated across cell populations.

Signature scores will be aggregated at the donor level and examined in relation to chronological age and psychiatric diagnosis.

The goal is to identify **accelerated-aging-like transcriptional profiles**, without claiming to measure biological age directly or establish causal acceleration of aging.

## 5. Repository Structure

```text
PsychAD_Project/
│
├── 01_DATA/
│   └── Data references, metadata and gene sets
│
├── 02_CODE/
│   └── R scripts and analysis notebooks
│
├── 03_RESULTS/
│   └── Statistical results and output tables
│
├── 04_FIGURES/
│   └── Plots and scientific visualizations
│
├── 05_DOCS/
│   └── Methodological notes and documentation
│
└── README.md
```

The directories will be populated progressively as the analyses are implemented.

Large raw sequencing datasets and intermediate files may be stored externally when they exceed GitHub storage limits.

## 6. Software and Reproducibility

The project is primarily developed in **R**, using methods introduced during the single-cell transcriptomics course.

Main planned packages include:

- **Seurat:** single-nucleus data processing and visualization.
- **DESeq2:** pseudobulk differential expression analysis.
- **clusterProfiler:** pathway enrichment analysis.
- **CellChat:** inference of cell–cell communication networks.
- **ggplot2 / dplyr:** data manipulation and visualization.

Additional packages may be introduced depending on the requirements of each analysis.

Exact software versions, execution instructions, required datasets and reproducibility information will be documented as the workflow is finalized.

## 7. Project Organization

The project is conducted collaboratively, with responsibilities distributed across the main analyses.

| Team member | Responsibilities |
|---|---|
| Myriam | Data preprocessing, donor selection and Q1 |
| Mariam | Differential expression (Q2) and cell–cell communication (Q4) |
| Ann-Solène | Functional enrichment (Q3) and transcriptomic aging (Q5) |

All analyses are subject to collaborative review, methodological discussion and cross-validation.

## 8. Scientific Limitations

Several considerations will guide the interpretation of the results:

- Donors, rather than individual nuclei, must be considered the biological replicates for disease comparisons.
- Psychiatric diagnoses may be associated with differences in age, treatment exposure and other confounding factors.
- Post-mortem transcriptomic profiles do not directly establish causal mechanisms.
- Predicted intercellular communication does not demonstrate functional signaling.
- Aging-associated gene expression patterns do not necessarily indicate accelerated biological aging.

These limitations will be addressed through statistical adjustment, sensitivity analyses and transparent interpretation.

## 9. References

[1] Center for Disease Neurogenomics @ Mount Sinai. (2026). DiseaseNeuroGenomics/psychAD_SciData [Logiciel]. https://github.com/DiseaseNeuroGenomics/psychAD_SciData

[2] Fullard, J. F. et al. PsychAD Consortium, Girdhar, K., … Roussos, P. (2025). Population-scale cross-disorder atlas of the human prefrontal cortex at single-cell resolution. Scientific Data, 12(1), 954. https://doi.org/10.1038/s41597-025-04687-5

[3] Lee, D., Koutrouli, M., Masse, N.Y. et al. Single-cell atlas of transcriptomic vulnerability across brain disorders. Nature 657, 988–1002 (2026). https://doi.org/10.1038/s41586-025-09573-z

[4] Center for Disease Neurogenomics @ Mount Sinai. (2026). Diseaseneurogenomics/psychadxd [Jupyter Notebook]. https://github.com/DiseaseNeuroGenomics/PsychADxD (Édition originale 2024)

[5] Yang, H., Clarence, T., Scott, M.R. et al. Lifespan single-cell transcriptomic atlas of the human prefrontal cortex. Nature 657, 1003–1015 (2026). https://doi.org/10.1038/s41586-026-10271-7

[6] Yu, G. Biomedical Knowledge Mining using GOSemSim and clusterProfiler (official book of clusterProfiler).

[7] Harvard Chan Bioinformatics Core. Functional analysis of pseudobulk DE (HBC Training). https://hbctraining.github.io/Pseudobulk-for-scRNAseq-Quarto/lessons/10_functional_analysis_pseudobulk.html

[8] Wu, T. et al. (2021). clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. The Innovation, 2, 100141. https://doi.org/10.1016/j.xinn.2021.100141

[9] Subramanian, A. et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. PNAS, 102, 15545–15550. https://doi.org/10.1073/pnas.0506580102

[10] Liberzon, A. et al. (2015). The Molecular Signatures Database hallmark gene set collection. Cell Systems, 1, 417–425. https://doi.org/10.1016/j.cels.2015.12.004

[11] Koopmans, F. et al. (2019). SynGO: an evidence-based, expert-curated knowledge base for the synapse. Neuron, 103, 217–234. https://doi.org/10.1016/j.neuron.2019.05.002

[12] Ruzicka, W. B. et al. (2024). Single-cell multi-cohort dissection of the schizophrenia transcriptome. Science, 384, eadg5136. https://doi.org/10.1126/science.adg5136

[13] Jin S, Plikus MV, Nie Q. CellChat for systematic analysis of cell-cell communication from single-cell transcriptomics. Nat Protoc. 2025 Jan;20(1):180-219. doi: 10.1038/s41596-024-01045-4. Epub 2024 Sep 16. PMID: 39289562.

[14] Zhao, W., Johnston, K. G., Ren, H., Xu, X., & Nie, Q. (2023). Inferring neuron-neuron communications from single-cell transcriptomics through NeuronChat. Nature Communications, 14(1), 1128. https://doi.org/10.1038/s41467-023-36800-w

[15] Wei-BioMath. (2026). Wei-biomath/neuronchat [R]. https://github.com/Wei-BioMath/NeuronChat (Édition originale 2022)

[16] López-Otín, C. et al. (2023). Hallmarks of aging: an expanding universe. Cell, 186, 243–278. https://doi.org/10.1016/j.cell.2022.11.001

[17] Mattson, M. P. & Arumugam, T. V. (2018). Hallmarks of brain aging: adaptive and pathological modification by metabolic states. Cell Metabolism, 27, 1176–1199. https://doi.org/10.1016/j.cmet.2018.05.011 → hallmarks du vieillissement cérébral

[18] Kaufmann, T. et al. (2019). Common brain disorders are associated with heritable patterns of apparent aging of the brain. Nature Neuroscience, 22, 1617–1623. https://doi.org/10.1038/s41593-019-0471-7

[19] McKinney, B. C. et al. (2017). DNA methylation evidence against the accelerated aging hypothesis of schizophrenia. npj Schizophrenia, 3, 13. https://doi.org/10.1038/s41537-017-0026-4

[20] Fröhlich, A. S. et al. (2024). Single-nucleus transcriptomic profiling of human orbitofrontal cortex reveals convergent effects of aging and psychiatric disease. Nature Neuroscience, 27, 2021–2032. https://doi.org/10.1038/s41593-024-01742-z

[21] Yang, H., Clarence, T., Scott, M. R. et al. (2026). Lifespan single-cell transcriptomic atlas of the human prefrontal cortex. Nature, 657, 1003–1015. https://doi.org/10.1038/s41586-026-10271-7

[22] Saul, D. et al. (2022). A new gene set identifies senescent cells and predicts senescence-associated pathways across tissues. Nature Communications, 13, 4827. https://doi.org/10.1038/s41467-022-32552-1

[23] Galatro, T. F. et al. (2017). Transcriptomic analysis of purified human cortical microglia reveals age-associated changes. Nature Neuroscience, 20, 1162–1171. https://doi.org/10.1038/nn.4597

[24] The Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. Genetics, 224, iyad031. https://doi.org/10.1093/genetics/iyad031

[25] Hao, Y. et al. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology, 42, 293–304. https://doi.org/10.1038/s41587-023-01767-y

[26] Gonzales, M. M. et al. (2023). Senolytic therapy in mild Alzheimer’s disease: a phase 1 feasibility trial. Nature Medicine, 29, 2481–2488. https://doi.org/10.1038/s41591-023-02543-w

---

**Project status:** In development — analysis pipelines and results will be added progressively.
