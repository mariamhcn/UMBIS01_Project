# CellChat v2 — From Ligand–Receptor Interactions to Comparative Cell–Cell Communication Networks
***An introductory, reproducible R Markdown tutorial for inferring and visualizing cell–cell communication with CellChat, using a publicly available human skin single-cell RNA-seq dataset.***

The tutorial introduces the core CellChat workflow on non-lesional (NL) skin, then extends it to an exploratory comparison with lesional (LS) skin.

## Objectives

By the end of this tutorial, you will be able to :
- Prepare an annotated single-cell expression matrix for CellChat,
- Select and explore the human ligand–receptor interaction database,
- Infer ligand–receptor interactions and pathway-level communication networks,
- Identify prominent signaling senders, receivers, and pathways,
- Visualize and descriptively compare communication networks across two conditions,
- Understand the methodological limitations of transcriptomics-based communication inference.
  
## Dataset

We use the human skin scRNA-seq example dataset provided by the CellChat developers, containing normalized gene expression data and annotated cell metadata.
| Condition | Cells |
| --- | ---: |
| Non-lesional skin (NL) | 2,552 |
| Lesional skin (LS) | 5,011 |
| **Total** | **7,563** |

The expression matrix contains 17,328 genes, and cells are assigned to 12 annotated populations, including fibroblast, dendritic-cell, and T-cell populations.

**Original dataset:** [CellChat human skin dataset on Figshare](https://ndownloader.figshare.com/files/25950872)

The R Markdown script automatically downloads the .rda file to 01_DATA/ if it is absent. Internet access is required for the initial download. The dataset has already been preprocessed, normalized, and annotated, this tutorial does not repeat upstream single-cell quality control or clustering.


## Tutorial contents

**Part 1 — Core CellChat workflow (NL skin)**

1. Load libraries and the dataset.
2. Create a CellChat object.
3. Explore and select CellChatDB.human.
4. Identify overexpressed signaling genes and ligand–receptor interactions.
5. Infer and filter communication probabilities.
6. Compute signaling pathways and aggregate communication networks.
7. Visualize networks with circle plots, heatmaps, and bubble plots.
8. Explore sender/receiver roles and save results.
   
**Part 2 — Comparative analysis (NL vs LS skin)**
1. Reconstruct the lesional-skin communication network using comparable parameters.
2. Merge the independently analyzed NL and LS CellChat objects.
3. Compare inferred interaction counts and strengths.
4. Examine differences in signaling pathways, with a focus on CXCL.
5. Summarize findings and discuss limitations.
   
*Part 2 is an optional extension of the main tutorial.*


## Requirements
- R and RStudio (or another environment able to render R Markdown).
- Internet access to install packages and download the dataset.
- Sufficient memory and disk space to process the example dataset.
- The R packages specified in renv.lock, including CellChat and the packages used to render the notebook.
This project uses renv to record package versions. Restoration may require compatible versions of R, system libraries, and build tools. Packages installed from GitHub or Bioconductor may have platform-specific requirements.


## How to reproduce the analysis

1. Obtain the repository. Clone it with Git, or select Code → Download ZIP on GitHub and extract it.
```bash
git clone <YOUR-REPOSITORY-URL>
cd <YOUR-REPOSITORY-FOLDER>
```

2. Open the project in RStudio with the included .Rproj file. Make sure the working directory is the repository root (where the .Rmd and renv.lock files are located).

3. Restore the R environment. In the R console:
```r
install.packages("renv") # only if renv is not yet available
renv::restore()
renv::status()
```
If `renv::restore()` cannot install a GitHub/Bioconductor dependency, consult the package's official installation instructions and check that the lockfile records the correct source.

4. Open CellChat_Tutorial_UMBIS01.Rmd and click Knit.

The script downloads the dataset if necessary, executes the analysis, and produces the HTML tutorial. Its relative paths assume execution from the repository root.


## Repository structure

```text
CellChat_Tutorial/
├── CellChat_Tutorial.Rproj          # RStudio project
├── CellChat_Tutorial_UMBIS01.Rmd   # Executable tutorial source
├── CellChat_Tutorial_UMBIS01.html  # Rendered tutorial (read without R)
├── README.md                       # Project description and instructions
├── renv.lock                       # Locked package versions (root directory)
├── .Rprofile                       # Automatic renv activation
├── .gitignore                      # Excludes temporary files and local libraries
├── renv/                           # renv activation/configuration scripts
├── 01_DATA/                        # Input dataset
├── 02_RESULTS/                     # Saved RDS objects and CSV tables
├── 03_FIGURES/                     # Generated figures (organized by part)
└── 04_REPRODUCIBILITY/             # Reproducibility notes and session information
```

## Interpretation and limitations
CellChat infers potential communication based on expression patterns and a curated ligand–receptor database. Its predictions do not directly demonstrate ligand secretion, receptor activation, physical proximity, or functional signaling.

The NL and LS datasets contain different cell numbers and population distributions. Additionally, populations below the `min.cells = 10` threshold are excluded from inferred communications: in NL, this affects cDC1 and inflammatory dendritic cells. Consequently, the cross-condition comparisons are descriptive and exploratory, not donor-level statistical tests of disease effects. Biological conclusions require independent validation.


## References :
[1] Jin S, Plikus MV, Nie Q. CellChat for systematic analysis of cell-cell communication from single-cell transcriptomics. Nat Protoc. 2025 Jan;20(1):180-219. doi: 10.1038/s41596-024-01045-4. Epub 2024 Sep 16. PMID: 39289562.

[2] Jin, S., Guerrero-Juarez, C.F., Zhang, L. et al. Inference and analysis of cell-cell communication using CellChat. Nat Commun 12, 1088 (2021). https://doi.org/10.1038/s41467-021-21246-9

[3] Official CellChat repository : Interface between mathematics, AI and biomedicine. (2026). Jinworks/cellchat [R]. https://github.com/jinworks/CellChat (Original Edition 2023)

[4] Original example dataset : https://ndownloader.figshare.com/files/25950872
