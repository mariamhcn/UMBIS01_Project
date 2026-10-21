# Reproducibility Guide — CellChat v2 Tutorial

## 1. Overview

This directory contains supplementary information to support the reproducibility of the CellChat v2 tutorial.

All analyses were performed in R using CellChat v2. The same main inference parameters were applied to non-lesional (NL) and lesional (LS) skin.

## 2. Software Requirements and Environment Restoration

This project was developed with **R 4.6.0**, Bioconductor 3.23, and CellChat v2 on macOS Apple Silicon. Successful execution has also been tested with R 4.6.1 on another macOS computer.

The R package environment is managed through `renv`. The `renv.lock` file records the required package versions, including `presto`, which is used by CellChat for fast differential expression testing, and `yaml`, which supports R Markdown dependency parsing.

**Installation procedure**

1. Install R 4.6.x and RStudio.
2. Clone or download the complete GitHub repository.
3. Open `CellChat_Tutorial.Rproj` in RStudio.
4. If required, install `renv` using `install.packages("renv")`.
5. Restore the recorded dependencies using `renv::restore()`.
6. Verify the environment using `renv::status()`.
7. Render the tutorial using `rmarkdown::render("CellChat_Tutorial_UMBIS01.Rmd")`.

The project uses the CRAN release of `presto`, avoiding the need for a separate GitHub installation.

If package installation fails, check the R version, network access, and system compilation requirements. `renv` manages R packages but does not install R itself or external compilation tools.

**Important:** Only update `renv.lock` from a validated development environment. New users should restore the recorded dependencies rather than creating a new snapshot.

## 3. Computational Environment

The project uses **`renv`** to manage R package dependencies and record their versions.

The following files are located at the root of the repository, not inside this directory :

- `renv.lock`: records the R package dependencies, versions, and available installation sources,
- `.Rprofile`: activates the project environment when opening the project,
- `renv/activate.R`: provides the project environment activation script,
- `CellChat_Tutorial.Rproj`: RStudio project configuration.

The file `04_REPRODUCIBILITY/session_info.txt` provides supplementary information about the R version, platform, and packages loaded during the analysis.

## 4. Reproducing the Analysis

**Step 1 — Obtain the project**

Clone or download the complete GitHub repository.

**Step 2 — Open the RStudio project**

Open `CellChat_Tutorial.Rproj` from the root directory.

**Step 3 — Restore the R environment**

If necessary, install `renv` :

```r
install.packages("renv")
```

Restore the recorded package dependencies :

```r
renv::restore()
```

GitHub and Bioconductor dependencies, including CellChat, may require an internet connection and compatible system tools.

**Step 4 — Run the tutorial**

Open `CellChat_Tutorial_UMBIS01.Rmd` and click **Knit**, or run :

```r
rmarkdown::render("CellChat_Tutorial_UMBIS01.Rmd")
```

The analysis uses relative paths and downloads the official human skin dataset from Figshare if it is missing locally.

**Step 5 — Verify the outputs**

After successful execution, verify that the project generates the expected outputs in :

- `02_RESULTS/`: inferred communication results,
- `03_FIGURES/PART_1_NL/`: non-lesional skin figures,
- `03_FIGURES/PART_2_NL_vs_LS/`: comparative analysis figures,
- `04_REPRODUCIBILITY/session_info.txt`: computational session information.

The HTML output should reproduce the figures, tables, and interpretations presented in the notebook.

## 5. Reproducibility Considerations

The analysis uses CellChat inference settings across NL and LS, including `triMean` averaging, 100 bootstrap iterations, `population.size = FALSE`, and `min.cells = 10`.

Fixed random seeds are used where applicable.

Nevertheless, reproducibility of the computational workflow does not imply biological validation of predicted interactions. CellChat infers potential communication using transcriptomic expression and curated ligand–receptor information.

The unequal cell numbers and population composition between NL and LS, together with the absence of donor-level statistical comparisons, mean that condition-specific differences should be interpreted as exploratory.

## 6. Troubleshooting

If environment restoration fails, check that `renv` is activated and consult :

```r
renv::status()
```

If an installation error concerns CellChat or another GitHub/Bioconductor dependency, inspect the installation error and the corresponding package entry before changing the project environment.

If the Knit fails with an object-not-found error, verify that all chunks are executed in order and that the analysis does not depend on objects previously created in an interactive R session.
