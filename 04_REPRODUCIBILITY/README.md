# Reproducibility

This tutorial was developed using R and CellChat v2. Package dependencies and versions are recorded in `renv.lock` to facilitate reproducibility.

## Requirements

	•	R and RStudio
	•	An internet connection to download the dataset and install dependencies
	•	A compatible R environment (see renv.lock for the original R version)
  
## Instructions

	1.	Clone or download this GitHub repository,
	2.	Open the project in RStudio,
	3.	Install renv if necessary : `install.packages("renv")`, 
	4.	Restore the project dependencies: `renv::restore()`,
	5.	Open CellChat_Tutorial_UMBIS01.Rmd and click Knit to generate the HTML tutorial.
  
The human skin dataset is automatically downloaded from the official CellChat Figshare resource if it is not already available locally.
The analysis generates communication network results and figures in the corresponding project directories.

NB : The workflow infers potential cell–cell communications from transcriptomic data. Comparisons between non-lesional and lesional skin are exploratory and should not be interpreted as statistically validated biological differences.
