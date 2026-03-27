# Pediatric B-ALL Survival Analysis (Peru)
This repository contains the R scripts used in the statistical analysis of my undergraduate thesis.
- 📄 Thesis language: Spanish  
- 💻 Code and documentation: English
# Survival analysis
This repository contains the source code developed for the survival analysis of the thesis entitled "Cytogenetic and molecular factors associated with overall and event-free survival in pediatric B-cell acute lymphoblastic leukemia".
## Content
- **Descriptive graphs**: To describe the sample in terms of summary indicators
- **Kaplan-Meier & Log-Rank**: To estimate survival functions of each variable (categorical)
- **Cox models**: To identify independent association with the survival among the covariates.
### Data Visualization
The following visualization tools were used in the analysis:
1. Interactive tables created with the `DT` package.
2. Box-and-whisker plots to summarize the distribution of continuous variables.
3. Kaplan–Meier survival curves customized using the `survminer` package.
## Reproducibility
All analyses were conducted in R. Required packages include:
- survival
- survminer
- DT
Session information can be provided upon request.
## Data Availability
Due to patient privacy restrictions, the dataset used in this study is not publicly available.
## Citation
If you use this repository, please cite:
Caro Retamozo, S. (2026). *Survival Analysis Pipeline for Clinical Data (Kaplan-Meier and Cox Models)* (Version 1.0.2). Zenodo. https://doi.org/10.5281/zenodo.19260913
