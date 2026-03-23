# Beyond Food Cues: Safety Signals Facilitate Scavenger Mutualism

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19177396.svg)](https://doi.org/10.5281/zenodo.19177396)

This repository contains the raw data, metadata, and R scripts required to reproduce the statistical analyses and figures for the manuscript:

**Beyond food cues: safety signals facilitate scavenger mutualism under human disturbance and food abundance**  
*Ruitian Xu, Hongfang Wang, Nan Lyu, Yu Guan, Lei Bao, Jianping Ge*  

## 📂 Repository Structure

To ensure full reproducibility, the repository is organized into three main directories:

- **`Dataset/`**: Contains all raw and processed data necessary for the analyses.
- **`R_scripts/`**: Contains sequentially numbered R scripts corresponding to the methodological steps described in the manuscript.
- **`results/`**: The default output directory for all generated figures, statistical tables, and intermediate data (automatically created when running the scripts).

## 📊 Data Description (`Dataset/`)

- **`Dataset_S1_Experimental_Design_and_Metadata.csv`**: Metadata for the Latin square experimental design, including daily food abundance (0, 6, 18, 36 items) and sampling point locations (distance from the road).
- **`Dataset_S2_Time_Series_Individual_Counts/`**: Folder containing daily CSV files with raw time-series counts of crows and specific behavioral states of white-tailed eagles (Watching vs. Feeding) per 30-second frame.
- **`Dataset_S3_Crow_Feeding_Behavior/`**: Folder containing daily CSV files logging the total number of visible crows alongside the number of crows actively engaged in head-down feeding.
- **`Dataset_S4_Eagle_Foraging_Decisions.csv`**: Trial-level behavioral data summarizing eagle foraging success, specific timing metrics (arrival, feeding, and latency), and crow signal presence.
- **`Dataset_S5_Temporal_Sequences.csv`**: Daily arrival order sequence (AOS) and peak order sequence (EOS) metrics across different experimental locations.

## 💻 Code Description (`R_scripts/`)

The analysis workflow is divided into 7 standardized R scripts. Please run them in numerical order. Scripts will automatically read from `../Dataset/` and write outputs to `../results/`.

1. **`01_temporal_correlation.R`**: Performs log-transformations and linear regressions to extract residuals, followed by non-parametric Spearman rank correlation analysis of temporal behaviors.
2. **`02_parallel_analysis.R`**: Conducts causal mediation analysis and sensitivity tests to reveal the two-stage information transfer pathway.
3. **`03_willingness_decisions_and_efficiency.R`**: Uses Firth's penalized logistic regression and Cox proportional hazards survival analysis to assess eagle feeding willingness and decision latency.
4. **`04_mutualistic_foraging_behavior.R`**: Processes long-format data and fits binomial Generalized Linear Mixed Models (GLMMs) to test changes in crow feeding investment before and after eagle intervention.
5. **`05_crow_and_eagle_temporal_sequence.R`**: Applies LOESS smoothing to population time-series data to visualize temporal synchrony between species.
6. **`06_max_count_analysis_for_env_factors.R`**: Fits GLM with Poisson distributions to analyze the interactive effects of human disturbance and food abundance on maximum aggregation scales.
7. **`07_clm_analysis_for_env_factors.R`**: Performs Cumulative Link Models (CLM) and Bayesian ordinal regression to evaluate the regulatory effects of environmental factors on sequence metrics.

## 🛠️ Software & Dependencies

All analyses were conducted in **R (version 4.4.1 or higher)**. The following major R packages are required to run the scripts:

```R
# Data manipulation and visualization
install.packages(c("tidyverse", "ggplot2", "cowplot", "patchwork", "readr", "stringr"))

# Statistical modeling
install.packages(c("lme4", "survival", "survminer", "logistf", "mediation", "ordinal", "car", "emmeans"))

# Model diagnostics and Bayesian regression
install.packages(c("DHARMa", "brms", "tidybayes", "posterior"))

# Formatting and networks
install.packages(c("broom", "broom.mixed", "DiagrammeR", "DiagrammeRsvg", "rsvg"))
````

## 🚀 How to Run

1. Clone or download this repository to your local machine.
2. Open your R IDE (e.g., RStudio).
3. **Crucial**: Set your working directory to the `R_scripts/` folder.
   *(e.g., `setwd("path/to/Heterospecific-Information/R_scripts")`)*
4. Run the scripts sequentially from `01` to `07`. All high-resolution figures (PNG & PDF) and statistical tables will be automatically saved in the newly created `results/` folder.

## 📝 License & Contact

This project is licensed under the MIT License. If you have any questions regarding the data or code, please contact Ruitian Xu (<ruitian_xu@163.com>) or the corresponding authors listed in the manuscript.
