# GLP-1 Receptor Agonist Weight Regain Analysis

This repository contains the R code used for *Trajectory of weight regain after cessation of GLP-1 receptor agonists: a systematic review and nonlinear meta-regression*. The code fits nonlinear mixed-effects recovery models and generates reproducible plots of weight regain with fixed effect and study-specific curves.

## Contents
- **R/** — Functions for data preparation, model fitting, and plotting
- **scripts/** — Analysis scripts calling the functions above 
- **data/** — Contains the dataset (`data_sheet.csv`) used in the analysis
- **figures/** — Output figures (main analysis plot, etc.)

## Requirements
- R 4.3.2
- Packages:  
  - tidyverse  
  - nlme  
  - ggplot2  
  - ggsci  
  - metafor  
- Package versions are managed with [`renv`](https://rstudio.github.io/renv/).  

## Quick start

Clone this repository, then open it in R (RStudio recommended).  
To set up the environment and install all required packages:

```r
# Install renv if you don’t already have it
install.packages("renv")

# Restore the project packages from renv.lock
renv::restore()
```

Run the main analysis script: this will fit the basic non-linear model and save the figure to `figures/graph_main.png`.  
```r
source("scripts/main.R")
```

## Data

The dataset (`data/data_sheet.csv`) contains extracted summary statistics from published clinical trial reports of GLP-1 receptor agonists. It is provided solely to reproduce the analyses described in our manuscript.  

The dataset does not contain individual patient-level data. Original trial data are owned by the respective sponsors and are not redistributed here.

## License
This code is released under the MIT License. See the LICENSE file for details.
"# glp1_weight_regain_analysis" 
