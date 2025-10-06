# Applied Bioinformatics Assessment

This repository contains solutions to the Applied Bioinformatics Assessment using **R** and **GitHub**.

## Files

### 1. `gene_expression_analysis.R`
- **Input:** `gene_expression.tsv`  
- **Outputs:** 
  - Console summaries (first rows, top genes, counts)  
  - `gene_expression_histogram.png` (histogram of mean gene expression)  

### 2. `growth_analysis.R`
- **Input:** `growth_data.csv`  
- **Outputs:**  
  - Console summaries of mean & SD tree circumference  
  - `circumference_start_end_boxplot.png` (start vs end boxplot)  
  - `mean_10yr_growth_by_site.csv` (mean growth summary)  
  - `ttest_10yr_growth_result.txt` (t-test results)  

### 3. `README.md`
This file (documentation).

## How to run
1. Clone the repository:
   ```bash
   git clone https://github.com/<USERNAME>/<REPO>.git
   ```
2. Open RStudio in the project folder.  
3. Run each script (`gene_expression_analysis.R`, `growth_analysis.R`).  
4. Outputs will be saved in the project directory.  

## Dependencies
- R (≥4.0)  
- Packages: `dplyr`, `tidyr`, `ggplot2`
