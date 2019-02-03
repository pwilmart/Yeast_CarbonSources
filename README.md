# Yeast_CarbonSources
Yeast grown with galactose, glucose, or raffinose carbon sources from the Gygi lab labeled with 10-plex TMT reagents, in triplicate.

These are analyses of a public dataset (PRIDE [PXD002875](https://www.ebi.ac.uk/pride/archive/projects/PXD002875)) from Paulo, O'Connell, Gaun, and Gygi:

> Paulo, J.A., O’Connell, J.D., Gaun, A. and Gygi, S.P., 2015. Proteome-wide quantitative multiplexed profiling of protein expression: carbon-source dependency in Saccharomyces cerevisiae. Molecular biology of the cell, 26(22), pp.4063-4074.

There were 24 RAW files of yeast grown in three different carbon sources. It was a 3x3 (9-plex) TMT experiment done with the SPS MS3 (MultiNotch) method on a Thermo Fusion instrument.

Analyses:
- [PAW pipeline](https://github.com/pwilmart/PAW_pipeline.git)
  - [CarbonSources_part-1](https://pwilmart.github.io/TMT_analysis_examples/CarbonSources_part-1.html)
  - CarbonSources_part-2
- MaxQuant
  - [CarbonSources_MQ](https://pwilmart.github.io/TMT_analysis_examples/CarbonSources_MQ.html)

---
## `PAW` folder contents

**File types:**
  - `*.ipynb` - Jupyter notebooks

  - `*.r` - code cells from notebooks

  - `*.html` - notebooks rendered in html

  - `results_files` folder:
    - `*.log` - console output log files from pipeline steps
    - `*.txt` - tab-delimited text results_files
      - protein summaries
      - peptide summaries
    - `*.xlsx` - Excel files
    - `R-input.txt` - prepped table of TMT data for importing into r
    - `CarbonSources_results.txt` - statistical testing results from r

---

## `MQ` folder contents

**Files:**
- CarbonSources_MQ.ipynb - Jupyter notebook for statistical analysis
- CarbonSources_MQ.html - html rendering of notebook
- CarbonSources_MQ.r - code cells from notebook
- CarbonSources_results.txt - statistical testing results
- parameters.txt - summary of MQ parameter settings
- proteinGroups.txt - main protein-level results file from MQ
- proteinGroups.xlsx - Modified Excel file (for table prepping)
- R-input.txt - prepped table for import into r
- summary.txt - summary file from MQ (LC run stats)

---

### R input table prep

Basic steps were similar for both pipelines:
1. flag proteins to exclude
  - common contaminants
  - decoys
  - proteins with no reporter ion signals
1. sort excluded proteins to bottom of table
1. make new tab
  - add column of protein accessions
  - add columns of the TMT channels
1. export the new tab contents to text files
  - table should be well-formed (single header line) and rectangular
1. read text file into R

Statistical test results in R are collected into a data frame in the same order as the imported proteins. At the end of the notebook, the results file is saved as a text file for adding back to the main protein results spreadsheet file.

Eventually, there needs to be a coherent, comprehensive summary file that contains the proteomics results, the statistical testing results, and any other information to aid biological interpretation (rich annotations, etc.). An Excel file is a good format for this since adding descriptive text and formatting are easy. A basic Excel sheet can be easily distributed in Supplemental files and opened in Open Office applications.
