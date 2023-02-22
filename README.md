# Cellular energy regulates mRNA translation and degradation in a codon-specific manner 
Code for the data processing and analyses present on the manuscript: "Cellular energy regulates mRNA translation and degradation in a codon-specific manner". 

## Repository Structure and Description

This repository is divided into figure generating scripts (found in `figures`) and analysis/data processing scripts (found in `scripts`).

* The `scripts` folder is divided into data processing and analysis. The analysis contains notebooks named with the numbers of the figures they belong to in the paper. The data processing is divided into 5PSeq, GTEx and Mouse (Tabula Muris). Each folder provides pipelines to process and generate data required for the respective analyses.

* The `figures` folder contains the scripts to generate paper figures. All data required to generate them is created from the analysis notebooks, which store it in `figures/figure_data`.


## Data availability
* The GTEx data used to calculate the exonic-to-intronic read count ratio and obtain certain sample annotations is protected. It can only be accessed with a suitable license. Therefore, the analysis and data processing notebooks deriving from such data do not display their outputs and only display their code.

* Publicly shareable data required to execute the analysis can be found here: https://zenodo.org/record/7662464.



