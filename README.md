# Code to reproduce analysis published in:
"A global survey of intramolecular isopeptide bonds"

**Preprint available at:** 
[]()

This repository contains the scripts, notebooks, and data structure used to perform the analyses reported in the manuscript. It is intended to enable reproducibility and transparency.

## Requirements
The analysis was performed using Python 3.12.2 with the following packages:
```bash
biotite==1.0.1  
pandas==2.1.1  
seaborn==0.13.2  
bio==1.83  
matplotlib==3.8.0  
numpy==1.26.4
```

## Project structure

```bash
.
├── analysis/        # Subprojects with inputs, outputs, scripts and notebooks
├── bin/             # Shared Python scripts
├── data/            # Input datasets (see data/tables/README.md)
├── notebooks/       # Jupyter notebooks for data analysis
├── .env             # Project-wide configuration variables
├── figures          # Programatically generated figures
├── README.md        # This file

```

Each subproject within analysis/ follows a consistent structure:

- input/, output/, and bin/ or sh/ folders

- A README.md with subproject-specific notes

- Optional Jupyter notebooks

Naming convention for files:

```bash
<YYYYMMDD>_<name>_<version>
```

## Reproducing the Analysis

To reproduce the key analyses:

1. Clone this repository.

2. Install the required packages (see above).

3. Download large files from [Zenodo](https://doi.org/10.5281/zenodo.15024939)

4. Follow the instructions in each subproject’s README.md or use the Jupyter notebooks provided in notebooks/ or analysis/*.

## Data Availability

All non-generated input data is stored in data/ and described in data/README.md.
Output files are reproducible and not backed up here but downloadable from [Zenodo](https://doi.org/10.5281/zenodo.15024939).


## Contact

For questions related to the analysis, contact:
Francesco Costa
Email: fcosta@ebi.ac.uk
