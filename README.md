# samseg_cross_sectional
contains utilities and scripts for samseg cross sectional segmentation (with MS lesion segmentation)

### Requirements

Reuired to install the following dependencies:

- Freesurfer v. 7.3.2, c.f.: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall


Please create a conda environment with the following packages:
- 

### Introduction

The database should be BIDS conform, since scripts require a database as input that is in BIDS fromat.

### Dataset structure
BIDS-compliant datastructure:
Output Data Structure (not 100% BIDS, we are working on it!). Important files highlighted!:

### Processing

1. To run the samseg longitudinal pipeline + fsl-based pbvc calculation, please install freesurfer and run the following command:

```
python3 xxx/xxx.py --input_directory /path/to/bids --number_of_workers 32 --freesurfer_path /path/to/fs/installation
```

2. To aggregate all results into a single csv tabel for analysis please run the following command:

```
python3 xxx/run_analysis.py --input_directory /path/to/processed/cohort --output_directory /path/to/output
```

### Any questions?

Please open an issue :)
