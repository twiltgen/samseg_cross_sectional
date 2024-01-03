# samseg_cross_sectional
contains utilities and scripts for samseg cross sectional segmentation (with MS lesion segmentation)

### Requirements

Reuired to install the following dependencies:

- Freesurfer v. 7.3.2, c.f.: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall


### Introduction

The database should be BIDS conform, since scripts require a database as input that is in BIDS fromat.


### Processing

1. To run the samseg cross sectional pipeline, please install freesurfer and run the following command:

```
python3 xxx/run_samseg_cross.py -i /path/to/bids/database -n 15 --coregister --remove_temp
```

2. To aggregate all results into a single csv table for analysis please run the following command:

```
python3 xxx/run_analysis.py -i /path/to/bids/database -o /path/to/output/folder
```

### optional
You can check which cases have been processed by using the following command:

```
python3 xxx/check_processed.py -i /path/to/bids/database -o /path/to/output/folder
```
The script will generate two csv files: one with a list of processed images and one with a list of imaged that have npot been processed yet


### Any questions?

Please open an issue :)
