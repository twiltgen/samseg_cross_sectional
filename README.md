# samseg_cross_sectional
contains utilities and scripts for samseg cross sectional segmentation (with MS lesion segmentation)

### Requirements

Reuired to install the following dependencies:

- Freesurfer v. 7.3.2, c.f.: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall

I am using an environment with the packages that are listed in requirements.txt. 
Since this environment is used for other processing purposes as well, it might contain packages that are not required for the processing pipeline in this repository
You can install the packages by running the following command:
```
pip install -r requirements.txt
```


### Introduction

The database should be BIDS conform, since scripts require a database as input that is in BIDS format.


### Processing

1. To run the samseg cross sectional pipeline, please install freesurfer and run the following command:

```
python3 xxx/run_samseg_cross.py -i /path/to/bids/database -n 15 --coregister --remove_temp
```
Options:
- "-n 15": -n defines the number of workers that should be used 
(ATTENTION: in the script we set the number of threads used by SAMSEG to 4, i.e., 4x15 = 60) 

- "--coregister": use this flag if you want the T1w and FLAIR image to be coregistered before running SAMSEG, which is a requirement for SAMSEG. (The script also double checks if the coregistered files exist and calls mri_coreg if they are missing, so in theory the flag is not critical, but better safe than sorry :) )

- "--remove_temp": use this flag if you want to delete all files in the template folder. During processing, all output files of SAMSEG are stored in the same folder. In order to be BIDS conform, we need to rename and copy the files to their correct location (i.e., the "anat" folder in the derivatives). The original output files can then be kept or be deleted by using the --remove_temp flag.


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
