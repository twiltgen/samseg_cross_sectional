import argparse
import os
from pathlib import Path
import multiprocessing
from tqdm import tqdm

def split_list(alist, splits=1):
    length = len(alist)
    return [alist[i * length // splits: (i + 1) * length // splits]
            for i in range(splits)]

parser = argparse.ArgumentParser(description='Segment t1/f2 of brains in given database.')
parser.add_argument('-i', '--input_directory', help='Folder of database.', required=True)
parser.add_argument('-n', '--number_of_workers', help='Number of parallel processing cores.', type=int, default=os.cpu_count())

# read the arguments
args = parser.parse_args()
print("Number of CPU Processes:",args.number_of_workers)

# make the derivatives directory if it is not available
from pathlib import Path
Path(os.path.join(args.input_directory,"derivatives")).mkdir(parents=True, exist_ok=True)

derivatives_dir = os.path.join(args.input_directory,"derivatives", "labels", "samseg")

# TO Do - incorporate the bids validator!

t1w_label = 'T1w.nii.gz'
t2w_label = 'T2w.nii.gz'
flair_label = 'FLAIR.nii.gz'
co_reg_flair = 'reg-FLAIR.lta'
co_reg_flair_nifti = 'reg-FLAIR.nii'

# define the common convention for t1.nii / t2.nii / etc here
# @TODO change this to the BIDS file-format structure in the long run
# assuming the old file tree structure that was used in the lab from 2017-2022


# gather the whole directory list and split it among the cpu cores
# instead of looping over patients only that have a different number of scans!

# required files are t1 and f2
# all folders including a t1 brain scan
t1_list = [str(path.parent) for path in Path(args.input_directory).rglob(f'*{t1w_label}')]
# get unique folders
t1_set = set(t1_list)
# similarly for f2
f2_list = [str(path.parent) for path in Path(args.input_directory).rglob(f'*{flair_label}')]
f2_set = set(f2_list)

# only use sessions that have both, T1w and Flair!
paths = list(set.intersection(t1_set,f2_set))
files = split_list(paths, args.number_of_workers)

print(files)

def getSubjectID(path):
    stringList = path.split("/")
    indices = [i for i, s in enumerate(stringList) if 'sub-' in s]
    return stringList[indices[0]]

def getSessionID(path):
    stringList = path.split("/")
    indices = [i for i, s in enumerate(stringList) if 'ses-' in s]
    return stringList[indices[0]]

def run_samseg(file_list,number):

    print("Running process:",number)
    for session_path in tqdm(file_list):

        # extract patient ids
        patient_id = getSubjectID(session_path)
        # extract session_ids
        session_id = getSessionID(session_path)
        T1w_bids_path = os.path.join(session_path, f'{patient_id}_{session_id}_{t1w_label}')
        print(T1w_bids_path)
        FLAIR_bids_path = os.path.join(session_path, f'{patient_id}_{session_id}_{flair_label}')
        print(FLAIR_bids_path)

        # make the directories for the derivative folders

        CoReg_Flair_path_native = os.path.join(derivatives_dir,patient_id, session_id, f'{patient_id}_{session_id}_{co_reg_flair}')
        CoReg_Flair_path_nifti= os.path.join(derivatives_dir,patient_id, session_id,f'{patient_id}_{session_id}_{co_reg_flair_nifti}')

        derivatives_path = os.path.join(derivatives_dir,patient_id, session_id)


        try:
            # assert if scans exist
            assert os.path.isfile(T1w_bids_path), 'T1 not available'
            assert os.path.isfile(FLAIR_bids_path), 'F2 not available'

            print("running it")

            # do samseg processing
            os.system(f'export FREESURFER_HOME=$HOME/freesurfer ; \
                        mri_coreg --mov {FLAIR_bids_path} --ref {T1w_bids_path} --reg {CoReg_Flair_path_native} ; \
                        mri_vol2vol --mov {FLAIR_bids_path} --reg {CoReg_Flair_path_native} --o {CoReg_Flair_path_nifti} --targ {T1w_bids_path} ; \
                        run_samseg --input {T1w_bids_path} {CoReg_Flair_path_nifti} --pallidum-separate --output {derivatives_path} ;\
                        ')

        except AssertionError:
            print('Sequences not available. Proceeding to next patient or patient scan.')
            continue

# initialize multithreading
pool = multiprocessing.Pool(processes=args.number_of_workers)

# creation, initialisation and launch of the different processes
for x in range(0, args.number_of_workers):
    pool.apply_async(run_samseg, args=(files[x], x))

pool.close()
pool.join()
