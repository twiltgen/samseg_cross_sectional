import argparse
import os
#from tqdm import tqdm
from pathlib import Path
import multiprocessing

def split_list(alist, splits=1):
    length = len(alist)
    return [alist[i * length // splits: (i + 1) * length // splits]
            for i in range(splits)]

parser = argparse.ArgumentParser(description='Segment t1/f2 of brains in given database.')
parser.add_argument('-i', '--input_directory', help='Folder of database.', required=True)
parser.add_argument('-n', '--number_of_workers', help='Number of parallel processing cores.', type=int, default=os.cpu_count())

t1w_label = 't1.nii'
t2w_label = 't2.nii'
flair_label = 'f2.nii'
co_reg_flair = 'f2_to_t1.lta'
co_reg_flair_nifti = 'f2_reg.nii'

# define the common convention for t1.nii / t2.nii / etc here
# @TODO change this to the BIDS file-format structure in the long run
# assuming the old file tree structure that was used in the lab from 2017-2022

# read the arguments
args = parser.parse_args()
print("Number of CPU Processes:",args.number_of_workers)

# gather the whole directory list and split it among the cpu cores
# instead of looping over patients only that have a different number of scans!

# required files are t1 and f2
# all folders including a t1 brain scan
t1_list = [str(path.parent) for path in Path(args.input_directory).rglob('*t1.nii')]
# get unique folders
t1_set = set(t1_list)
# similarly for f2
f2_list = [str(path.parent) for path in Path(args.input_directory).rglob('*f2.nii')]
f2_set = set(f2_list)

paths = list(set.intersection(t1_set,f2_set))
files = split_list(paths, args.number_of_workers)

print(files)


def run_samseg(file_list,number):

    print("Running process:",number)
    for session_path in tqdm(file_list):
        t1 = os.path.join(session_path, t1w_label)
        f2 = os.path.join(session_path, flair_label)

        try:
            # assert if scans exist
            assert os.path.isfile(t1), 'T1 not available'
            assert os.path.isfile(f2), 'F2 not available'

            print("running it")

            # do samseg processing
            os.system(f'export FREESURFER_HOME=$HOME/freesurfer ; \
                        cd {session_path} ; \
                        mri_coreg --mov {flair_label} --ref {t1w_label} --reg {co_reg_flair} ; \
                        mri_vol2vol --mov {flair_label} --reg {co_reg_flair} --o {co_reg_flair_nifti} --targ {t1w_label} ; \
                        run_samseg --input {t1w_label} {co_reg_flair_nifti} --pallidum-separate --output samsegAtlasOutput/ ; \
                        run_samseg --input {t1w_label} {co_reg_flair_nifti} --pallidum-separate --lesion --lesion-mask-pattern 0 1 --output samsegLesionOutput/ ; \
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
