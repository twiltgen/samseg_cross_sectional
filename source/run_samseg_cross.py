import argparse
import os
import shutil
import datetime
from pathlib import Path
import multiprocessing
from utils import getSessionID, getSubjectID, CopyandCheck, split_list, getfileList

def coreg_T1_FLAIR(derivatives_dir, im_t1, im_flair, im_flair_reg, output_dir, freesurfer_path):
    """
    This function coregisters the T1w and FLAIR images (FLAIR -> T1w), which is necessary for proceesing with samseg. 
    We use the mri_coreg to generate the transformation and mri_vol2vol to apply the transformation to the FLAIR image.
    All resulting files are saved to the output folder. In order to be compliant with BIDS convention, we copy these files from the output folder 
    to the corresponding location in the BIDS databse. 

    Parameters:
    -----------
    derivatives_dir : str
        Path of the SAMSEG derivatives folder in the BIDS database
    im_t1 : str
        Path of the T1w image in original subject space
    im_flair : str
        Path of the FLAIR image in original subject space
    output_dir : str
        Path to the output folder where the SAMSEG files should be stored
    freesurfer_path : str
        Path to freesurfer binaries
    
    Returns:
    --------
    None 
        This function produces FLAIR images registered to T1w images. Resulting files are copied to BIDS database.
    """
    # get subject and session IDs
    subID = getSubjectID(path = im_t1)
    sesID = getSessionID(path = im_t1)


    # pre-define paths of registered image(s) 
    flair_reg_field = str(Path(im_flair_reg).name).replace(".mgz", ".lta")


    # create output folder and BIDS target folder if they do not exist
    if not os.path.exists(output_dir):
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    deriv_ses = os.path.join(derivatives_dir, f'sub-{subID}', f'ses-{sesID}', 'anat')
    if not os.path.exists(deriv_ses):
        Path(deriv_ses).mkdir(parents=True, exist_ok=True)


    # run mri_coreg and get transformation
    print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: start FLAIR->T1w registration...')
    os.system(f'export FREESURFER_HOME={freesurfer_path} ; \
                cd {output_dir}; \
                mri_coreg --mov {im_flair} --ref {im_t1} --reg {flair_reg_field};\
                ')
    # run mri_vol2vol and apply transformation to FLAIR
    os.system(f'export FREESURFER_HOME={freesurfer_path} ; \
                cd {output_dir}; \
                mri_vol2vol --mov {im_flair} --reg {flair_reg_field} --o {im_flair_reg} --targ {im_t1};\
                ')
    print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: FLAIR->T1w registration DONE!')


    # copy the FLAIR transformation file
    print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: copy SAMSEG output files to BIDS...')
    reg_field_temp_location = os.path.join(output_dir, flair_reg_field)
    reg_field_target_location = os.path.join(deriv_ses, flair_reg_field)
    CopyandCheck(src = reg_field_temp_location, 
                 dst = reg_field_target_location)
    # copy the registered FLAIR image file
    flair_temp_location = os.path.join(output_dir, im_flair_reg)
    flair_target_location = os.path.join(deriv_ses, im_flair_reg)
    CopyandCheck(src = flair_temp_location, 
                 dst = flair_target_location)

def process_samseg(dirs, derivatives_dir, freesurfer_path, remove_temp=False, coregister=False):
    """
    This function applies SAMSEG segmentation and also applies required pre-processing steps of the T1w and FLAIR images if necessary. 
    Pre-processing includes mri_coreg to generate the transformation and mri_vol2vol to apply the transformation to the FLAIR image.
    Next, run SAMSEG cross sectional segmentation (including lesion segmentation). We use the T1w image and the registered FLAIR image as input.
    Next, we check if SAMSEG segmentation was successful by making sure that the seg.mgz file was generated.
    All resulting files are saved to the output folder. In order to be compliant with BIDS convention, we copy these files from the output folder 
    to the corresponding location in the BIDS databse. 
    Optionally, the output folder can be deleted (e.g., to clean up if it is not needed anymore)

    Parameters:
    -----------
    dirs : list
        List with all subject IDs for which we want to generate SAMSEG segmentations
    derivatives_dir : str
        Path of the SAMSEG derivatives folder in the BIDS database
    freesurfer_path : str
        Path to freesurfer binaries
    remove_temp : bool
        Boolean variable indicating if the output folder should be removed after segmentation files were generated and copied to BIDS database
    coregister : bool
        Boolean variable indicating if the FLAIR image needs to be registered to the T1w image
    
    Returns:
    --------
    None 
        This function produces SAMSEG segmentation files, inlcuding brain volume segmentation files (with lesion segmentation) 
    """
    # iterate through all subject folders
    for dir in dirs:

        # assemble T1w file lists
        # since we need T1w AND FLAIR images, we can first list all T1w images and then check if FLAIR image also exists
        t1w = getfileList(path = dir, 
                          suffix = '*T1w.*')
        t1w = [str(x) for x in t1w if (('.nii.gz' in str(x)) and 
                                       (not 'GADOLINIUM' in str(x)) and
                                       (not 'derivatives' in str(x)))]
        
        if len(t1w)>0:
            # get subject ID of current subject
            subID = getSubjectID(path = t1w[0])
        else:
            print(f'{datetime.datetime.now()} sub-{subID}: No T1w image available, proceed to next case...')
            continue


        # iterate over all session with T1w images and check if FLAIR files are available and if segmentation already exists
        for i in range(len(t1w)):
            try:

                # get session ID of current T1w image
                sesID = getSessionID(path = t1w[i])

                # check availability of files and folders (create folders if necessary)
                flair = str(t1w[i]).replace('_T1w.nii.gz', '_FLAIR.nii.gz')
                if not os.path.exists(flair):
                    raise ValueError(f'sub-{subID}_ses-{sesID}: FLAIR image not available!!')
                
                temp_dir = os.path.join(derivatives_dir, f'sub-{subID}', f'ses-{sesID}', 'temp')
                if not os.path.exists(temp_dir):
                    Path(temp_dir).mkdir(parents=True, exist_ok=True)
                
                deriv_ses = os.path.join(derivatives_dir, f'sub-{subID}', f'ses-{sesID}', 'anat')
                if not os.path.exists(deriv_ses):
                    Path(deriv_ses).mkdir(parents=True, exist_ok=True)
                
                # skip to next case if segmentation and samseg stats already exist
                seg_file = os.path.join(deriv_ses, f'sub-{subID}_ses-{sesID}_space-T1w_seg.mgz')
                stats_file = os.path.join(deriv_ses, f'sub-{subID}_ses-{sesID}_space-T1w_samseg.stats')
                tiv_file = os.path.join(deriv_ses, f'sub-{subID}_ses-{sesID}_space-T1w_sbtiv.stats')
                if os.path.exists(seg_file) and os.path.exists(stats_file) and os.path.exists(tiv_file):
                    print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: SAMSEG segmentation already exists, skip and proceed to next case...')
                    continue

                # register FLAIR->T1w if --coregister flag was used or if registered FLAIR does not exist
                flair_reg = os.path.join(deriv_ses, str(Path(flair).name).replace("FLAIR.nii.gz", "space-T1w_FLAIR.mgz"))
                if coregister or (not os.path.exists(flair_reg)):
                    # register FLAIR->T1w and copy files to target location in BIDS databse
                    coreg_T1_FLAIR(derivatives_dir = derivatives_dir,
                                   im_t1 = t1w[i],
                                   im_flair = flair,
                                   im_flair_reg = Path(flair_reg).name,
                                   output_dir = temp_dir,
                                   freesurfer_path = freesurfer_path)


                # run SAMSEG cross sectional segmentation 
                print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: start SAMSEG segmentation...')
                os.system(f'export FREESURFER_HOME={freesurfer_path} ; \
                            cd {temp_dir}; \
                            run_samseg --input {t1w[i]} {flair_reg} --threads 4 --pallidum-separate --lesion --lesion-mask-pattern 0 1 -o .\
                            ')
                print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: SAMSEG segmentation DONE!')


                # check if folder contains seg.mgz file, indicating that SAMSEG successfully finished
                output_files = os.listdir(temp_dir)
                output_files = [str(x) for x in output_files if ('_space-T1w_FLAIR.' not in str(x))]
                if ('seg.mgz' in output_files):
                    print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: copy SAMSEG files to BIDS database...')
                    # iterate over all output files and copy them to derivatives anat folder
                    for filename in output_files:
                        # rename to BIDS
                        file_bids = f'sub-{subID}_ses-{sesID}_space-T1w_' + filename
                        # define location of file in temporary fodler
                        filename_temp_location = os.path.join(temp_dir, filename)
                        # define target location of file
                        filename_target_location = os.path.join(deriv_ses, file_bids)
                        # copy file
                        CopyandCheck(src = filename_temp_location, 
                                     dst = filename_target_location)
                    print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: copied all SAMSEG files to BIDS database!')
                else:
                    print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: failed to generate segmentation, delete ouput folder...')
                    shutil.rmtree(temp_dir)
                    if os.path.exists(temp_dir):
                        raise ValueError(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: failed to delete the template folder: {temp_dir}')
                    else:
                        print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: successfully deleted the template folder: {temp_dir}')
                        continue
                
                
                # convert .mgz FLAIR and segmentation files to nifti
                print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: convert seg.mgz to nifti file...')
                seg_file_nii = str(seg_file).replace('.mgz', '.nii.gz')
                os.system(f'export FREESURFER_HOME={freesurfer_path} ; \
                            mri_convert {seg_file} {seg_file_nii}')
                print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: convert seg.mgz to nifti DONE!')
                print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: convert space-T1w_FLAIR.mgz to nifti file...')
                flair_reg_nii = str(flair_reg).replace('.mgz', '.nii.gz')
                os.system(f'export FREESURFER_HOME={freesurfer_path} ; \
                            mri_convert {flair_reg} {flair_reg_nii}')
                print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: convert space-T1w_FLAIR.mgz to nifti DONE!')


                # delete the temp folder if --remove_temp flag was used
                if remove_temp:
                    shutil.rmtree(temp_dir)
                    if os.path.exists(temp_dir):
                        raise ValueError(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: failed to delete the template folder.')
                    else:
                        print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: successfully deleted the template folder.')

            except:
                print(f'{datetime.datetime.now()} sub-{subID}_ses-{sesID}: Error occured during processing, proceeding with next case.')
            
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run SAMSEG cross sectional Pipeline on cohort.')
    parser.add_argument('-i', '--input_directory', help='Folder of derivatives in BIDS database.', required=True)
    parser.add_argument('-n', '--number_of_workers', help='Number of parallel processing cores.', type=int, default=os.cpu_count()-1)
    parser.add_argument('-f', '--freesurfer_path', help='Path to freesurfer binaries.', default='/home/twiltgen/Tun_software/Freesurfer/FS_7.3.2/freesurfer')
    parser.add_argument('--coregister', action='store_true')
    parser.add_argument('--remove_temp', action='store_true')

    # read the arguments
    args = parser.parse_args()

    if args.remove_temp:
        remove_temp = True
    else:
        remove_temp = False
    
    if args.coregister:
        coregister = True
    else:
        coregister = False


    # generate derivatives/labels/
    derivatives_dir = os.path.join(args.input_directory, "derivatives/samseg-7.3.2")
    if not os.path.exists(derivatives_dir):
        Path(derivatives_dir).mkdir(parents=True, exist_ok=True)

    
    # generate list with subject folders for multiprocessing
    data_root = Path(os.path.join(args.input_directory))
    dirs = sorted(list(data_root.glob('*')))
    dirs = [str(x) for x in dirs]
    dirs = [x for x in dirs if "sub-" in x]
    files = split_list(alist = dirs, 
                       splits = args.number_of_workers)


    # initialize multithreading
    pool = multiprocessing.Pool(processes=args.number_of_workers)
    # call samseg processing function in multiprocessing setting
    for x in range(0, args.number_of_workers):
        pool.apply_async(process_samseg, args=(files[x], derivatives_dir, args.freesurfer_path, remove_temp, coregister))

    pool.close()
    pool.join()

