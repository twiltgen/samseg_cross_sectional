import argparse
import os
import pandas as pd

from utils import getfileList, getSessionID, getSubjectID

def combineStats(path, subID, sesID):
    '''
    This function reads the _samseg.stats file and the _sbtiv.stats file of a session and merges both into a dataframe

    :param path: path to BIDS derviatives database
    :param subID: the subject ID of the session
    :param sesID: the session ID of the session
    :return: return the volumes in a dataframe
    '''
    # define path of stat files and load them as dataframe
    samseg_path = os.path.join(path, "sub-"+subID, "ses-"+sesID, "anat", "sub-"+subID+"_ses-"+sesID+"_space-T1w_samseg.stats")
    tiv_path = os.path.join(path, "sub-"+subID, "ses-"+sesID, "anat", "sub-"+subID+"_ses-"+sesID+"_space-T1w_sbtiv.stats")
    df_samseg_stat = pd.read_csv(samseg_path, header=None, names=["ROI", "volume", "unit"])
    df_tiv_stat = pd.read_csv(tiv_path, header=None, names=["ROI", "volume", "unit"])

    # combine _samseg and _sbtiv and clean ROI names
    df = pd.concat([df_samseg_stat, df_tiv_stat])
    df["ROI"]=df["ROI"].str.replace("# Measure ","")
    # transpose dataframe and clean indices and column names
    df = df.loc[:,["ROI", "volume"]].reset_index().drop("index", axis=1)
    df = df.transpose()
    df.columns = list(df.iloc[0,0:])
    df = df.drop(index = 'ROI').reset_index().drop("index", axis=1)
    # add subject ID and session ID to the dataframe
    df["sub-ID"] = subID
    df["ses-ID"]= sesID
    df_IDs=df[["sub-ID", "ses-ID"]]
    df.drop(labels=["sub-ID", "ses-ID"], axis=1, inplace=True)
    df.insert(0, "ses-ID", df_IDs["ses-ID"])
    df.insert(0, "sub-ID", df_IDs["sub-ID"])

    return df

####################################################
# main script

parser = argparse.ArgumentParser(description='Read Volumes of SAMSEG Longitudinal Segmentation.')
parser.add_argument('-i', '--input_directory', help='Folder of derivatives in BIDS database.', required=True)
parser.add_argument('-o', '--output_directory', help='Destination folder for the output table with volume stats.', default='/home/twiltgen/media/twiltgen/raid3/Tun/MR_database/Data/Database')

# read the arguments
args = parser.parse_args()

# define path of the derivatives folder
derivatives_dir = os.path.join(args.input_directory, "derivatives/samseg-7.3.2")

# get a list with the paths of the _seg.mgz files 
# (we do this because we only want subject- and session-IDs for the cases that have been successfully segmented by SAMSEG)
seg_list = getfileList(path = derivatives_dir, 
                       suffix = '*_seg.mgz')
# initialize empty dataframe into which we will write the stats data of all cases
df_stat = pd.DataFrame()
for i in range(len(seg_list)):
    # get subject and session ID
    subjectID = getSubjectID(seg_list[i])
    sessionID = getSessionID(seg_list[i])
    # get stats of current session
    loop_stat = combineStats(derivatives_dir, subjectID, sessionID)
    # write stats in the final dataframe
    df_stat = pd.concat([df_stat, loop_stat])
    print(f'sub-{subjectID}_ses-{sessionID}: stats added.')

# write stats table to .csv file in chosen output directory
df_stat.to_csv(os.path.join(args.output_directory, "volume_samseg_raw.csv"), index=False)
