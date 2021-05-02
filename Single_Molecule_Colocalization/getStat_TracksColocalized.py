#!/Users/roy/anaconda3/bin/python3

import argparse
import pandas as pd

__author__ = "Ankit Roy"
__copyright__ = "Copyright 2021, Bieling Lab, Max Planck Institute of Molecular Physiology"
__license__ = "GPL"
__maintainer__ = "Ankit Roy"
__status__ = "Development"

#--- Fetch arguments
def get_args():
    parser = argparse.ArgumentParser()

    # Required arguments group
    required_args = parser.add_argument_group(title = "Required arguments")

    # Spot colocalization file
    required_args.add_argument("-cf", "--colocalization_file",
                                help = "Spot colocalization file name.")
    
    args = parser.parse_args()

    return args

#--- Get data
def dataIN(filename):
	data = pd.read_csv(filename)
	return data

#--- Modify field
def modifyField(data, condition, name, value):
    data.loc[condition, name] = value
    return data

#--- Get colocalized tracks
def get_ColocalizedTracks(data):
    data = data.copy()
    data["COLOCALIZED_TRACK"] = False

    # group data with pseudo track IDs and channels
    grouping = data.groupby(["PSEUDO_TRACK_ID", "CHANNEL"])
    
    # identify colocalized tracks
    # any track with at least one colocalized spot is considered a colocalized track
    for gid, group in grouping:
        if group["COLOCALIZED_SPOT"].sum() > 0:
            # group condition
            condition = (data["PSEUDO_TRACK_ID"] == gid[0]) & (data["CHANNEL"] == gid[1])
            # modify group field "COLOCALIZED_TRACK" to True
            data = modifyField(data, condition, "COLOCALIZED_TRACK", True)

    return data

#--- Subset colocalized tracks
def subsetData(data):
    data = data.copy()
    data = data[data["COLOCALIZED_TRACK"]]
    return data

#--- Count the number of colocalized frames
def count_ColocalizedFrames(data):
    data = data.copy()
    data["TOTAL_FRAME_COUNT"] = 0
    data["COLOCALIZED_FRAME_COUNT"] = 0

    # group data with pseudo track IDs and channels
    grouping = data.groupby(["PSEUDO_TRACK_ID", "CHANNEL"])

    for gid, group in grouping:
        f_count = len(group)                        # total frame count
        cf_count = group["COLOCALIZED_SPOT"].sum()  # colocalized frame count

        # group condition
        condition = (data["PSEUDO_TRACK_ID"] == gid[0]) & (data["CHANNEL"] == gid[1])

        # add total frame count
        data = modifyField(data, condition, "TOTAL_FRAME_COUNT", f_count)
        # add colocalized frame count
        data = modifyField(data, condition, "COLOCALIZED_FRAME_COUNT", cf_count)

    # add free frame count
    data["FREE_FRAME_COUNT"] = data.apply(lambda x: x["TOTAL_FRAME_COUNT"] - x["COLOCALIZED_FRAME_COUNT"], axis=1)
    # add colocalized frame fraction
    data["COLOCALIZED_FRAME_FRACTION"] = data.apply(lambda x: x["COLOCALIZED_FRAME_COUNT"]/x["TOTAL_FRAME_COUNT"], axis=1)

    return data

#--- Main function
def main():
    global args

    args = get_args()                               # input arguments
    coloc_data = dataIN(args.colocalization_file)   # Colocalization data

    # Get colocalized tracks
    coloc_data = get_ColocalizedTracks(coloc_data)

    # Subset of colocalized tracks
    sub_coloc_data = subsetData(coloc_data)

    # Get number of colocalized frames
    sub_coloc_data = count_ColocalizedFrames(sub_coloc_data)
    
#--- Run main function
if __name__ == '__main__':
    main()

# Ankit Roy
# 29th April, 2021
