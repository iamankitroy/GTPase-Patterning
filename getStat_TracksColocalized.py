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
                                help = "Spot colocalization file name.",
                                required = True)

    # Filteration options
    parser.add_argument("--filter",
                                help = "Filter based on 1) FREE_FRAME_COUNT or 2) COLOCALIZED_FRAME_FRACTION",
                                choices = ['1', '2'],
                                default = '0')

    # Filter by free frames in GDI channel
    parser.add_argument("--free_frames",
                                help = "Filter by the maximum number of frames GDI spots can remain un-colocalized.",
                                default = 0)

    # Filter by colocalization fraction of GDI channel
    parser.add_argument("--colocalization_fraction",
                                help = "Filter by the minimum fraction of frames that GDI spots should remain be colocalized.",
                                default = 1)

    # Field of view cutoff
    parser.add_argument("-fov", "--field",
                                help = "(default = 1.00) Field of view",
                                default = 1.00,
                                type = float)

    # Pixel size
    parser.add_argument("-ps", "--pixel_size",
								help = "(default = 0.178 µm) Pixel size",
								default = 0.178,
								type = float)

    # Image size
    parser.add_argument("-is", "--image_size",
								help = "(default = 512 px) Image size",
								default = 512,
								type = int)
    
    args = parser.parse_args()

    return args

#--- Get data
def dataIN(filename):
	data = pd.read_csv(filename, comment="#")
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

#--- Channel based filter
def filter_channel(data, channel, pseudo_ids):
    data = data.copy()
    channel_filter = data["CHANNEL"] == channel                 # filter CHANNEL
    pseudo_id_filter = data["PSEUDO_TRACK_ID"].isin(pseudo_ids) # filter PSEUDO_TRACK_ID

    data = data.loc[channel_filter & pseudo_id_filter, ]        # combine filters

    return data

#--- Filter GDI tracks based on un-colocalized frames
def filter_freeFrames(data, count):
    data = data.copy()

    # COLOCALIZATION_ID of GDI tracks that remain un-colocalized for <= user defined number of frames
    coloc_ids = set(data.loc[(data["FREE_FRAME_COUNT"] <= count) & (data["CHANNEL"] == "GDI"), "COLOCALIZATION_ID"])
    
    # GTPase and GDI ids of desired COLOCALIZATION_ID
    gtpase_ids = [i.split('-')[0] for i in coloc_ids]
    gdi_ids = [i.split('-')[1] for i in coloc_ids]

    # Apply filters and combine
    gtpase_filtered = filter_channel(data, "GTPase", gtpase_ids)
    gdi_filtered = filter_channel(data, "GDI", gdi_ids)
    data_filtered = pd.concat([gtpase_filtered, gdi_filtered])
    
    return data_filtered

#--- Filter GDI tracks based on colocalization fraction
def filter_colocFrameFraction(data, fraction):
    data = data.copy()

    # COLOCALIZATION_ID of GDI tracks with colocalization fraction >= user defined fraction
    coloc_ids = set(data.loc[(data["COLOCALIZED_FRAME_FRACTION"] >= fraction) & (data["CHANNEL"] == "GDI"), "COLOCALIZATION_ID"])

    # GTPase and GDI ids of desired COLOCALIZATION_ID
    gtpase_ids = [i.split('-')[0] for i in coloc_ids]
    gdi_ids = [i.split('-')[1] for i in coloc_ids]

    # Apply filters and combine
    gtpase_filtered = filter_channel(data, "GTPase", gtpase_ids)
    gdi_filtered = filter_channel(data, "GDI", gdi_ids)
    data_filtered = pd.concat([gtpase_filtered, gdi_filtered])

    return data_filtered

#--- Calculate landing rate
def calcLandingRate(data):
    fov = args.field                # field of view
    image_size = args.image_size    # image size in pixels
    pixel_size = args.pixel_size    # pixel size in µm

    # Total number of frames
    total_frames = max(set(data["FRAME"]))
    
    # Number of landing events
    gtpase_landing = len(set(data.loc[data["CHANNEL"] == "GTPase", "PSEUDO_TRACK_ID"]))
    gdi_landing = len(set(data.loc[data["CHANNEL"] == "GDI", "PSEUDO_TRACK_ID"]))

    # Landing rate
    gtpase_landing_rate = gtpase_landing/(total_frames * image_size * pixel_size * fov)
    gdi_landing_rate = gdi_landing/(total_frames * image_size * pixel_size * fov)
    
    return (gtpase_landing_rate, gdi_landing_rate)

#--- All colocalization statistics:
def getStat(all_data, subset_data, input_file):
    # Count all GTPase tracks
    total_gtpase_tracks = len(set(all_data.loc[all_data["CHANNEL"] == "GTPase", "PSEUDO_TRACK_ID"]))
    # Count all GDI tracks
    total_gdi_tracks = len(set(all_data.loc[all_data["CHANNEL"] == "GDI", "PSEUDO_TRACK_ID"]))
    # Count colocalized GTPase tracks
    coloc_gtpase_tracks = len(set(subset_data.loc[subset_data["CHANNEL"] == "GTPase", "PSEUDO_TRACK_ID"]))
    # Count colocalized GDI tracks
    coloc_gdi_tracks = len(set(subset_data.loc[subset_data["CHANNEL"] == "GDI", "PSEUDO_TRACK_ID"]))
    # Percentage of colocalized GTPase tracks
    percent_coloc_gtpase = coloc_gtpase_tracks/total_gtpase_tracks * 100
    # Percentage of colocalized GDI tracks
    percent_coloc_gdi = coloc_gdi_tracks/total_gdi_tracks * 100

    # Calculate landing rate
    gtpase_landing_rate, gdi_landing_rate = calcLandingRate(all_data)

    # Display stats
    header = "# {},{},{},{},{},{},{},{},{}\n".format("File_Name", "GTPase_Count", "GDI_Count", "Colocalized_GTPase_Count", "Colocalized_GDI_Count", "Percent_Colocalized_GTPase", "Percent_Colocalized_GDI", "GTPase_Landing_Rate", "GDI_Landing_Rate")
    stat_line = "# {},{},{},{},{},{:.2f},{:.2f},{:.2E},{:.2E}\n".format(input_file, total_gtpase_tracks, total_gdi_tracks, coloc_gtpase_tracks, coloc_gdi_tracks, percent_coloc_gtpase, percent_coloc_gdi, gtpase_landing_rate, gdi_landing_rate)

    return (header, stat_line)

#--- Write output files
def dataOUT(data_frame, outname, header, stat_line):
    outname = "{}_subset.csv".format(outname.split('.csv')[0])
    
    # write header and statistics
    with open(outname, 'w') as fh:
        fh.write("# {:=^40}\n".format(" Meta-data lines "))
        for arg in vars(args):
            fh.write("# {}: {}\n".format(arg, getattr(args, arg)))
        fh.write("# {:=^40}\n".format(" Summary lines "))
        fh.write(header)
        fh.write(stat_line)
        fh.write("# {:=^40}\n".format(" Colocalization data lines "))

    # write CSV file for colocalization events
    data_frame.to_csv(outname, index=False, float_format="%.3f", mode = 'a')

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

    # This is still experimental!!    
    if int(args.filter) == 1:
        # Filter GDI tracks based on the number of uncolocalized frames
        sub_coloc_data = filter_freeFrames(sub_coloc_data, args.free_frames)
    elif int(args.filter) == 2:
        # Filter GDI tracks based on fraction of colocalized frames
        sub_coloc_data = filter_colocFrameFraction(sub_coloc_data, args.colocalization_fraction)

    # Show colocalization statistics
    header, stat_line = getStat(coloc_data, sub_coloc_data, args.colocalization_file)

    # Wtite output file
    dataOUT(sub_coloc_data, args.colocalization_file, header, stat_line)

#--- Run main function
if __name__ == '__main__':
    main()

# Ankit Roy
# 29th April, 2021
# 3rd May, 2021
#   --> Added ability to filter based on the number of free frames in GDI tracks
#   --> Added ability to filter based on the fraction of colocalized frames in GDI tracks
#   --> Displays track based colocalization statistics
# 6th April, 2021
#   --> Now calculates landing rate in terms of tracks per µm^2
#   --> Stores colocalization statistics along with data in output CSV file
# 9th August, 2021
#   --> Colocalization file argurent set to a required argument
#   --> Ignores comment lines from colocalization file
#   --> Now adds input arguments as meta-data to output files.