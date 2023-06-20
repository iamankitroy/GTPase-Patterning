#!/bin/zsh

# List image files
function list_files()
{
	ls *Lipid.tif > filenames.temp
}

# Run segmentation script
function segment_images()
{
	for filename in `cat filenames.temp`
	do
		./autoSegmentation_LipidPatch.py $filename
	done
}

# Get file prefixes
function get_file_prefixes()
{
	sed "s/C=.*//g" filenames.temp > prefixes.temp
}

# Analyze pattern statistics
function analyse_patterns()
{
	get_file_prefixes

	for prefix in `cat prefixes.temp`
	do
		./GTPase_patterning_analysis.py $prefix\C=Lipid_segmented.tif $prefix\C=Lipid.tif $prefix\C=Total.tif $prefix\C=Active.tif background_dict.p
		echo Analysed $prefix
	done
}

# Remove temporary files
function clean_up()
{
	rm *.temp
}

# Automate steps
function automate()
{
	list_files

	segment_images

	analyse_patterns

	clean_up
}

# Run automation
automate

# Ankit Roy
# 28th April, 2023
# 20th June, 2023			>>		Updated script to accomodate changes to the GTPase pattern analysis script.
#							>>		Added comments.
