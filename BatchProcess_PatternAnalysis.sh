#!/bin/zsh

function list_files()
{
	ls *Lipid.tif > filenames.temp
}

function segment_images()
{
	for filename in `cat filenames.temp`
	do
		./autoSegmentation_LipidPatch.py $filename
	done
}

function get_file_prefixes()
{
	sed "s/C=.*//g" filenames.temp > prefixes.temp
}

function analyse_patterns()
{
	get_file_prefixes

	for prefix in `cat prefixes.temp`
	do
		./GTPase_patterning_analysis.py $prefix\C=Lipid_segmented.tif $prefix\C=Lipid.tif $prefix\C=Total.tif $prefix\C=Active.tif
		echo Analysed $prefix
	done
}

function clean_up()
{
	rm *.temp
}

function automate()
{
	list_files

	segment_images

	analyse_patterns

	clean_up
}

automate
