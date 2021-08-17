#!/usr/local/bin/zsh

# !! Dependecies !!
# The script makes use of associative arrays. Please make sure you use Bash 4.0 or higher OR Zsh to run the script.

# This script is used to batch process the colocalization analysis

#--- Read arguments
function read_arguments()
{
	param_file="input_parameters.txt"		# input parameters file
	args=`head -n 1 $param_file`			# input arguments
}

#--- Rename files
function rename_files()
{
	IFS=$'\n'

	for filename in `ls *.csv`
	do
		better_name=`echo $filename | tr ' ' '_'`	# replace whitespaces with _
		mv $filename $better_name					# rename file
	done
}

#--- Clean parameter file
function clean_parameter_file()
{
	# remove trailing spaces in parameter fields
	# replace internal whitespaces with _
	sed 's/ ,/,/g' $param_file | sed 's/, /,/g' | sed 's/ /_/g' > temp.txt
	mv temp.txt $param_file
}

#--- Run SpotColocalization
function run_SpotColocalization()
{
	echo -e "\nRunning SpotColocalization...\n"
	./SpotColocalization.py -d ${arg_list[-d]} -fov ${arg_list[-fov]} -ps ${arg_list[-ps]} -is ${arg_list[-is]} --first_frame ${arg_list[--first_frame]} --last_frame ${arg_list[--last_frame]} -gp ${arg_list[-gp]} -gd ${arg_list[-gd]} --outfile ${arg_list[--outfile]}
	echo -e "\nComplete!\n"
}

#--- Run getStat_TrackColocalized
function run_getStat_TrackColocalized()
{
	echo -e "\nRunning getStat_TrackColocalized...\n"
	python getStat_TracksColocalized.py -fov ${arg_list[-fov]} -ps ${arg_list[-ps]} -is ${arg_list[-is]} -cf ${arg_list[--outfile]}
	echo -e "\nComplete!\n"
}

#--- Set default arguments
function set_default_arguments()
{
	arg_list[-d]=0.5
	arg_list[-fov]=1.00
	arg_list[-ps]=0.178
	arg_list[-is]=512
	arg_list[--first_frame]=1
	arg_list[--last_frame]=1000
	arg_list[--outfile]="Colocalization.csv"
}

#--- Construct argument list and run sub process
function run_processes()
{
	args=(`echo $args | tr ',' '\n'`)				# argument name array
	all_params=(`tail -n+2 $param_file`)			# array of all parameters

	total_subjobs=${#all_params[@]}

	# combine parameter with argument name
	for param_line in ${all_params[@]}
	do
		params=(`echo $param_line | tr ',' '\n'`)	# parameters
		declare -A arg_list							# argument list
		set_default_arguments						# setting defaults

		for field_num in $(seq 0 ${#args[@]})
		do
			key=${args[$field_num]}					# argument name
			value=${params[$field_num]}				# argument value
			arg_list[$key]=$value					# store argument name and value pair
		done

		# run SpotColocalization
		run_SpotColocalization

		# run getStat_TrackColocalized
		run_getStat_TrackColocalized

		# progress status
		total_subjobs=$(( $total_subjobs - 1 ))
		echo -e "\n========== $total_subjobs REMAINING! ==========\n\n"

	done
}

#--- Main function
function automate()
{
	read_arguments
	echo "Arguments read."

	rename_files
	echo "Files renamed."

	clean_parameter_file
	echo "Parameter file cleaned."

	echo "\n========== Batch processing started! =========="
	run_processes
	echo "Batch processing complete!"
}

#--- Automate batch processing
automate

# Ankit Roy
# 13th August, 2021
# 16th August, 2021
#	-->	Updated script to use associative arrays to parse argparse arguments
