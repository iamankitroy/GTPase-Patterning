#!/bin/bash

# This script renames all CSV files in the working directory
# Files are renamed such that all whitespaces are replaced by underscores

IFS=$'\n'

for file in `ls *.csv`
do
	better_name=`echo $file | sed 's/ /_/g'`
	mv $file $better_name
done

echo "Files renamed."

# Ankit Roy
# 12th August, 2021
