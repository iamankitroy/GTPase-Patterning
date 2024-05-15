#!/bin/zsh

# Track counter
counter()
{
	count=`sed '/^#/d' $1 | sed -n "/,$2$/p" | cut -d',' -f3 | sort | uniq | wc -l`
}

# Cleaner
cleanup()
{
	sed -i "s/ //g" $1
}

# Main function
main()
{
	echo "Filename,n_GTPase,n_GDI" | tee -a ncounts.csv

	for i in `ls *colocalisation.csv`
	do
		counter "$i" "GTPase"
		gtpase_count=$count			# GTPase track count

		counter "$i" "GDI"
		gdi_count=$count			# GDI track count

		echo "$i,$gtpase_count,$gdi_count" | tee -a ncounts.csv
	done

	cleanup ncounts.csv
}

# Run main
main

# Ankit Roy
# 10th May, 2024
