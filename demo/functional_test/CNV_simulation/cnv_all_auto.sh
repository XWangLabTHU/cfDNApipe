#!/bin/bash

mkdir ./IC15/CNV/Bin_auto
mkdir ./IC17/CNV/Bin_auto
mkdir ./IC20/CNV/Bin_auto
mkdir ./IC35/CNV/Bin_auto
mkdir ./IC37/CNV/Bin_auto

arr=("01" "05" "1" "2" "3" "4" "5" "6" "7" "8" "9")

for ratio in ${arr[@]}
do
	echo $ratio
	for i in $(seq 1 1 10)
	do
		mkdir "./IC15/CNV/Bin_auto/ds_"$ratio"_"$i
		mkdir "./IC17/CNV/Bin_auto/ds_"$ratio"_"$i
		mkdir "./IC20/CNV/Bin_auto/ds_"$ratio"_"$i
		mkdir "./IC35/CNV/Bin_auto/ds_"$ratio"_"$i
		mkdir "./IC37/CNV/Bin_auto/ds_"$ratio"_"$i
		nohup python WGS_hg19_pe_without_control_auto.py "./IC15/CNV/Bin_auto/ds_"$ratio"_"$i "./IC15/IC15_"$ratio"_"$i".bam" > "./IC15/CNV/Bin_auto/ds_"$ratio"_"$i".log" 2>&1 &
		nohup python WGS_hg19_pe_without_control_auto.py "./IC17/CNV/Bin_auto/ds_"$ratio"_"$i "./IC17/IC17_"$ratio"_"$i".bam" > "./IC17/CNV/Bin_auto/ds_"$ratio"_"$i".log" 2>&1 &
		nohup python WGS_hg19_pe_without_control_auto.py "./IC20/CNV/Bin_auto/ds_"$ratio"_"$i "./IC20/IC20_"$ratio"_"$i".bam" > "./IC20/CNV/Bin_auto/ds_"$ratio"_"$i".log" 2>&1 &
		nohup python WGS_hg19_pe_without_control_auto.py "./IC35/CNV/Bin_auto/ds_"$ratio"_"$i "./IC35/IC35_"$ratio"_"$i".bam" > "./IC35/CNV/Bin_auto/ds_"$ratio"_"$i".log" 2>&1 &
		nohup python WGS_hg19_pe_without_control_auto.py "./IC37/CNV/Bin_auto/ds_"$ratio"_"$i "./IC37/IC37_"$ratio"_"$i".bam" > "./IC37/CNV/Bin_auto/ds_"$ratio"_"$i".log" 2>&1 &
		sleep 1s
	done
	# sleep 1200s
done

