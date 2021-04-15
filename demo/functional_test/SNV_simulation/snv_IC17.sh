#!/bin/bash

# arr=("01" "05" "1" "2" "3" "4" "5" "6" "7" "8" "9")
arr=("01" "05")

for ratio in ${arr[@]}
do
	echo $ratio
	for i in $(seq 1 1 10)
	do
		mkdir "./IC17/SNV/ds_"$ratio"_"$i
		nohup python WGS_hg19_pe_without_control_snv.py "./IC17/SNV/ds_"$ratio"_"$i "./IC17/IC17_"$ratio"_"$i".bam" > "./IC17/SNV/ds_"$ratio"_"$i".log" 2>&1 &
		sleep 1s
	done
	# sleep 30m
done

