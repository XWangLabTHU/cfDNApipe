#!/bin/bash

arr=("01" "05" "1" "2" "3" "4" "5" "6" "7" "8" "9")

for ratio in ${arr[@]}
do
	for i in $(seq 1 1 10)
	do
		mkdir "ds_"$ratio"_"$i
		nohup python WGS_hg19_pe_without_control.py "ds_"$ratio"_"$i "IC17_"$ratio"_"$i".bam" > "IC17_"$ratio"_"$i".log" 2>&1 &
	done
	sleep 40m
done

