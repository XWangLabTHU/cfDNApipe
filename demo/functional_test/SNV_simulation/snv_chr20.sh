#!/bin/bash

arr=("IC15" "IC20" "IC35" "IC37")

for sample in ${arr[@]}
do
	nohup python WGS_hg19_pe_without_control_snv.py "./"$sample"/SNV/chr20" $sample"_chr20.bam" > "./"$sample"/SNV/chr20.log" 2>&1 &
	sleep 1s
done

