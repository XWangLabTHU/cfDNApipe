#!/bin/bash

arr=("IC15" "IC17" "IC20" "IC35" "IC37")
bin=("auto" "500" "1000" "2000" "5000" "7000" "10000")

for sample in ${arr[@]}
do
	echo "Now, processing "$sample
	for i in ${bin[@]}
	do
		mkdir "./"$sample"/CNV/Bin_"$i"/chr20"
		nohup python "WGS_hg19_pe_without_control_"$i".py" "./"$sample"/CNV/Bin_"$i"/chr20" $sample"_chr20.bam" > "./"$sample"/CNV/Bin_"$i"/chr20.log" 2>&1 &
		sleep 1s
	done
done

