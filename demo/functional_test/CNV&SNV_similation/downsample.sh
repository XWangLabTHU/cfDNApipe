#!/bin/bash

arr=("01" "05" "1" "2" "3" "4" "5" "6" "7" "8" "9")

for ratio in ${arr[@]}
do
	for i in $(seq 1 1 10)
	do
		nohup samtools view -bhs $i"."$ratio -o "IC17_"$ratio"_"$i".bam" IC17_chr20.bam > "IC17_"$ratio"_"$i".log" 2>&1 &
	done
done