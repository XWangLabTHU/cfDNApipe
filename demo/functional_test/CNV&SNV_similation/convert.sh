#!/bin/bash

arr=("01" "05" "1" "2" "3" "4" "5" "6" "7" "8" "9")

for ratio in ${arr[@]}
do
	for i in $(seq 1 1 10)
	do
		cnvkit.py export bed "./ds_"$ratio"_"$i"/intermediate_result/step_CNV01_cnvbatch/IC17_"$ratio"_"$i".cnr" -o "./bed/IC17_"$ratio"_"$i".bed"
	done
done