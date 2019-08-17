# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 10:03:52 2019

@author: zhang
"""
import re
adapter = []
with open('test.log', 'r') as f:
    lines = f.readlines()
    for line in lines:
        if re.search(r'Consensus:', line):
            tmp_line = line.split()
            adapter.append(tmp_line[1])

print(adapter)
            


def getAdapetrFromFile(file):
    adapter = []
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if re.search(r'Consensus:', line):
                tmp_line = line.split()
                adapter.append(tmp_line[1])
                
cmd = "bowtie2 -x /home/wzhang/genome/atacflowdata/hg19 -1 /home/wzhang/test/inputs/test2_1.fq -2 /home/wzhang/test/inputs/test2_2.fq -q -N 1 -X 2000 --no-mixed --no-discordant --dovetail --time --un-conc-gz /home/wzhang/test/outputs/test2.gz -p 20 | samtools view -b -S -@ 20 - | samtools sort -@ 20 -o /home/wzhang/test/outputs/test2.bam - > ./111.log"

import subprocess
import sys
with open('test.log', 'wb') as f:  # replace 'w' with 'wb' for Python 3
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)
    for c in iter(lambda: process.stdout.read(1), b''):  # replace '' with b'' for Python 3
        sys.stdout.write(c)
        f.write(c)




import subprocess
import sys
sys.stdout = open('111.log', 'w')
print('test')
grepout = subprocess.check_output(cmd, shell = True)
print('test')



nohup bowtie2 -x /home/wzhang/genome/atacflowdata/hg19 -1 /home/wzhang/test/inputs/test2_1.fq -2 /home/wzhang/test/inputs/test2_2.fq \
-q -N 1 -X 2000 --no-mixed --no-discordant --dovetail --time --un-conc-gz /home/wzhang/test/outputs/test2.gz -p 20 \
| samtools view -b -S -@ 20 - | samtools sort -@ 20 -o /home/wzhang/test/outputs/test2.bam - > ./111.log






















