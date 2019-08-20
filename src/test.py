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



import subprocess, shlex

cmd = "bowtie2 -x /home/wzhang/genome/atacflowdata/hg19 -1 /home/wzhang/test/inputs/test2_1.fq -2 /home/wzhang/test/inputs/test2_2.fq -q -N 1 -X 2000 --no-mixed --no-discordant --dovetail --time --un-conc-gz /home/wzhang/test/outputs/test2.gz -p 20 | samtools view -bS -@ 20 - | samtools sort -@ 20 -o /home/wzhang/test/outputs/test2.bam -"

result = subprocess.run(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)


subprocess.run("ls -l", shell = True)

subprocess.run(cmd, shell = True)





import shlex
from subprocess import Popen, PIPE

def get_exitcode_stdout_stderr(cmd):
    """
    Execute the external command and get its exitcode, stdout and stderr.
    """
    args = shlex.split(cmd)

    proc = Popen(args, stdout=PIPE, stderr=PIPE, shell = True)
    out, err = proc.communicate()
    exitcode = proc.returncode
    #
    return exitcode, out, err

cmd = "ls -al"
exitcode, out, err = get_exitcode_stdout_stderr(cmd)


cmd = "bowtie2 -x /home/wzhang/genome/atacflowdata/hg19 -1 /home/wzhang/test/inputs/test2_1.fq -2 /home/wzhang/test/inputs/test2_2.fq -q -N 1 -X 2000 --no-mixed --no-discordant --dovetail --time --un-conc-gz /home/wzhang/test/outputs/test2.gz -p 20"



import subprocess

def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

# Example
for path in execute(["locate", "a"]):
    print(path, end="")



cmd1 = "bowtie2 -x /home/wzhang/genome/atacflowdata/hg19 -1 /home/wzhang/test/inputs/test2_1.fq -2 /home/wzhang/test/inputs/test2_2.fq -q -N 1 -X 2000 --no-mixed --no-discordant --dovetail --time --un-conc-gz /home/wzhang/test/outputs/test2.gz -p 20"
cmd2 = "samtools view -bS -@ 20 -o /home/wzhang/test/outputs/test2.bam"

import subprocess
p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell = True)
p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE, shell = True)
p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
output,err = p2.communicate()







p1 = Popen(["dmesg"], stdout=PIPE)
p2 = Popen(["grep", "nvidia"], stdin=p1.stdout, stdout=PIPE)
p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
output = p2.communicate()[0]


ps = subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell = True)
output = subprocess.check_output(cmd2, stdin=ps.stdout, shell = True)



111111111111111111111111111
cmd = "bowtie2 -x /home/wzhang/genome/atacflowdata/hg19 -1 /home/wzhang/test/inputs/test2_1.fq -2 /home/wzhang/test/inputs/test2_2.fq -q -N 1 -X 2000 --no-mixed --no-discordant --dovetail --time --un-conc-gz /home/wzhang/test/outputs/test2.gz -p 20 | samtools view -b -S -@ 20 - | samtools sort -@ 20 -o /home/wzhang/test/outputs/test2.bam - "

proc = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = proc.communicate()
errcode = proc.returncode

mess = stderr.decode(sys.getfilesystemencoding())
print(mess)

1111111111111111111111111111


cmd1 = "bowtie2 -x /home/wzhang/genome/atacflowdata/hg19 -1 /home/wzhang/test/inputs/test2_1.fq -2 /home/wzhang/test/inputs/test2_2.fq -q -N 1 -X 2000 --no-mixed --no-discordant --dovetail --time --un-conc-gz /home/wzhang/test/outputs/test2.gz -p 20"
cmd2 = "samtools view -bS -@ 20 -o /home/wzhang/test/outputs/test2.bam"


import subprocess
p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE, shell = True)
p1.stdout.close()
output,err = p2.communicate()

p1.stdout.readlines()
p1.stderr.readlines()


cmd = "bowtie2 -x /home/wzhang/genome/atacflowdata/hg19 -1 /home/wzhang/test/inputs/test2_1.fq -2 /home/wzhang/test/inputs/test2_2.fq -q -N 1 -X 2000 --no-mixed --no-discordant --dovetail --time --un-conc-gz /home/wzhang/test/outputs/test2.gz -p 20 | samtools view -b -S -@ 20 - | samtools sort -@ 20 -o /home/wzhang/test/outputs/test2.bam - "



def execute(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)

    # Poll process for new output until finished
    while True:
        nextline = process.stdout.readline()
        if (nextline == '') and (process.poll() is not None):
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()

    output, error = process.communicate()
    exitCode = process.returncode

    if exitCode == 0:
        return output, error, exitCode
    else:
        raise Exception(command, exitCode, error)



def execute(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Poll process for new output until finished
    while True:
        nextline = process.stderr.readline()
        sys.stdout.write(nextline)
        sys.stdout.flush()
        print(type(nextline))

    if exitCode == 0:
        return output, error
    else:
        output = output.decode(sys.getfilesystemencoding())
        raise Exception(command, exitCode, error)











