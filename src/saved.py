# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 09:55:10 2019

@author: zhang-wei
"""

__metaclass__ = type

class Entity():
    def __init__(self, size, x, y):
        self.x, self.y = x, y
        self.size = size

    def run(self):
        self.res1 = self.x + self.y


class a1(Entity):
    def __init__(self, a, b, c, d, e):
        super(a1, self).__init__(c, d, e)
        self.a = a
        self.b = b
        
    
    def walk(self):
        self.res2 = self.a * self.b



def plus(x, y, size):
    tmp = Entity(size, x, y)
    tmp.run()
    return(tmp)



from subprocess import STDOUT, check_output, CalledProcessError

try:
    grepOut = check_output("bowtie2 -h", shell = True)
    grepOut = grepOut.decode("utf-8")
except CalledProcessError as grepexc:
    print("error code", grepexc.returncode, grepexc.output)
    print(grepexc.output)
    
    
import sys
import subprocess

try:
    grepout = subprocess.check_output("fastqc -h", 
                                      shell = True)
    print(11111111111)
    print(grepout.decode(sys.getfilesystemencoding()))
    print(11111111111)
except subprocess.CalledProcessError as e:
    print('exit code: {}'.format(e.returncode))
    
    
    print('stdout: {}'.format(e.output.decode(sys.getfilesystemencoding())))
    print('stderr: {}'.format(e.stderr.decode(sys.getfilesystemencoding())))
    
    
from subprocess import Popen, PIPE

process = Popen('ls -al', stdout=PIPE, stderr=PIPE)
output, err = process.communicate()
    
    
    
from __future__ import print_function # Only Python 2.x
import subprocess

def execute(cmd):
    popen = subprocess.Popen(cmd, stdout = subprocess.PIPE, universal_newlines = True, shell = True)
    for stdout_line in iter(popen.stdout.readline, ""):
        print(1)
        print(stdout_line)
        print(1)
        yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

try:
    for path in execute("fastqc --outdir /home/zhangwei/test --threads 5  /home/zhangwei/test/test_1.fq /home/zhangwei/test/test_2.fq"):
        print(2)
        print(path, end = "")
        print(2)
except subprocess.CalledProcessError as e:
    print('exit code: {}'.format(e.returncode))
    print(e.cmd)
    print(e.output)





import subprocess

cmd1 = "fastqc --outdir /home/zhangwei/test --threads 5  /home/zhangwei/test/test_1.fq /home/zhangwei/test/test_2.fq"
cmd2 = "bowtie2 -h"
subprocess.run(cmd1, shell=True,stderr=subprocess.PIPE, stdout=None)



















    
    
    