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