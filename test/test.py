# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 19:32:30 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""

import os 

def isAlphaOrDigit(x):
    if x.isalpha() or x.isdigit():
        return True
    else:
        return False
    
def rmEndString(x, y):
    for item in y:
        if x.endswith(item):
            x = x.replace(item, '')
    
    return x

def getFileNamePrefix(fileName):
    return os.path.split(fileName)[-1]

def getMaxFileNamePrefix(file1, file2):
    file1 = getFileNamePrefix(file1)
    file2 = getFileNamePrefix(file2)
    len1 = len(file1)
    len2 = len(file2)
    for i in range(min(len1, len2)):
        if file1[i] != file2[i]:
            break
    if i == 0:
        return ''
    elif i == (min(len1, len2) - 1):
        tmp_str = file1[: i + 1]
    else:
        tmp_str = file1[: i]
    
    for k in reversed(range(len(tmp_str))):
        if isAlphaOrDigit(tmp_str[k]):
            final_name = tmp_str[: k + 1]
            final_name = rmEndString(final_name, ['.pair', '.read'])
            return final_name
        else:
            k = k - 1


def getMaxFileNamePrefixV2(file):
    final_name = os.path.splitext(os.path.basename(file))[0]
    final_name = rmEndString(final_name, ['_pe', '-sorted', '-rmdup', '.fq.gz', '.fq'])
    return final_name

            
            