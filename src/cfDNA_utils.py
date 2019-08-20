# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 11:37:46 2019

@author: zhang
"""


from collections import Iterable


__metaclass__ = type


class commonError(Exception):
    def __init__(self, message):
        self.message = message
        

# https://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists
def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el
            
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





