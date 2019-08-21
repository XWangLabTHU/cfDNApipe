# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 10:03:52 2019

@author: zhang
"""


import numpy as np 
from matplotlib import pyplot as plt 

x = np.arange(1,11) 
y =  np.sin(x) +  5 

fig = plt.figure()
plt.plot(x, y)
plt.savefig(r'C:\\Users\\Zhang\\Desktop\\MyFig.jpg')
plt.close(fig)