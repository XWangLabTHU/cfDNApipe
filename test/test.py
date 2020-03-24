# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 19:32:30 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""

from matplotlib import pyplot as plt
plt.plot(range(10))

# 知道的染色质在横轴的位置，假如1-2，那么中间就在1.5, 具体位置你再调整一下
plt.plot([1, 2], [1, 1], 'k-', lw=2)
plt.text(x=(1+2) / 2, y=1, s="chr1p", horizontalalignment='center')
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)  # labels along the bottom edge are off
plt.show()
plt.savefig('plot')
plt.clf()




finishFlag = self.stepInit(caseupstream)  # need to be checked

if not finishFlag:
    function1(1)
    self.run(cmd1)
    self.run(cmd2)
    function2(2)
    self.run(cmd3)

self.stepInfoRec(cmds=[cmd1, cmd2, cmd3], finishFlag=finishFlag)
