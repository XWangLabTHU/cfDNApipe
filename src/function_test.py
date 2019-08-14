# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:38:07 2019

@author: zhang
"""

from StepBase import StepBase
from Configure import Configure
import Fun_inputProcess, Fun_fastqc, Fun_identifyAdapter
from importlib import reload
reload(Fun_inputProcess)
reload(Fun_fastqc)
reload(Fun_identifyAdapter)

Configure.setGenome("hg19")
Configure.setRefDir(r'D:\pipeCheck\bismark_genome')
Configure.setThreads(5)
Configure.setOutDir(r'D:\pipeCheck')
Configure.pipeFolderInit()

res1 = Fun_inputProcess.inputprocess(inputFolder = r"D:\pipeCheck\inputs")

#res2 = fastqc(upstream = res1)





import StepBase
import Configure
import Fun_inputProcess, Fun_fastqc, Fun_identifyAdapter, Fun_adapterremoval
from importlib import reload
reload(StepBase)
reload(Fun_inputProcess)
reload(Fun_fastqc)
reload(Fun_identifyAdapter)
reload(Fun_adapterremoval)

Configure.Configure.setGenome("hg19")
Configure.Configure.setRefDir(r'/home/zhangwei/test/bismark_genome')
Configure.Configure.setThreads(5)
Configure.Configure.setOutDir(r'/home/zhangwei/test')
Configure.Configure.pipeFolderInit()

#a = inputprocess(fqInput1 = [r'D:\pipeCheck\test1_1.fq', r'D:\pipeCheck\test2_1.fq'],
#                 fqInput2 = [r'D:\pipeCheck\test1_2.fq', r'D:\pipeCheck\test2_2.fq'])

res1 = Fun_inputProcess.inputprocess(inputFolder = r"/home/zhangwei/test/inputs")
res2 = Fun_fastqc.fastqc(upstream = res1)

res3 = Fun_identifyAdapter.identifyAdapter(upstream = res1, formerrun = res2)

res4 = Fun_adapterremoval.adapterremoval(upstream = res3)















