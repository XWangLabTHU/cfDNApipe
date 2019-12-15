# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 19:32:30 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""


Configure2.setData('WGBS')
Configure2.setType('single')
Configure2.setThreads(20)
Configure2.setGenome('hg19')
Configure2.setRefDir(r'/home/wzhang/genome/hg19_bismark')
Configure2.setCase('Normal')
Configure2.setCtrl('Tumor')
Configure2.setOutDir(r'/data/wzhang/pipeline-for-comp')
Configure.pipeFolderInit()
Configure.refCheck(build = False)



