# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:38:07 2019

@author: zhang
"""


from cfDNApipe import *
Configure.setData('WGBS')
Configure.setThreads(20)
Configure.setGenome("hg19")
Configure.setRefDir(r'/home/wzhang/genome/hg19_bismark')
Configure.setOutDir(r'/data/wzhang/pipeline-for-paired-WGBS')
Configure.pipeFolderInit()
Configure.refCheck(build = True)


res1 = inputprocess(inputFolder = r"/data/wzhang/pipeline-for-paired-WGBS/raw")
res2 = fastqc(upstream = res1)
res3 = identifyAdapter(upstream = res1, formerrun = res2)
res4 = adapterremoval(upstream = res3)
res5 = bismark(upstream = res4)
res6 = bamsort(upstream = res5)
res7 = rmduplicate(upstream = res6)
res8 = bam2bed(upstream = res7)
res9 = fraglenplot(upstream = res8)
res10 = computemethyl(upstream = res7, formerrun = res9)
res11 = addRG(upstream = res7, formerrun = res10)




# single end WGBS

from cfDNApipe import *
Configure.setData('WGBS')
Configure.setType('single')
Configure.setThreads(20)
Configure.setGenome("hg19")
Configure.setRefDir(r'/home/wzhang/genome/hg19_bismark')
Configure.setOutDir(r'/data/wzhang/pipeline-for-single-WGBS')
Configure.pipeFolderInit()
Configure.refCheck(build = True)

res1 = inputprocess(inputFolder = r"/data/wzhang/pipeline-for-single-WGBS/raw")
res2 = fastqc(upstream = res1)
res3 = adapterremoval(upstream = res1, formerrun = res2, other_params = {'--qualitybase': 64, '--gzip': True})
res4 = bismark(upstream = res3, other_params = {'-q': True, '--phred64-quals': True,  '-N': 1, '--bowtie2': True, '--un': True})
res5 = bamsort(upstream = res4)
res6 = rmduplicate(upstream = res5)
res7 = bam2bed(upstream = res6)
res8 = computemethyl(upstream = res6, formerrun = res7)
res9 = addRG(upstream = res6, formerrun = res8)













