# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:38:07 2019

@author: zhang
"""

from cfDNApipe import *

Configure.setData('WGS')
Configure.setThreads(10)
Configure.setGenome("hg19")
Configure.setRefDir(r'/home/wzhang/genome/hg19')
Configure.setOutDir(r'/home/wzhang/test')
Configure.pipeFolderInit()
Configure.refCheck(build = False)


res1 = inputprocess(inputFolder = r"/home/wzhang/test/inputs")
res2 = fastqc(upstream = res1)
res3 = identifyAdapter(upstream = res1, formerrun = res2)
res4 = adapterremoval(upstream = res3)
res5 = bowtie2(upstream = res4)
res6 = bamsort(upstream = res5)
res7 = rmduplicate(upstream = res6)
res8 = bam2bed(upstream = res7)
res9 = fraglenplot(upstream = res8)


from cfDNApipe import *
Configure.setData('WGBS')
Configure.setThreads(20)
Configure.setGenome("hg19")
Configure.setRefDir(r'/home/wzhang/genome/hg19_bismark')
Configure.setOutDir(r'/data/wzhang/pipeline-test')
Configure.pipeFolderInit()
Configure.refCheck(build = True)


res1 = inputprocess(inputFolder = r"/data/wzhang/pipeline-test/raw-data")
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














