# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:38:07 2019

@author: zhang
"""




import cfDNA_utils
import StepBase
import Configure
import Fun_inputProcess, Fun_fastqc, Fun_identifyAdapter, Fun_adapterremoval
import Fun_bowtie2, Fun_bismark, Fun_bamsort, Fun_rmDuplicate, Fun_bam2bed
import Fun_fragLen, Fun_OCF
from importlib import reload
reload(cfDNA_utils)
reload(StepBase)
reload(Fun_inputProcess)
reload(Fun_fastqc)
reload(Fun_identifyAdapter)
reload(Fun_adapterremoval)
reload(Fun_bowtie2)
reload(Fun_bamsort)
reload(Fun_rmDuplicate)
reload(Fun_bam2bed)
reload(Fun_fragLen)
reload(Fun_OCF)

Configure.Configure.setGenome("hg19")
Configure.Configure.setRefDir(r'/home/wzhang/genome/hg19')
Configure.Configure.setThreads(5)
Configure.Configure.setOutDir(r'/home/wzhang/test')
Configure.Configure.pipeFolderInit()


res1 = Fun_inputProcess.inputprocess(inputFolder = r"/home/wzhang/test/inputs")
res2 = Fun_fastqc.fastqc(upstream = res1)
res3 = Fun_identifyAdapter.identifyAdapter(upstream = res1, formerrun = res2)
res4 = Fun_adapterremoval.adapterremoval(upstream = res3)
res5 = Fun_bowtie2.bowtie2(upstream = res4)
res6 = Fun_bamsort.bamsort(upstream = res5)
res7 = Fun_rmDuplicate.rmduplicate(upstream = res6)
res8 = Fun_bam2bed.bam2bed(upstream = res7)
res9 = Fun_fragLen.fraglenplot(upstream = res8)


Configure.Configure.setGenome("hg19")
Configure.Configure.setRefDir(r'/home/wzhang/genome/hg19_bismark')
Configure.Configure.setThreads(20)
Configure.Configure.setOutDir(r'/data/wzhang/pipeline-test')
Configure.Configure.pipeFolderInit()


res1 = Fun_inputProcess.inputprocess(inputFolder = r"/data/wzhang/pipeline-test/raw-data")
res2 = Fun_fastqc.fastqc(upstream = res1)
res3 = Fun_identifyAdapter.identifyAdapter(upstream = res1, formerrun = res2)
res4 = Fun_adapterremoval.adapterremoval(upstream = res3)
res5 = Fun_bismark.bismark(upstream = res4)
res6 = Fun_bamsort.bamsort(upstream = res5)
res7 = Fun_rmDuplicate.rmduplicate(upstream = res6)
res8 = Fun_bam2bed.bam2bed(upstream = res7)
res9 = Fun_fragLen.fraglenplot(upstream = res8)
res10 = Fun_OCF.computeOCF(upstream = res8, formerrun = res9, refRegInput = '/data/wzhang/OCF-test/OCF-regions.bed')
















