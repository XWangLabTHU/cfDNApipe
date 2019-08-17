# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:54:29 2019

@author: zhang
"""

import StepBase
reload(StepBase)
import Configure
import Fun_inputProcess, Fun_fastqc, Fun_identifyAdapter, Fun_adapterremoval
import Fun_bowtie2
from importlib import reload
reload(Fun_inputProcess)
reload(Fun_fastqc)
reload(Fun_identifyAdapter)
reload(Fun_adapterremoval)
reload(Fun_bowtie2)


res = Fun_fastqc.fastqc(fastqInput=['/home/zhangwei/test/inputs/test1_1.fq', '/home/zhangwei/test/inputs/test2_1.fq',  
                                    '/home/zhangwei/test/inputs/test1_2.fq', '/home/zhangwei/test/inputs/test2_2.fq'],
                        fastqcOutputDir = '/home/zhangwei/test/outputs')


res = Fun_identifyAdapter.identifyAdapter(fqInput1 = ['/home/zhangwei/test/inputs/test1_1.fq', '/home/zhangwei/test/inputs/test2_1.fq'],
                                          fqInput2 = ['/home/zhangwei/test/inputs/test1_2.fq', '/home/zhangwei/test/inputs/test2_2.fq'],
                                          outputdir = '/home/zhangwei/test/outputs')


res = Fun_adapterremoval.adapterremoval(fqInput1 = ['/home/zhangwei/test/inputs/test1_1.fq', '/home/zhangwei/test/inputs/test2_1.fq'],
                                        fqInput2 = ['/home/zhangwei/test/inputs/test1_2.fq', '/home/zhangwei/test/inputs/test2_2.fq'],
                                        outputdir = '/home/zhangwei/test/outputs',
                                        adapter1 = ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'],
                                        adapter2 = ['AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
                                        threads = 24)

res = Fun_bowtie2.bowtie2(seqInput1 = ['/home/wzhang/test/inputs/test1_1.fq', '/home/wzhang/test/inputs/test2_1.fq'],
                          seqInput2 = ['/home/wzhang/test/inputs/test1_2.fq', '/home/wzhang/test/inputs/test2_2.fq'], 
                          ref = r'/home/wzhang/genome/atacflowdata/hg19',
                          outputdir = '/home/wzhang/test/outputs',
                          threads = 20)

