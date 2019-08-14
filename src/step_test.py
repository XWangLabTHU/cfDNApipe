# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:54:29 2019

@author: zhang
"""

import StepBase
reload(StepBase)
import Configure
import Fun_inputProcess, Fun_fastqc, Fun_identifyAdapter, Fun_adapterremoval
from importlib import reload
reload(Fun_inputProcess)
reload(Fun_fastqc)
reload(Fun_identifyAdapter)
reload(Fun_adapterremoval)


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



AdapterRemoval --file1 /home/zhangwei/test/inputs/test1_1.fq --file2 /home/zhangwei/test/inputs/test1_2.fq \
--adapter1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--basename /home/zhangwei/test/outputs/test1 --threads 24