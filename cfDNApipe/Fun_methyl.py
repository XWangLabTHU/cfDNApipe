# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 18:27:32 2019

@author: Jiqi Huang, Chang Li
"""


import os, pysam
from .StepBase import StepBase
from .Configure import Configure
from .cfDNA_utils import commonError


__metaclass__ = type


class computemethyl(StepBase):
    def __init__(self, 
         bamInput = None, # list
         bedInput = None,
         outputdir = None, # str
         threads = 1,
         upstream = None,
         formerrun = None,
         **kwargs):
        if upstream is None:
            super(computemethyl, self).__init__()
            self.setInput('bamInput', bamInput)
            self.setInput('bedInput', bedInput)
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('bamInput')[0])))
            else:
                self.setOutput('outputdir', outputdir)
            
            self.setParam('threads', threads)
            
        else:
            if formerrun is None:
                super(computemethyl, self).__init__(upstream.getStepID())
            else:
                super(computemethyl, self).__init__(formerrun.getStepID())
                
            Configure.configureCheck()
            upstream.checkFilePath()
            
            if upstream.__class__.__name__ == 'rmduplicate':
                self.setInput('bamInput', upstream.getOutput('bamOutput'))
            else:
                raise commonError('Parameter upstream must from rmduplicate.')
            
            self.setInput('bedInput', bedInput)
            self.setOutput('outputdir', self.getStepFolderPath())
            self.setParam('threads', Configure.getThreads())
        
        self.setOutput('txtOutput', [os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x)) + '-methyl.txt' for x in self.getInput('bamInput')])
        
        count_data = [[0, 0, 0, 0, 0, 0, 0, 0]]  #8 columns represent: CpG_counts, CHG_counts, CHH_counts, unknown_C_counts, methylated_CpG_counts, methylated_CHG_counts, methylated_CHH_counts, methylated_unknown_C_counts
        lines = 0
        bed_input = open(self.getInput('bedInput'), 'r')
        bed_lines = bed_input.readlines()
        multi_run_len = len(self.getInput('bamInput'))
        for i in range(multi_run_len):
            bam_input = pysam.AlignmentFile(self.getInput('bamInput')[i], 'rb')
            txt_output = open(self.getOutput('txtOutput')[i], 'w')
            for j in range(len(bed_lines)):
                bed_line = bed_lines[j]
                bed_intv = bed_line.strip().split('\t')
                bam_intv = bam_input.fetch(reference = bed_intv[0], start = int(bed_intv[1]), end = int(bed_intv[2]))
                for r in bam_intv:
                    for char in r.get_tag('XM'):
                        if char == 'z':    #CpG unmethylated
                            count_data[lines][0] += 1
                        elif char == 'Z':    #CpG methylated
                            count_data[lines][0] += 1
                            count_data[lines][4] += 1
                        elif char == 'x':    #CHG unmethylated
                            count_data[lines][1] += 1
                        elif char == 'X':    #CHG methylated
                            count_data[lines][1] += 1
                            count_data[lines][5] += 1
                        elif char == 'h':    #CHH unmethylated
                            count_data[lines][2] += 1
                        elif char == 'H':    #CHH methylated
                            count_data[lines][2] += 1
                            count_data[lines][6] += 1
                        elif char == 'u':    #unknown C unmethylated
                            count_data[lines][3] += 1
                        elif char == 'U':    #unknown C methylated
                            count_data[lines][3] += 1
                            count_data[lines][7] += 1
                txt_output.write('Methylation ratio in ' + bed_intv[0] + ', ' + bed_intv[1] + ' - ' + bed_intv[2] + ':\n')
                txt_output.write('C methylated in CpG context: ' + str(count_data[lines][4] / (count_data[lines][0] + 1e-20)) + '\n')
                txt_output.write('C methylated in CHG context: ' + str(count_data[lines][5] / (count_data[lines][1] + 1e-20)) + '\n')
                txt_output.write('C methylated in CHH context: ' + str(count_data[lines][6] / (count_data[lines][2] + 1e-20)) + '\n')
                txt_output.write('C methylated in unknown context: ' + str(count_data[lines][7] / (count_data[lines][3] + 1e-20)) + '\n')
                txt_output.write('C total methylation ratio: ' + str((count_data[lines][4] + count_data[lines][5] + count_data[lines][6] + count_data[lines][7]) / (count_data[lines][0] + count_data[lines][1] + count_data[lines][2] + count_data[lines][3] + 1e-20)) + '\n')
                txt_output.write('--------------------------\n\n')
                lines += 1
                count_data.append([0, 0, 0, 0, 0, 0, 0, 0])
            bam_input.close()
            txt_output.close()
            
        bed_input.close()
        
        finishFlag = self.stepInit(upstream)
        
        self.excute(finishFlag, runFlag = False)
