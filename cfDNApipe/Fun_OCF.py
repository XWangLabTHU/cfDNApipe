# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:51:54 2019

@author: zhang
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
import pybedtools
from .Configure import Configure


__metaclass__ = type


class computeOCF(StepBase):
    def __init__(self, 
         bedInput = None, # list
         refRegInput = None,
         outputdir = None, # str
         upstream = None,
         initStep = False,
         **kwargs):
        super(computeOCF, self).__init__(initStep)
        if upstream is None:
            self.setInput('bedInput', bedInput)
            self.setInput('refRegInput', refRegInput)
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('bedInput')[0])))
            else:
                self.setOutput('outputdir', outputdir)
            
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            
            if upstream.__class__.__name__ == 'bam2bed':
                self.setInput('bedInput', upstream.getOutput("bedOutput"))
            else:
                raise commonError('Parameter upstream must from bam2bed.')
            
            self.setInput('refRegInput', refRegInput)
            self.setOutput('outputdir', self.getStepFolderPath())
        
        self.setOutput('txtOutput', [os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x)) + '.txt' for x in self.getInput('bedInput')])
        
        finishFlag = self.stepInit(upstream)
        
        if finishFlag:
            self.excute(finishFlag)
        else:
            multi_run_len = len(self.getInput('bedInput'))
            
            for i in range(multi_run_len):
                print("Now, processing file: " + self.getInput('bedInput')[i])
                a = pybedtools.BedTool(self.getInput('bedInput')[i])
                b = pybedtools.BedTool(self.getInput('refRegInput'))
                a.intersect(b, wo = True, sorted = True, output=self.getOutput('txtOutput')[i])
            
            self.excute(finishFlag, runFlag = False)
            
            
            
            
            