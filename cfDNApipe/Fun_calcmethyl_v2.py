# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 15:13:42 2020

@author: Jiaqi Huang
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError, calcMethylV2
from .Configure import Configure
import os


__metaclass__ = type


class calculate_methyl(StepBase):
    def __init__(self, 
         tbxInput = None, # list
         bedInput = None,
         outputdir = None, # str
         threads = 1,
         upstream = None,
         initStep = False,
         **kwargs):
        
        super(calculate_methyl, self).__init__(initStep)
        
        if upstream is None:
            self.setInput('tbxInput', tbxInput)
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('covgzInput')[0])))
            else:
                self.setOutput('outputdir', outputdir)
            
            self.setParam('threads', threads)
            
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            
            if upstream.__class__.__name__ == 'compress_methyl':
                self.setInput('tbxInput', upstream.getOutput('tbxOutput'))
            else:
                raise commonError('Parameter upstream must from compress_methyl.')
                
            self.setOutput('outputdir', self.getStepFolderPath())
            self.setParam('threads', Configure.getThreads())
        
        if bedInput is not None:
            self.setInput('bedInput', bedInput)
        else:
            raise commonError('No bed file input!')
        
        self.setOutput('txtOutput', [os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x)) + '-result.txt' for x in self.getInput('tbxInput')])
        
        finishFlag = self.stepInit(upstream)
        
        if finishFlag:
            self.excute(finishFlag)
        else:
            multi_run_len = len(self.getInput('tbxInput'))
            for i in range(multi_run_len):
                calcMethylV2(tbxInput = self.getInput('tbxInput')[i], bedInput = self.getInput('bedInput'), txtOutput = self.getOutput('txtOutput')[i])
               
            self.excute(finishFlag, runFlag = False)        

