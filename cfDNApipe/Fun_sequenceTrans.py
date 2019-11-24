# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 18:27:32 2019

@author: Jiqi Huang, Chang Li
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError
from .Configure import Configure
from collections import defaultdict
import os


__metaclass__ = type


class sequencetransfer(StepBase):
    def __init__(self, 
         bamInput = None, # list
         bedInput = None,
         outputdir = None, # str
         threads = 1,
         upstream = None,
         formerrun = None,
         **kwargs):
        
        bedflag = False
        
        if upstream is None:
            super(sequencetransfer, self).__init__()
            self.setInput('bamInput', bamInput)
            
            if bedInput is not None:
                bedflag = True
                self.setInput('bedInput', bedInput)
                
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('bamInput')[0])))
            else:
                self.setOutput('outputdir', outputdir)
            
            self.setParam('threads', threads)
        else:
            if formerrun is None:
                super(sequencetransfer, self).__init__(upstream.getStepID())
            else:
                super(sequencetransfer, self).__init__(formerrun.getStepID())
                
            Configure.configureCheck()
            upstream.checkFilePath()
            
            if upstream.__class__.__name__ == 'rmduplicate':
                self.setInput('bamInput', upstream.getOutput('bamOutput'))
            else:
                raise commonError('Parameter upstream must from inputprocess or adapterremoval.')
            
            if bedInput is not None:
                bedflag = True
                self.setInput('bedInput', bedInput)
                
            self.setOutput('outputdir', self.getStepFolderPath())
            self.setParam('threads', Configure.getThreads())
        
        self.setOutput('txtOutput', [os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x)) + '-seq.txt' for x in self.getInput('bamInput')])
        
        finishFlag = self.stepInit(upstream)
        
        if finishFlag:
            self.excute(finishFlag)
        else:
            multi_run_len = len(self.getInput('bamInput'))
            for i in range(multi_run_len):
                if bedflag:
                    self.seqTrans(bamInput = self.getInput('bamInput')[i], txtOutput = self.getOutput('txtOutput')[i], bedInput = self.getInput('bedInput'))
                else:
                    self.seqTrans(bamInput = self.getInput('bamInput')[i], txtOutput = self.getOutput('txtOutput')[i])
               
            self.excute(finishFlag, runFlag = False)        

