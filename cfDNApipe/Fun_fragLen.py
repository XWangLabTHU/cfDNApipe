# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 10:51:10 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError, fraglendistribution
import os
from .Configure import Configure


__metaclass__ = type


class fraglenplot(StepBase):
    def __init__(self, 
         bedInput = None, # list
         outputdir = None, # str
         maxLimit = 250,
         upstream = None,
         initStep = False,
         **kwargs):
        super(fraglenplot, self).__init__(initStep)
        if upstream is None:
            self.setInput('bedInput', bedInput)
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('bamInput')[1])))
            else:
                self.setOutput('outputdir', outputdir)
            
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            
            if upstream.__class__.__name__ == 'bam2bed':
                self.setInput('bedInput', upstream.getOutput('bedOutput'))
            else:
                raise commonError('Parameter upstream must from bam2bed.')
            
            self.setOutput('outputdir', self.getStepFolderPath())
        
        self.setParam('maxLimit', maxLimit)
        self.setOutput('plotOutput', [os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x)) + '_fraglen.png' for x in self.getInput('bedInput')])
        self.setOutput('npyOutput', [os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x)) + '_fraglen.npy' for x in self.getInput('bedInput')])
        
        finishFlag = self.stepInit(upstream)
        
        if finishFlag:
            self.excute(finishFlag)
        else:
            multi_run_len = len(self.getInput('bedInput'))
            
            for i in range(multi_run_len):
                print("Now, ploting fragment length distribution for " + self.getInput('bedInput')[i])
                fraglendistribution(bedInput = self.getInput('bedInput')[i], plotOutput = self.getOutput('plotOutput')[i], 
                                    binOutput = self.getOutput('npyOutput')[i], maxLimit = self.getParam('maxLimit'))
            
            self.excute(finishFlag, runFlag = False)























