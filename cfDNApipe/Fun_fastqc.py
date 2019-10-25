# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 18:27:32 2019

@author: zhang
"""


from .StepBase import StepBase
from .cfDNA_utils import flatten
import os
from .Configure import Configure


__metaclass__ = type


class fastqc(StepBase):
    def __init__(self, 
                 fastqInput = None,
                 fastqcOutputDir  = None,
                 threads = 1,
                 other_params = None,
                 upstream = None,
                 formerrun = None,
                 **kwargs):
        if upstream is None:
            # In this situation, input file and output path should be checked
            super(fastqc, self).__init__()
            self.setInput('fastqInputs', fastqInput)
            self.checkInputFilePath()
            if fastqcOutputDir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('fastqInputs')[0])))
            else:
                self.setOutput('outputdir', fastqcOutputDir)
            
            self.setParam('threads', threads)
            
        else:
            # here to check upstream source!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            # super init
            if formerrun is None:
                super(fastqc, self).__init__(upstream.getStepID())
            else:
                super(fastqc, self).__init__(formerrun.getStepID())
            
            # check Configure for running pipeline
            Configure.configureCheck()
            
            upstream.checkFilePath()
            
            self.setInput('fastqInputs', list(flatten([upstream.getOutput('fq1'), upstream.getOutput('fq2')])))
            self.checkInputFilePath()
            
            self.setOutput('outputdir', self.getStepFolderPath())
            
            self.setParam('threads', Configure.getThreads())
            
        if other_params is None:
            self.setParam('other_params', '')
        else:
            self.setParam('other_params', other_params)
            
        # create cmd
        cmd = self.cmdCreate(["fastqc",
                              '--outdir', self.getOutput('outputdir'),
                              '--threads', self.getParam('threads'),
                              self.getParam('other_params'),
                              self.inputs['fastqInputs']])

        self.setParam('cmd', cmd)
            
        finishFlag = self.stepInit(upstream)
        
        self.excute(finishFlag)
































