# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 19:44:37 2019

@author: zhang
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure


__metaclass__ = type


class identifyAdapter(StepBase):
    def __init__(self, 
                 fqInput1 = None,
                 fqInput2 = None,
                 outputdir = None,
                 threads = 1,
                 other_params = None,
                 upstream = None,
                 initStep = False,
                 **kwargs):
        super(identifyAdapter, self).__init__(initStep)
        if upstream is None:
            self.setInput('fq1', fqInput1)
            self.setInput('fq2', fqInput2)
            self.checkInputFilePath()
            
            self.setOutput('fq1', fqInput1)
            self.setOutput('fq2', fqInput2)
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('fq1')[0])))
            else:
                self.setOutput('outputdir', outputdir)
            
            self.setParam('threads', threads)
            
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            
            self.setInput('fq1', upstream.getOutput('fq1'))
            self.setInput('fq2', upstream.getOutput('fq2'))
            self.checkInputFilePath()

            self.setOutput('fq1', upstream.getOutput('fq1'))
            self.setOutput('fq2', upstream.getOutput('fq2'))
            
            self.setOutput('outputdir', self.getStepFolderPath())
            
            self.setParam('threads', Configure.getThreads())

        if other_params is None:
            self.setParam('other_params', '')
        else:
            self.setParam('other_params', other_params)
        
        if len(self.getOutput('fq1')) == len(self.getOutput('fq2')):
            multi_run_len = len(self.getOutput('fq1'))
        else:
            raise commonError('Paired end Input files are not consistent.')
        
        all_cmd = []
        
        for i in range(multi_run_len):
            tmp_fq1, tmp_fq2  = self.getInput('fq1')[i], self.getInput('fq2')[i]
            tmp_basename = self.getMaxFileNamePrefix(tmp_fq1, tmp_fq2)
            self.setOutput(tmp_basename + '-adapterFile', os.path.join(self.getOutput('outputdir'), tmp_basename + '-adapters.log'))
        
            tmp_cmd = self.cmdCreate(["AdapterRemoval --identify-adapters", 
                                       '--threads', self.getParam('threads'),
                                       '--file1', tmp_fq1,
                                       '--file2', tmp_fq2,
                                       self.getParam('other_params'),
                                       '>', self.getOutput(tmp_basename + '-adapterFile')])
            all_cmd.append(tmp_cmd)
        
        self.setParam('cmd', all_cmd)
        
        finishFlag = self.stepInit(upstream)
        
        self.excute(finishFlag)







