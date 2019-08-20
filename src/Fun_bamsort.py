# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 11:45:21 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""

from StepBase import StepBase
from cfDNA_utils import commonError
import os, re
from Configure import Configure


__metaclass__ = type


class bamsort(StepBase):
    def __init__(self, 
             bamInput = None, # list
             outputdir = None, # str
             threads = 1,
             upstream = None,
             formerrun = None,
             **kwargs):
        if upstream is None:
            super(bamsort, self).__init__()
            self.setInput('bamInput', bamInput)
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('bamInput')[1])))
            else:
                self.setOutput('outputdir', outputdir)
            
            self.setParam('threads', threads)
            
        else:
            if formerrun is None:
                super(bamsort, self).__init__(upstream.getStepID())
            else:
                super(bamsort, self).__init__(formerrun.getStepID())
                
            # check Configure for running pipeline
            Configure.configureCheck()
            upstream.checkFilePath()
            
            if upstream.__class__.__name__ == 'bowtie2':
                self.setInput('bamInput', upstream.getOutput('bamOutput'))
            elif upstream.__class__.__name__ == 'bismark':
                self.setInput('bamInput', upstream.getOutput('bamOutput'))
            else:
                raise commonError('Parameter upstream must from inputprocess or adapterremoval.')
            
            self.setOutput('outputdir', self.getStepFolderPath())
            self.setParam('threads', Configure.getThreads())
            
        
        self.setOutput('bamOutput', [os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x)) + '-sorted.bam' for x in self.getInput('bamInput')])
        
        self.setOutput('baiOutput', [x + '.bai' for x in self.getOutput('bamOutput')])
        
        multi_run_len = len(self.getInput('bamInput'))
        
        all_cmd = []
        
        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(['samtools sort',
                                      '-@', self.getParam('threads'),
                                      '-o', self.getOutput('bamOutput')[i],
                                      self.getInput('bamInput')[i],
                                      ';',
                                      'samtools index',
                                      '-@', self.getParam('threads'),
                                      self.getOutput('bamOutput')[i]])
            all_cmd.append(tmp_cmd)
        
        self.setParam('cmd', all_cmd)
        
        finishFlag = self.stepInit(upstream)
        
        self.excute(finishFlag)

















































