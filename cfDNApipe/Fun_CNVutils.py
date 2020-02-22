# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 11:45:21 2019

@author: Jiaqi Huang

"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure


__metaclass__ = type


class CNV_utils(StepBase):
    def __init__(self, 
             bigwigInput = None,
             fastaInput = None,
             bamInput = None, # list
             outputdir = None, # str
             threads = 1,
             upstream = None,
             initStep = False,
             **kwargs):
        super(CNV_utils, self).__init__(initStep)
        if upstream is None:
            self.setInput('bamInput', bamInput)
            self.setInput('bigwigInput', bigwigInput)
            self.setInput('fastaInput', fastaInput)
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('bamInput')[1])))
            else:
                self.setOutput('outputdir', outputdir)
            
            self.setParam('threads', threads)
            
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            
            if upstream.__class__.__name__ == 'bamsort':
                self.setInput('bamInput', upstream.getOutput('bamOutput'))
            else:
                raise commonError('Parameter upstream must from bamsort.')
            
            self.setInput('bigwigInput', bigwigInput)
            self.setInput('fastaInput', fastaInput)
            self.checkInputFilePath()
        
            self.setOutput('outputdir', self.getStepFolderPath())
            self.setParam('threads', Configure.getThreads())
              
        self.setOutput('mapOutput', os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(self.getInput('bigwigInput')) + '.map.wig')
        self.setOutput('gcOutput', os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(self.getInput('fastaInput')) + '.gc.wig')
        self.setOutput('readOutput', [os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x) + '.read.wig') for x in self.getInput('bamInput')])
        
        multi_run_len = len(self.getInput('bamInput'))
        
        all_cmd = []
        map_tmp_cmd = self.cmdCreate(['./mapCounter',
                                      '-w', 100000,
                                      self.getInput('bigwigInput'),
                                      '>', self.getOutput('mapOutput')])
        all_cmd.append(map_tmp_cmd)
        gc_tmp_cmd = self.cmdCreate(['./gcCounter',
                                      '-w', 100000,
                                      self.getInput('fastaInput'),
                                      '>', self.getOutput('gcOutput')])
        all_cmd.append(gc_tmp_cmd)
        for i in range(multi_run_len):
            read_tmp_cmd = self.cmdCreate(['./readCounter',
                                           '-w', 100000,
                                           self.getInput('bamInput')[i],
                                           '>', self.getOutput('readOutput')[i]])
            all_cmd.append(read_tmp_cmd)
        
        self.setParam('cmd', all_cmd)
        
        finishFlag = self.stepInit(upstream)
        
        self.excute(finishFlag)

















































