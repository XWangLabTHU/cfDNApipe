# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 11:17:26 2019

@author: Jiaqi Huang

E-mail: huangjq16@mails.tsinghua.edu.cn
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError
import os, re
from .Configure import Configure


__metaclass__ = type


class bismark_deduplicate(StepBase):
    def __init__(self, 
             bamInput = None, # list
             outputdir = None, # str
             threads = 1,
             other_params = {},
             upstream = None,
             formerrun = None,
             **kwargs):
        if upstream is None:
            super(bismark_deduplicate, self).__init__()
            self.setInput('bamInput', bamInput)
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('bamInput')[1])))
            else:
                self.setOutput('outputdir', outputdir)
                
            self.setParam('threads', threads)
            
        else:
            # here to check upstream source!!!!!!!!!!!!!!!!!!!!!!!!!!
            # upstream can come from bismark
            
            if formerrun is None:
                super(bismark_deduplicate, self).__init__(upstream.getStepID())
            else:
                super(bismark_deduplicate, self).__init__(formerrun.getStepID())
                
            # check Configure for running pipeline
            Configure.configureCheck()
            
            upstream.checkFilePath()
            
            self.setParam('type', Configure.getType())
            if self.getParam('type') == 'paired':
                other_params.update{--paired = True}
            elif self.getParam('type') == 'single':
                other_params.update{--single = True}
            else:
                commonError("Wrong data type, must be 'single' or 'paired'!")
            
            if upstream.__class__.__name__ == 'bismark':
                self.setInput('bamInput', upstream.getOutput('bamOutput'))
            else:
                raise commonError('Parameter upstream must from bismark.')
            
            self.setOutput('outputdir', self.getStepFolderPath())
            self.setParam('threads', Configure.getThreads())
        
    
        if other_params is None:
            self.setParam('other_params', '')
        else:
            self.setParam('other_params',  other_params)
        
        self.setOutput('bamOutput', [self.getOutput('outputdir') + '/' + x.split('/')[-1][ : -3] + 'deduplicated.bam' for x in self.getInput('bamInput')])
        
        multi_run_len = len(self.getInput('bamInput'))
        
        all_cmd = []
        
        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(["deduplicate_bismark", 
                                       self.getParam('other_params'),
                                       '--bam', self.getInput('bamInput'),
                                       '--output_dir', self.getOutput('outputdir')])
            all_cmd.append(tmp_cmd)
        
        self.setParam('cmd', all_cmd)
        
        finishFlag = self.stepInit(upstream)
        
        self.excute(finishFlag)



