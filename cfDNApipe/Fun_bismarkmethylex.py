# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 18:35:24 2019

@author: Jiaqi Huang

E-mail: huangjq16@mails.tsinghua.edu.cn
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os, re
from .Configure import Configure


__metaclass__ = type


class bismark_methylation_extractor(StepBase):
    def __init__(self, 
             bamInput = None, # list
             outputdir = None, # str
             genomedir = None,
             threads = 1,
             other_params = {'--no_overlap': True, '--report': True, 
                             '--no_header': True, '--gzip': True, '--multicore': 4, 
                             '--cytosine_report': True, '--zero_based': True},
             upstream = None,
             formerrun = None,
             **kwargs):
        if upstream is None:
            super(bismark_methylation_extractor, self).__init__()
            self.setInput('bamInput', bamInput)
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('bamInput')[1])))
            else:
                self.setOutput('outputdir', outputdir)
                
            self.setParam('threads', threads)
            
        else:
            # here to check upstream source!!!!!!!!!!!!!!!!!!!!!!!!!!
            # upstream can come from bismark or bismark_deduplicate
            
            if formerrun is None:
                super(bismark_methylation_extractor, self).__init__(upstream.getStepID())
            else:
                super(bismark_methylation_extractor, self).__init__(formerrun.getStepID())
                
            # check Configure for running pipeline
            Configure.configureCheck()
            
            upstream.checkFilePath()
           
            self.setParam('type', Configure.getType())
            
            if self.getParam('type') == 'paired':
                other_params.update({'--paired': True})
            elif self.setParam('type') == 'single':
                other_params.update({'--single': True})
            else:
                commonError("Wrong data type, must be 'single' or 'paired'!")
            
            
            if upstream.__class__.__name__ == 'bismark' or 'bismark_deduplicate':
                self.setInput('bamInput', upstream.getOutput('bamOutput'))
            else:
                raise commonError('Parameter upstream must from bismark or bismark_deduplicate.')
            
            self.setOutput('outputdir', self.getStepFolderPath())
            self.setParam('threads', Configure.getThreads())
        
    
        if other_params is None:
            self.setParam('other_params', '')
        else:
            self.setParam('other_params',  other_params)
            
        self.setOutput('gzOutput', [self.getOutput('outputdir') + '/' + x.split('/')[-1][ : -3] + 'CpG_report.txt.gz' for x in self.getInput('bamInput')])
        
        self.setInput('genomedir', Configure.getRefDir())
        
        multi_run_len = len(self.getInput('bamInput'))
        
        all_cmd = []
        
        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(["bismark_methylation_extractor", 
                                       self.getParam('other_params'),
                                       '--output', self.getOutput('outputdir'),
                                       '--genome_folder', self.getInput('genomedir'),
                                       self.getInput('bamInput')[i]])
            all_cmd.append(tmp_cmd)
        
        self.setParam('cmd', all_cmd)
        
        finishFlag = self.stepInit(upstream)
        
        self.excute(finishFlag)



