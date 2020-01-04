# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 18:35:24 2019

@author: Jiaqi Huang

E-mail: huangjq16@mails.tsinghua.edu.cn
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os, math
from .Configure import Configure


__metaclass__ = type


class bismark_methylation_extractor(StepBase):
    def __init__(self, 
             bamInput = None, # list
             outputdir = None, # str
             threads = 1,
             other_params = {'--no_overlap': True, '--report': True,
                             '--no_header': True, '--gzip': True,
                             '--bedGraph': True, '--zero_based': True},
             paired = True,
             upstream = None,
             initStep = False,
             **kwargs):
        '''
        '''
        super(bismark_methylation_extractor, self).__init__(initStep)
        if upstream is None:
            self.setInput('bamInput', bamInput)
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('bamInput')[1])))
            else:
                self.setOutput('outputdir', outputdir)
                
            self.setParam('threads', threads)
            
            if paired:
                self.setParam('type', 'paired')
            else:
                self.setParam('type', 'single')
            
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
           
            self.setParam('type', Configure.getType())
            
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
        
        if self.getParam('type') == 'paired':
            other_params.update({'--paired': True})
        elif self.setParam('type') == 'single':
            other_params.update({'--single': True})
        else:
            commonError("Wrong data type, must be 'single' or 'paired'!")
            
        other_params.update({'--multicore': math.ceil(self.getParam('threads') / 3)})
        other_params.update({'--output': self.getOutput('outputdir')})
            
        self.setOutput('covOutput', [os.path.join(self.getOutput('outputdir'), os.path.splitext(os.path.basename(x))[0]) + '.bedGraph.gz.bismark.zero.cov' for x in self.getInput('bamInput')])
        self.setOutput('bedGraphOutput', [os.path.join(self.getOutput('outputdir'), os.path.splitext(os.path.basename(x))[0]) + '.bedGraph.gz' for x in self.getInput('bamInput')])
        self.setOutput('covgzOutput', [os.path.join(self.getOutput('outputdir'), os.path.splitext(os.path.basename(x))[0]) + '.bismark.cov.gz' for x in self.getInput('bamInput')])
        self.setOutput('reportOutput', [os.path.join(self.getOutput('outputdir'), os.path.splitext(os.path.basename(x))[0]) + '_splitting_report.txt' for x in self.getInput('bamInput')])
        
        all_cmd = []
        
        tmp_cmd = self.cmdCreate(["bismark_methylation_extractor", 
                                  self.getParam('other_params'),
                                  self.getInput('bamInput')])
        all_cmd.append(tmp_cmd)
        
        self.setParam('cmd', all_cmd)
        
        finishFlag = self.stepInit(upstream)
        
        self.excute(finishFlag)



