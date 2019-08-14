# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:26:39 2019

@author: zhang
"""

from Configure import Configure
from StepBase import StepBase
from cfDNA_utils import commonError
import os


__metaclass__ = type


class inputprocess(StepBase):
    def __init__(self,
                 fqInput1 = None,
                 fqInput2 = None,
                 inputFolder = None,
                 paired = True):
        super(inputprocess, self).__init__()
        
        # check Configure for running pipeline
        Configure.configureCheck()
        
        if not paired:
            raise commonError("Not support single end now, we will update as quick as possible.")
        
        if inputFolder is not None:  # using folder first
            all_files = os.listdir(inputFolder)
            all_files.sort()
            all_files = list(map(lambda x: os.path.join(inputFolder, x), all_files))
            if paired:
                fqInput1 = []
                fqInput2 = []
                for i in range(len(all_files)):
                    if i % 2:
                        fqInput2.append(all_files[i])
                    else:
                        fqInput1.append(all_files[i])
                self.setInput('fq1', fqInput1)
                self.setInput('fq2', fqInput2)
        else:
            if fqInput2 is None:
                raise commonError("Not support single end now, we will update as quick as possible.")
            else:
                self.setInput('fq1', fqInput1)
                self.setInput('fq2', fqInput2)
        
        self.checkInputFilePath()
        
        self.setOutput('fq1', self.getInput('fq1'))
        self.setOutput('fq2', self.getInput('fq2'))
        self.setOutput('outputdir', self.getStepFolderPath())
        
        finishFlag = self.stepInit(upstream = True)
        
        self.excute(finishFlag, runFlag = False)
        
        






















