# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:51:54 2019

@author: zhang
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError, computeCUE, OCFplot
import os
import pybedtools
from .Configure import Configure


__metaclass__ = type


class computeOCF(StepBase):
    def __init__(self, 
         casebedInput = None, # list
         ctrlbedInput = None,
         refRegInput = None,
         outputdir = None, # str
         caseupstream = None,
         ctrlupstream = None,
         labelInput = None,
         initStep = False,
         **kwargs):
        super(computeOCF, self).__init__(initStep)
        labelflag = False
        if caseupstream is None and ctrlupstream is None:
            self.setInput('casebedInput', casebedInput)
            self.setInput('ctrlbedInput', ctrlbedInput)
            self.setInput('refRegInput', refRegInput)
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('bedInput')[0])))
            else:
                self.setOutput('outputdir', outputdir)
                
        else:
            Configure.configureCheck()
            caseupstream.checkFilePath()
            ctrlupstream.checkFilePath()
            
            if caseupstream.__class__.__name__ == 'bam2bed':
                self.setInput('casebedInput', caseupstream.getOutput("bedOutput"))
            else:
                raise commonError('Case Parameter upstream must from bam2bed.')
            if ctrlupstream.__class__.__name__ == 'bam2bed':
                self.setInput('ctrlbedInput', ctrlupstream.getOutput("bedOutput"))
            else:
                raise commonError('Control Parameter upstream must from bam2bed.')
            
            self.setInput('refRegInput', refRegInput)
            self.setOutput('outputdir', self.getStepFolderPath())
        
        if labelInput is not None:
            self.setInput('labelInput', labelInput)
            labelflag = True
        
        self.setOutput('casetxtOutput', [os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x)) + '.txt' for x in self.getInput('casebedInput')])
        self.setOutput('ctrltxtOutput', [os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x)) + '.txt' for x in self.getInput('ctrlbedInput')])
       
        save_flag = ["Tcell", "Liver", "Placenta", "Lung", "Breast", "Intestine", "Ovary"]
        casecudOutput = []
        ctrlcudOutput = []
        for x in self.getInput('casebedInput'):
            prefix = os.path.splitext(os.path.basename(x))[0]
            for flag in save_flag:
                casecudOutput.append(self.getOutput('outputdir') + '/' + prefix + '-' + flag + '-cud.txt')
        for x in self.getInput('ctrlbedInput'):
            prefix = os.path.splitext(os.path.basename(x))[0]
            for flag in save_flag:
                ctrlcudOutput.append(self.getOutput('outputdir') + '/' + prefix + '-' + flag + '-cud.txt')
        self.setOutput('casecudOutput', casecudOutput)
        self.setOutput('ctrlcudOutput', ctrlcudOutput)
        
        self.setOutput('ocfOutput', [os.path.join(self.getOutput('outputdir'), 'OCF-case.txt'), os.path.join(self.getOutput('outputdir'), 'OCF-control.txt')])
        self.setOutput('plotOutput', [os.path.join(self.getOutput('outputdir'), 'OCF-') + flag + '.png' for flag in save_flag])
        
        finishFlag = self.stepInit(caseupstream) #need to be checked
        
        if finishFlag:
            self.excute(finishFlag)
        else:
            case_multi_run_len = len(self.getInput('casebedInput'))
            ctrl_multi_run_len = len(self.getInput('ctrlbedInput'))
            ocf_case = [[] for i in range(case_multi_run_len)]
            ocf_ctrl = [[] for i in range(ctrl_multi_run_len)]
            case_fp = open(self.getOutput('ocfOutput')[0], 'w+')
            ctrl_fp = open(self.getOutput('ocfOutput')[1], 'w+')
            for i in range(case_multi_run_len):
                print("Now, processing file: " + self.getInput('casebedInput')[i])
                ocf_case[i] = computeCUE(inputFile = self.getInput('casebedInput')[i], refFile = self.getInput('refRegInput'), txtOutput = self.getOutput('casetxtOutput')[i], cudOutput = self.getOutput('casecudOutput')[7 * i : 7 * i + 7])
                for ocfvalue in ocf_case[i]:
                    case_fp.write(str(ocfvalue) + '\t')
                case_fp.write('\n')
            case_fp.close()
            for i in range(ctrl_multi_run_len):
                print("Now, processing file: " + self.getInput('ctrlbedInput')[i])
                ocf_ctrl[i] = computeCUE(inputFile = self.getInput('ctrlbedInput')[i], refFile = self.getInput('refRegInput'), txtOutput = self.getOutput('ctrltxtOutput')[i], cudOutput = self.getOutput('ctrlcudOutput')[7 * i : 7 * i + 7])
                for ocfvalue in ocf_ctrl[i]:
                    ctrl_fp.write(str(ocfvalue) + '\t')
                ctrl_fp.write('\n')
            ctrl_fp.close()
            if labelflag:
                OCFplot(ocf_case, ocf_ctrl, self.getOutput('plotOutput'), self.getInput('labelInput'))
            else:
                OCFplot(ocf_case, ocf_ctrl, self.getOutput('plotOutput'))
            
            self.excute(finishFlag, runFlag = False)
            
            
            
            
            