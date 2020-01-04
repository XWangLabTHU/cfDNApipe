# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
@author: LY
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure


__metaclass__ = type


class addRG(StepBase):
    def __init__(self,
                 addRGInput = None,
                 addRGOutputDir = None,
                 upstream = None,
                 initStep = False,
                 **kwargs):
        super(addRG, self).__init__(initStep)
        if upstream is None:
            self.setInput('addRGInput', addRGInput)
            self.checkInputFilePath()
            
            if addRGOutputDir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('addRGInput')[1])))
            else:
                self.setOutput('outputdir', addRGOutputDir)

        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ in ['rmduplicate', 'bismark', 'bismark_deduplicate', 'bowtie2']:
                self.setInput('addRGInput', upstream.getOutput('bamOutput'))
            else:
                raise commonError('Parameter upstream must from rmduplicate.')
                
            self.setOutput('outputdir', self.getStepFolderPath())

        self.setOutput('addRGOutput', [os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x)) + '-RG.bam' for x in self.getInput('addRGInput')])
        self.setParam('RGID', [self.getMaxFileNamePrefixV2(x) for x in self.getInput('addRGInput')])
        self.setParam('RGLB', [self.getMaxFileNamePrefixV2(x) for x in self.getInput('addRGInput')])
        self.setParam('RGPL', [self.getMaxFileNamePrefixV2(x) for x in self.getInput('addRGInput')])
        self.setParam('RGSM', [self.getMaxFileNamePrefixV2(x) for x in self.getInput('addRGInput')])
        self.setParam('RGPU', [self.getMaxFileNamePrefixV2(x) for x in self.getInput('addRGInput')])
        
        multi_run_len = len(self.getInput('addRGInput'))

        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(['picard AddOrReplaceReadGroups',
                                      'I=' + self.getInput('addRGInput')[i],
                                      'O=' + self.getOutput('addRGOutput')[i],
                                      'SORT_ORDER=coordinate',
                                      'RGID=' + self.getParam('RGID')[i],
                                      'RGLB=' + self.getParam('RGLB')[i],
                                      'RGPL=' + self.getParam('RGPL')[i],
                                      'RGSM=' + self.getParam('RGSM')[i],
                                      'RGPU=' + self.getParam('RGPU')[i],
                                      'CREATE_INDEX=True'])
            all_cmd.append(tmp_cmd)

        self.setParam('cmd', all_cmd)

        finishFlag = self.stepInit(upstream)

        self.excute(finishFlag)





# Functions to get file name prefix?
def GetAddRGName(str):
    addRGbegin = (str.rfind('/'))
    addRGend1 = (str[addRGbegin:].find('_'))
    addRGend2 = (str[addRGbegin:].find('.'))
    if addRGbegin < 0:
        raise commonError('Sample naming error')

    elif addRGend1 > 0:
        return (str[addRGbegin + 1:addRGbegin + addRGend1])
    elif addRGend2 > 0:
        return (str[addRGbegin + 1:addRGbegin + addRGend2])
    else:
        raise commonError('Sample naming error')
