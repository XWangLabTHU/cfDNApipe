# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:51:54 2019

@author: zhang, Huang
"""

from .StepBase2 import StepBase2
from .cfDNA_utils import commonError, wig2df, correctReadCount, chromarm_sum, compute_z_score
import os
import pandas as pd
import numpy as np
from .Configure import Configure

__metaclass__ = type


class computeCNV(StepBase2):
    def __init__(
            self,
            casereadInput=None,  # list
            casegcInput=None,
            ctrlreadInput=None,   # list
            ctrlgcInput=None,
            outputdir=None,  # str
            caseupstream=None,
            ctrlupstream=None,
            labelInput=None,
            stepNum=None,
            **kwargs):
        if (stepNum is None) and (caseupstream is not None) and (ctrlupstream is None):
            super(computeCNV, self).__init__(stepNum, caseupstream)
        elif ((stepNum is None) and (caseupstream is None) and (ctrlupstream is not None)):
            super(computeCNV, self).__init__(stepNum, ctrlupstream)
        elif ((stepNum is None) and (caseupstream is not None) and (ctrlupstream is not None)):
            if caseupstream.getStepID() >= ctrlupstream.getStepID():
                super(computeCNV, self).__init__(stepNum, caseupstream)
            else:
                super(computeCNV, self).__init__(stepNum, ctrlupstream)
        else:
            super(computeCNV, self).__init__(stepNum)
        
        if labelInput is not None:
            labels = labelInput
        else:
            labels = ["case", "control"]
        
        if caseupstream is None and ctrlupstream is None:
            self.setInput("casereadInput", casereadInput)
            self.setInput("casegcInput", casegcInput)
            self.setInput("ctrlreadInput", ctrlreadInput)
            self.setInput("ctrlgcInput", ctrlgcInput)
            self.checkInputFilePath()

            if outputdir is None:
                commonError("Parameter 'outputdir' cannot be None!")
            else:
                self.setOutput("outputdir", outputdir)

        else:
            Configure.configureCheck()
            caseupstream.checkFilePath()
            ctrlupstream.checkFilePath()

            if caseupstream.__class__.__name__ == "readCount":
                self.setInput("casereadInput",
                              caseupstream.getOutput("readOutput"))
                self.setInput("casegcInput",
                              caseupstream.getOutput("gcOutput"))
            else:
                raise commonError("Parameter caseupstream must from readCount.")
            if ctrlupstream.__class__.__name__ == "readCount":
                self.setInput("ctrlreadInput",
                              ctrlupstream.getOutput("readOutput"))
                self.setInput("ctrlgcInput",
                              ctrlupstream.getOutput("gcOutput"))
            else:
                raise commonError(
                    "Parameter ctrlupstream must from readCount.")

            self.setOutput("outputdir", self.getStepFolderPath())
            
        self.setInput("cytoBandInput", Configure.getConfig('cytoBand'))
        
        self.setOutput(
            "txtOutput",
            self.getOutput("outputdir") + "/Z-score.txt")
        
        self.setOutput(
            "casereadplotOutput",
            [self.getOutput("outputdir") + '/' + self.getMaxFileNamePrefixV2(x) + ".png" for x in self.getInput("casereadInput")]
        )
        
        self.setOutput(
            "ctrlreadplotOutput",
            [self.getOutput("outputdir") + '/' + self.getMaxFileNamePrefixV2(x) + ".png" for x in self.getInput("ctrlreadInput")]
        )
        
        self.setOutput(
            "plotOutput", 
            self.getOutput("outputdir") + "/CNV.png")

        finishFlag = self.stepInit(caseupstream)

        if finishFlag:
            self.excute(finishFlag)
        else:
            case_multi_run_len = len(self.getInput("casereadInput"))
            ctrl_multi_run_len = len(self.getInput("ctrlreadInput"))
            case_chrom = [[] for i in range(case_multi_run_len)]
            ctrl_chrom = [[] for i in range(ctrl_multi_run_len)]
            genes = []
            case_df_gc = wig2df(self.getInput("casegcInput"))
            for i in range(case_multi_run_len):
                case_df_read = wig2df(self.getInput("casereadInput")[i])
                case_read_correct = correctReadCount(case_df_read, case_df_gc, self.getOutput("casereadplotOutput")[i])
                case_chrom[i], genes = chromarm_sum(case_read_correct, self.getInput("cytoBandInput"))
            case_df = pd.DataFrame(
                np.transpose(case_chrom), 
                columns = [self.getMaxFileNamePrefixV2(x)[ : -5] for x in self.getInput("casereadInput")],
                index = genes
            )
            
            ctrl_df_gc = wig2df(self.getInput("ctrlgcInput"))
            for i in range(ctrl_multi_run_len):
                ctrl_df_read = wig2df(self.getInput("ctrlreadInput")[i])
                ctrl_read_correct = correctReadCount(ctrl_df_read, ctrl_df_gc, self.getOutput("ctrlreadplotOutput")[i])
                ctrl_chrom[i], genes = chromarm_sum(ctrl_read_correct, self.getInput("cytoBandInput"))
            ctrl_df = pd.DataFrame(
                np.transpose(ctrl_chrom), 
                columns = [self.getMaxFileNamePrefixV2(x)[ : -5] for x in self.getInput("ctrlreadInput")],
                index = genes
            )
                
            compute_z_score(case_df, ctrl_df, self.getOutput("txtOutput"), self.getOutput("plotOutput"))

            self.excute(finishFlag, runFlag=False)
