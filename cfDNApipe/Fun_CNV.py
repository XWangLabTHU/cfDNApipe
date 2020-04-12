# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:51:54 2019

@author: zhang, Huang
"""

from .StepBase2 import StepBase2
from .cfDNA_utils import commonError, sumChromarm, plotCNVheatmap
import pandas as pd
import numpy as np
import os
from .Configure2 import Configure2

__metaclass__ = type


class computeCNV(StepBase2):
    def __init__(
        self,
        casetxtInput=None,
        ctrltxtInput=None,
        outputdir=None,
        cytoBandInput=None,
        stepNum=None,
        caseupstream=None,
        ctrlupstream=None,
        **kwargs
    ):
        """
        This function is used for computing CNV z-scores.

        computeCNV(casetxtInput=None, ctrltxtInput=None, outputdir=None, cytoBandInput=None, stepNum=None, caseupstream=None, ctrlupstream=None,)
        {P}arameters:
            casetxtInput: list, paths of files of GC corrected read counts for case samples.
            ctrltxtInput: list, paths of files of GC corrected read counts for control samples.
            outputdir: str, output result folder, None means the same folder as input files.
            cytoBandInput: str, path of the cytoBand file.
            stepNum: Step number for folder name.
            caseupstream: Not used parameter, do not set this parameter.
            ctrlupstream: Not used parameter, do not set this parameter.
        """

        if (stepNum is None) and (caseupstream is not None) and (ctrlupstream is None):
            super(computeCNV, self).__init__(stepNum, caseupstream)
        elif (stepNum is None) and (caseupstream is None) and (ctrlupstream is not None):
            super(computeCNV, self).__init__(stepNum, ctrlupstream)
        elif (stepNum is None) and (caseupstream is not None) and (ctrlupstream is not None):
            if caseupstream.getStepID() >= ctrlupstream.getStepID():
                super(computeCNV, self).__init__(stepNum, caseupstream)
            else:
                super(computeCNV, self).__init__(stepNum, ctrlupstream)
        else:
            super(computeCNV, self).__init__(stepNum)

        if caseupstream is None and ctrlupstream is None:
            self.setInput("casetxtInput", casetxtInput)
            self.setInput("ctrltxtInput", ctrltxtInput)
            self.checkInputFilePath()

            if outputdir is None:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(self.getInput("casetxtInput")[1])),
                )
            else:
                self.setOutput("outputdir", outputdir)

        else:
            Configure2.configureCheck()
            caseupstream.checkFilePath()
            ctrlupstream.checkFilePath()

            if caseupstream.__class__.__name__ == "GCCorrect":
                self.setInput("casetxtInput", caseupstream.getOutput("txtOutput"))
            else:
                raise commonError("Parameter upstream must from GCCorrect.")
            if ctrlupstream.__class__.__name__ == "GCCorrect":
                self.setInput("ctrltxtInput", ctrlupstream.getOutput("txtOutput"))
            else:
                raise commonError("Parameter upstream must from GCCorrect.")
            self.checkInputFilePath()

            self.setOutput("outputdir", self.getStepFolderPath())

        if cytoBandInput is None:
            self.setInput("cytoBandInput", Configure2.getConfig("cytoBand"))
        else:
            self.setInput("cytoBandInput", cytoBandInput)

        self.setOutput("txtOutput", self.getOutput("outputdir") + "/Z-score.txt")

        self.setOutput("plotOutput", self.getOutput("outputdir") + "/CNV.png")

        finishFlag = self.stepInit(caseupstream)

        case_multi_run_len = len(self.getInput("casetxtInput"))
        ctrl_multi_run_len = len(self.getInput("ctrltxtInput"))

        if not finishFlag:
            # process CNV
            case_chrom = [[] for i in range(case_multi_run_len)]
            ctrl_chrom = [[] for i in range(ctrl_multi_run_len)]
            genes = []  # stores the names of the chromosome arms
            for i in range(case_multi_run_len):
                case_chrom[i], genes = sumChromarm(self.getInput("casetxtInput")[i], self.getInput("cytoBandInput"))

            # build the case chromarms-divided dataframe, columns are case samples, rows are chromosome arms
            case_chrom_df = pd.DataFrame(
                np.transpose(case_chrom),
                columns=[self.getMaxFileNamePrefixV2(x).split(".")[0] for x in self.getInput("casetxtInput")],
                index=genes,
            )

            for i in range(ctrl_multi_run_len):
                ctrl_chrom[i], genes = sumChromarm(self.getInput("ctrltxtInput")[i], self.getInput("cytoBandInput"))

            ctrl_chrom_df = pd.DataFrame(
                np.transpose(ctrl_chrom),
                columns=[self.getMaxFileNamePrefixV2(x).split(".")[0] for x in self.getInput("ctrltxtInput")],
                index=genes,
            )

            # compute z-score and plot heatmap
            plotCNVheatmap(
                case_chrom_df, ctrl_chrom_df, self.getOutput("txtOutput"), self.getOutput("plotOutput"),
            )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
