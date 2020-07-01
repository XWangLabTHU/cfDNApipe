# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:51:54 2019

@author: zhang, Huang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, OCF_boxplot
import os
from .Configure2 import Configure2

__metaclass__ = type


class OCFplot(StepBase):
    def __init__(
        self,
        caseocfInput=None,
        ctrlocfInput=None,
        outputdir=None,
        threads=1,
        upstream=None,
        labelInput=None,
        stepNum=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for plotting the boxplot of OCF values from computeOCF.

        OCFplot(caseocfInput=None, ctrlocfInput=None, outputdir=None, threads=1, upstream=None, labelInput=None, stepNum=None,  verbose=True)
        {P}arameters:
            caseocfInput: list, input files of OCF values of case samples.
            ctrlocfInput: list, input files of OCF values of control samples.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            upstream: upstream output results, used for pipeline.
            labelInput: list, [name_of_case, name_of_control](e.g. ["HCC", "CTR"])
            stepNum: int, step number for folder name.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """
        super(OCFplot, self).__init__(stepNum, upstream)

        # set caseocfInput and ctrlocfInput
        if (upstream is None) or (upstream is True):
            self.setInput("caseocfInput", caseocfInput)
            self.setInput("ctrlocfInput", ctrlocfInput)
        else:
            Configure2.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "computeOCF":
                self.setInput("caseocfInput", upstream.getOutput("caseocfOutput"))
                self.setInput("ctrlocfInput", upstream.getOutput("ctrlocfOutput"))
            else:
                raise commonError("Parameter upstream must from computeOCF.")

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("caseocfInput")[1])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set labelInput
        if labelInput is not None:
            self.setParam("label", labelInput)
        else:
            self.setParam("label", ["case", "control"])

        self.setOutput(
            "plotOutput", os.path.join(self.getOutput("outputdir"), "OCF-boxplot.png")
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            OCF_boxplot(
                self.getInput("caseocfInput"),
                self.getInput("ctrlocfInput"),
                self.getOutput("plotOutput"),
                self.getParam("label"),
            )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
