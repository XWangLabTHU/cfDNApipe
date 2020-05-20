# -*- coding: utf-8 -*-
"""
Created on Sun Mar 8 17:15:15 2020

@author: Jiaqi Huang
"""

from .StepBase2 import StepBase2
from .cfDNA_utils import commonError, processPCA, clusterplot
import os
from .Configure2 import Configure2

__metaclass__ = type


class PCAplot(StepBase2):
    def __init__(
        self,
        casetxtInput=None,
        ctrltxtInput=None,
        outputdir=None,
        caseupstream=None,
        ctrlupstream=None,
        labelInput=None,
        stepNum=None,
        **kwargs
    ):
        if (stepNum is None) and (caseupstream is not None) and (ctrlupstream is None):
            super(PCAplot, self).__init__(stepNum, caseupstream)
        elif (stepNum is None) and (caseupstream is None) and (ctrlupstream is not None):
            super(PCAplot, self).__init__(stepNum, ctrlupstream)
        elif (stepNum is None) and (caseupstream is not None) and (ctrlupstream is not None):
            if caseupstream.getStepID() >= ctrlupstream.getStepID():
                super(PCAplot, self).__init__(stepNum, caseupstream)
            else:
                super(PCAplot, self).__init__(stepNum, ctrlupstream)
        else:
            super(PCAplot, self).__init__(stepNum)

        labelflag = False

        # set casebedInput and ctrlbedInput
        if ((caseupstream is None) and (ctrlupstream is None)) or (caseupstream is True) or (ctrlupstream is True):
            self.setInput("casetxtInput", casetxtInput)
            self.setInput("ctrltxtInput", ctrltxtInput)
        else:
            Configure2.configureCheck()
            caseupstream.checkFilePath()
            ctrlupstream.checkFilePath()
            if caseupstream.__class__.__name__ == "calculate_methyl":
                self.setInput("casetxtInput", caseupstream.getOutput("txtOutput"))
            else:
                raise commonError("Parameter caseupstream must from calculate_methyl.")
            if ctrlupstream.__class__.__name__ == "calculate_methyl":
                self.setInput("ctrltxtInput", ctrlupstream.getOutput("txtOutput"))
            else:
                raise commonError("Parameter ctrlupstream must from calculate_methyl.")

        self.checkInputFilePath()

        # set outputdir
        if (caseupstream is None) and (ctrlupstream is None):
            if outputdir is None:
                commonError("Parameter 'outputdir' cannot be None!")
            else:
                self.setOutput("outputdir", outputdir)

        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set labelInput
        if labelInput is not None:
            self.setParam("label", labelInput)
            labelflag = True

        self.setOutput("plotOutput", os.path.join(self.getOutput("outputdir"), "cluster_map.png"))

        finishFlag = self.stepInit(caseupstream)  # need to be checked

        if not finishFlag:
            casedata, ctrldata = processPCA(self.getInput("casetxtInput"), self.getInput("ctrltxtInput"))
            if labelflag:
                clusterplot(
                    casedata, ctrldata, self.getOutput("plotOutput"), self.getParam("label"),
                )
            else:
                clusterplot(
                    casedata, ctrldata, self.getOutput("plotOutput"),
                )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
