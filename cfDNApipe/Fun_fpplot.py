# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 14:21:15 2020

@author: Jiaqi Huang

"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, fragProfileplot
from .Configure2 import Configure2
import os


__metaclass__ = type


class fragprofplot(StepBase):
    def __init__(
        self,
        casetxtInput=None,
        ctrltxtInput=None,
        cytoBandInput=None,
        labelInput=None,
        outputdir=None,
        stepNum=None,
        caseupstream=None,
        ctrlupstream=None,
        **kwargs
    ):
        """
        This function is used for computing and plotting the fragmentation profile plot.

        fragprofplot(casetxtInput=None, ctrltxtInput=None, outputdir=None, cytoBandInput=None, stepNum=None, caseupstream=None,
        ctrlupstream=None)
        {P}arameters:
            casetxtInput: list, paths of files of input short and long counts for case samples,
                          [sample1_short, sample1_long, sample2_short, sample2_long, etc.].
            ctrltxtInput: list, paths of files of input short and long counts for control samples(format is same as casetxtInput).
            cytoBandInput: str, path of the cytoBand file.
            labelInput: list, name of case and control(which will be marked on the output plot), [case_name, control_name]
            outputdir: str, output result folder, None means the same folder as input files.
            stepNum: Step number for folder name.
            caseupstream: Not used parameter, do not set this parameter.
            ctrlupstream: Not used parameter, do not set this parameter.
        """

        if (stepNum is None) and (caseupstream is not None) and (ctrlupstream is None):
            super(fragprofplot, self).__init__(stepNum, caseupstream)
        elif (
            (stepNum is None) and (caseupstream is None) and (ctrlupstream is not None)
        ):
            super(fragprofplot, self).__init__(stepNum, ctrlupstream)
        elif (
            (stepNum is None)
            and (caseupstream is not None)
            and (ctrlupstream is not None)
        ):
            if caseupstream.getStepID() >= ctrlupstream.getStepID():
                super(fragprofplot, self).__init__(stepNum, caseupstream)
            else:
                super(fragprofplot, self).__init__(stepNum, ctrlupstream)
        else:
            super(fragprofplot, self).__init__(stepNum)

        # set casetxtInput and ctrltxtInput
        if (
            ((caseupstream is None) and (ctrlupstream is None))
            or (caseupstream is True)
            or (ctrlupstream is True)
        ):
            self.setInput("casetxtInput", casetxtInput)
            self.setInput("ctrltxtInput", ctrltxtInput)
        else:
            Configure2.configureCheck()
            caseupstream.checkFilePath()
            ctrlupstream.checkFilePath()
            if (
                caseupstream.__class__.__name__ == "fpCounter"
                or caseupstream.__class__.__name__ == "GCCorrect"
            ):
                self.setInput("casetxtInput", caseupstream.getOutput("txtOutput"))
            else:
                raise commonError(
                    "Parameter upstream must from fpCounter or GCCorrect."
                )

            if (
                ctrlupstream.__class__.__name__ == "fpCounter"
                or ctrlupstream.__class__.__name__ == "GCCorrect"
            ):
                self.setInput("ctrltxtInput", ctrlupstream.getOutput("txtOutput"))
            else:
                raise commonError(
                    "Parameter upstream must from fpCounter or GCCorrect."
                )

        self.checkInputFilePath()

        # set outputdir
        if (caseupstream is None) and (ctrlupstream is None):
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("casetxtInput")[1])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        if cytoBandInput is None:
            self.setInput("cytoBandInput", Configure2.getConfig("cytoBand"))
        else:
            self.setInput("cytoBandInput", cytoBandInput)

        if labelInput is not None:
            self.setParam("label", labelInput)
        else:
            self.setParam("label", ["case", "control"])

        self.setOutput(
            "plotOutput",
            os.path.join(
                self.getOutput("outputdir"), "fragmentation_profile_plot.png",
            ),
        )
        
        self.setOutput(
            "txtOutput",
            os.path.join(
                self.getOutput("outputdir"), "fragmentation_profile.txt",
            ),
        )

        finishFlag = self.stepInit(caseupstream)

        if not finishFlag:
            fragProfileplot(
                self.getInput("casetxtInput"),
                self.getInput("ctrltxtInput"),
                self.getInput("cytoBandInput"),
                self.getOutput("plotOutput"),
                self.getOutput("txtOutput"),
                self.getParam("label"),
            )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
