# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 18:27:32 2019

@author: Jiaqi Huang, Chang Li
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError
from .Configure import Configure
import os


__metaclass__ = type


class sequencetransfer(StepBase):
    def __init__(self, bamInput=None, bedInput=None, outputdir=None, threads=1, stepNum=None, upstream=None, **kwargs):

        super(sequencetransfer, self).__init__(stepNum, upstream)

        # set bamInput
        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "rmduplicate":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from inputprocess or adapterremoval.")

        self.checkInputFilePath()

        bedflag = False

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set bedInput, threads
        if upstream is None:
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure.getThreads())

        if bedInput is not None:
            bedflag = True
            self.setInput("bedInput", bedInput)
        else:
            commonError("Parameter bedInput is None!")

        self.setOutput(
            "txtOutput",
            [
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + "-seq.txt"
                for x in self.getInput("bamInput")
            ],
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            multi_run_len = len(self.getInput("bamInput"))
            for i in range(multi_run_len):
                if bedflag:
                    self.seqTrans(
                        bamInput=self.getInput("bamInput")[i],
                        txtOutput=self.getOutput("txtOutput")[i],
                        bedInput=self.getInput("bedInput"),
                    )
                else:
                    self.seqTrans(
                        bamInput=self.getInput("bamInput")[i], txtOutput=self.getOutput("txtOutput")[i],
                    )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
