# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 18:27:32 2019

@author: Jiqi Huang, Chang Li
"""


import os
from .StepBase import StepBase
from .Configure import Configure
from .cfDNA_utils import commonError, calcMethyl


__metaclass__ = type


class computemethyl(StepBase):
    def __init__(
        self,
        bamInput=None,  # list
        bedInput=None,
        outputdir=None,  # str
        stepNum=None,
        upstream=None,
        **kwargs
    ):
        super(computemethyl, self).__init__(stepNum, upstream)
        if upstream is None:
            self.setInput("bamInput", bamInput)
            if bedInput is None:
                self.setInput(
                    "bedInput",
                    os.path.join(
                        Configure.getRefDir(),
                        Configure.getGenome() + "cpgIslandExt.bed",
                    ),
                )
            else:
                self.setInput("bedInput", bedInput)
            self.checkInputFilePath()

            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "rmduplicate":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
                if bedInput is None:
                    self.setInput(
                        "bedInput",
                        os.path.join(
                            Configure.getRefDir(),
                            Configure.getGenome() + "cpgIslandExt.bed",
                        ),
                    )
                else:
                    self.setInput("bedInput", bedInput)
            else:
                raise commonError("Parameter upstream must from rmduplicate.")

            self.setOutput("outputdir", self.getStepFolderPath())

        self.setOutput(
            "txtOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + ".txt"
                for x in self.getInput("bamInput")
            ],
        )

        finishFlag = self.stepInit(upstream)

        if finishFlag:
            self.excute(finishFlag)
        else:
            multi_run_len = len(self.getInput("bamInput"))

            for i in range(multi_run_len):
                print("Now, processing file: " + self.getInput("bamInput")[i])
                calcMethyl(
                    bamInput=self.getInput("bamInput")[i],
                    bedInput=self.getInput("bedInput"),
                    txtOutput=self.getOutput("txtOutput")[i],
                )

            self.excute(finishFlag, runFlag=False)
