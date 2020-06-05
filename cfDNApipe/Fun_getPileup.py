# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
Modify on Sun Apr 26 11:27:32 2020
@author: LY, Shuying He
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import math

__metaclass__ = type


class getPileup(StepBase):
    def __init__(
        self,
        bamInput=None,
        biallelicvcfInput=None,
        outputdir=None,
        stepNum=None,
        upstream=None,
        threads=1,
        verbose=False,
        **kwargs
    ):

        super(getPileup, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "addRG":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            elif upstream.__class__.__name__ == "BQSR":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from addRG or BQSR.")
        self.checkInputFilePath()

        if upstream is None:
            self.setParam("threads", threads)
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())

        self.setInput("biallelicvcfInput", biallelicvcfInput)
        self.setOutput(
            "getPileupOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + ".pileups.table"
                for x in self.getInput("bamInput")
            ],
        )

        self.setOutput("bamOutput", self.getInput("bamInput"))

        multi_run_len = len(self.getInput("bamInput"))
        all_cmd = []
        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(
                [
                    "gatk GetPileupSummaries",
                    "-I",
                    self.getInput("bamInput")[i],
                    "-V",
                    self.getInput("biallelicvcfInput"),
                    "-L",
                    self.getInput("biallelicvcfInput"),
                    "-O",
                    self.getOutput("getPileupOutput")[i],
                ]
            )
            all_cmd.append(tmp_cmd)

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            if verbose:
                self.run(all_cmd)
            else:
                self.multiRun(
                    args=all_cmd,
                    func=None,
                    nCore=math.ceil(self.getParam("threads") / 4),
                )

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)
