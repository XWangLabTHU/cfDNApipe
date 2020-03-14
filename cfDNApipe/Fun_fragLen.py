# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 10:51:10 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, fraglendistribution, fraglenmultiplot
import os
from .Configure import Configure

__metaclass__ = type


class fraglenplot(StepBase):
    def __init__(
            self,
            bedInput=None,  # list
            outputdir=None,  # str
            maxLimit=500,
            stepNum=None,
            upstream=None,
            **kwargs):
        super(fraglenplot, self).__init__(stepNum, upstream)
        if upstream is None:
            self.setInput("bedInput", bedInput)
            self.checkInputFilePath()

            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(
                        os.path.abspath(self.getInput("bedInput")[1])),
                )
            else:
                self.setOutput("outputdir", outputdir)

        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "bam2bed":
                self.setInput("bedInput", upstream.getOutput("bedOutput"))
            else:
                raise commonError("Parameter upstream must from bam2bed.")

            self.setOutput("outputdir", self.getStepFolderPath())

        self.setParam("maxLimit", maxLimit)
        self.setOutput(
            "singleplotOutput",
            [
                os.path.join(self.getOutput("outputdir"),
                             self.getMaxFileNamePrefixV2(x)) + "_fraglen.png"
                for x in self.getInput("bedInput")
            ],
        )
        self.setOutput(
            "multiplotOutput",
            self.getOutput("outputdir") + "/" + "length_distribution.png",
        )
        self.setOutput(
            "npyOutput",
            [
                os.path.join(self.getOutput("outputdir"),
                             self.getMaxFileNamePrefixV2(x)) + "_fraglen.npy"
                for x in self.getInput("bedInput")
            ],
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            multi_run_len = len(self.getInput("bedInput"))
            len_data = []
            for i in range(multi_run_len):
                print("Now, ploting fragment length distribution for " +
                      self.getInput("bedInput")[i])
                len_data.append(
                    fraglendistribution(
                        bedInput=self.getInput("bedInput")[i],
                        plotOutput=self.getOutput("singleplotOutput")[i],
                        binOutput=self.getOutput("npyOutput")[i],
                        maxLimit=self.getParam("maxLimit"),
                    ))

            fraglenmultiplot(
                dataInput=len_data,
                plotOutput=self.getOutput("multiplotOutput"),
            )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
