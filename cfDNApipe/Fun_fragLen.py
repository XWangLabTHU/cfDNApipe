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
import math

__metaclass__ = type


class fraglenplot(StepBase):
    def __init__(
        self,
        bedInput=None,
        outputdir=None,
        maxLimit=500,
        threads=None,
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for ploting fragment length distribution.

        fraglenplot(bedInput=None, outputdir=None, maxLimit=500, threads=None, stepNum=None, upstream=None, verbose=True)
        {P}arameters:
            bedInput: list, input bed files.
            ref: bismark reference path.
            threads: int, how many thread to use.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline. This parameter can be True, which means a new pipeline start.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """
        super(fraglenplot, self).__init__(stepNum, upstream)

        # set fastq input
        if (upstream is None) or (upstream is True):
            self.setInput("bedInput", bedInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "bam2bed":
                self.setInput("bedInput", upstream.getOutput("bedOutput"))
            else:
                raise commonError("Parameter upstream must from bam2bed.")

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(self.getInput("bedInput")[1])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set ref, threads, paired
        if upstream is None:
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure.getThreads())

        # set maxLimit
        self.setParam("maxLimit", maxLimit)

        self.setOutput(
            "singleplotOutput",
            [
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + "_fraglen.png"
                for x in self.getInput("bedInput")
            ],
        )
        self.setOutput(
            "multiplotOutput", os.path.join(self.getOutput("outputdir"), "length_distribution.png"),
        )
        self.setOutput(
            "npyOutput",
            [
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + "_fraglen.npy"
                for x in self.getInput("bedInput")
            ],
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            multi_run_len = len(self.getInput("bedInput"))
            if verbose:
                len_data = []
                for i in range(multi_run_len):
                    print("Now, ploting fragment length distribution for " + self.getInput("bedInput")[i])
                    len_data.append(
                        fraglendistribution(
                            bedInput=self.getInput("bedInput")[i],
                            plotOutput=self.getOutput("singleplotOutput")[i],
                            binOutput=self.getOutput("npyOutput")[i],
                            maxLimit=self.getParam("maxLimit"),
                        )
                    )

                fraglenmultiplot(
                    dataInput=len_data, plotOutput=self.getOutput("multiplotOutput"),
                )
            else:
                args = [
                    [
                        self.getInput("bedInput")[i],
                        self.getOutput("singleplotOutput")[i],
                        self.getOutput("npyOutput")[i],
                        self.getParam("maxLimit"),
                    ]
                    for i in range(multi_run_len)
                ]
                len_data = self.multiRun(args=args, func=fraglendistribution, nCore=math.ceil(self.getParam("threads") / 4))
                fraglenmultiplot(
                    dataInput=len_data, plotOutput=self.getOutput("multiplotOutput"),
                )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
