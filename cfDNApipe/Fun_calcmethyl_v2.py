# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 15:13:42 2020

@author: Jiaqi Huang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, calcMethylV2
from .Configure import Configure
import os
import math

__metaclass__ = type


class calculate_methyl(StepBase):
    def __init__(
        self,
        tbxInput=None,
        bedInput=None,
        outputdir=None,
        threads=1,
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for computing methylation level from indexed methylation coverage file.

        calculate_methyl(tbxInput=None, bedInput=None, outputdir=None, threads=1, stepNum=None, upstream=None,)
        {P}arameters:
            tbxInput: list, input indexed methylation coverage files.
            bedInput: str, bed file contains genome regions which will be computed for methylation level.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """
        super(calculate_methyl, self).__init__(stepNum, upstream)

        # set fastq input
        if (upstream is None) or (upstream is True):
            self.setInput("tbxInput", tbxInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "compress_methyl":
                self.setInput("tbxInput", upstream.getOutput("tbxOutput"))
            else:
                raise commonError("Parameter upstream must from compress_methyl.")

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(self.getInput("covgzInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set threads
        if upstream is None:
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure.getThreads())

        # set bedInput
        if bedInput is None:
            self.setInput("bedInput", Configure.getConfig("CpGisland"))
        else:
            self.setInput("bedInput", bedInput)

        self.setOutput(
            "txtOutput",
            [
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + "_Methy.txt"
                for x in self.getInput("tbxInput")
            ],
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            multi_run_len = len(self.getInput("tbxInput"))
            if verbose:
                for i in range(multi_run_len):
                    calcMethylV2(
                        tbxInput=self.getInput("tbxInput")[i],
                        bedInput=self.getInput("bedInput"),
                        txtOutput=self.getOutput("txtOutput")[i],
                    )
            else:
                multiArgs = [
                    [self.getInput("tbxInput")[i], self.getInput("bedInput"), self.getOutput("txtOutput")[i]]
                    for i in range(multi_run_len)
                ]
                self.multiRun(args=multiArgs, func=calcMethylV2, nCore=math.ceil(self.getParam("threads") / 4))

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
