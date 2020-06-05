# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 15:11:11 2020

@author: Jiaqi Huang
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError, compressMethy
from .Configure import Configure
import os
import math


__metaclass__ = type


class compress_methyl(StepBase):
    def __init__(
        self,
        covInput=None,
        outputdir=None,
        threads=1,
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for compressing and fast indexing methlation information from bismark_methylation_extractor.

        compress_methyl(covInput=None, outputdir=None, stepNum=None, upstream=None,)
        {P}arameters:
            covInput: list, input methylation coverage files.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """

        super(compress_methyl, self).__init__(stepNum, upstream)

        # set fastq input
        if (upstream is None) or (upstream is True):
            self.setInput("covInput", covInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "bismark_methylation_extractor":
                self.setInput("covInput", upstream.getOutput("covOutput"))
            else:
                raise commonError(
                    "Parameter upstream must from bismark_methylation_extractor."
                )

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("covInput")[0])),
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

        self.setOutput(
            "tbxOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + ".gz"
                for x in self.getInput("covInput")
            ],
        )
        self.setOutput(
            "tbiOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + ".gz.tbi"
                for x in self.getInput("covInput")
            ],
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            multi_run_len = len(self.getInput("covInput"))
            if verbose:
                for i in range(multi_run_len):
                    compressMethy(
                        InputFile=self.getInput("covInput")[i],
                        OutputFile=self.getOutput("tbxOutput")[i],
                    )
            else:
                args = [
                    [self.getInput("covInput")[i], self.getOutput("tbxOutput")[i]]
                    for i in range(multi_run_len)
                ]
                self.multiRun(
                    args=args,
                    func=compressMethy,
                    nCore=math.ceil(self.getParam("threads") / 4),
                )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
