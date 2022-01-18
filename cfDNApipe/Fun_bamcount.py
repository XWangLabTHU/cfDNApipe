# -*- coding: utf-8 -*-
"""
Created on Wed Apr 8 12:51:24 2020

@author: Jiaqi Huang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, count_bam, maxCore, divide_bin_2
import os
import math
from .Configure import Configure

__metaclass__ = type


class bamCounter(StepBase):
    def __init__(
        self,
        bamInput=None,
        chromsizeInput=None,
        binlen=None,
        outputdir=None,
        threads=1,
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for counting counts in each dividing bin of the input bam files.

        bamCounter(bamInput=None, chromsizeInput=None, binlen=None, outputdir=None, threads=1, stepNum=None, upstream=None, verbose=True,)
        {P}arameters:
            bamInput: list, paths of input bedgz files waiting to be processed.
            chromsizeInput: str, path of chromsize file.
            binlen: int, length of each bin; default is 5000000(5Mb) for fragmentation profile, or 100000(100kb) for CNV.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            stepNum: Step number for folder name.
            upstream: Not used parameter, do not set this parameter.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """

        super(bamCounter, self).__init__(stepNum, upstream)

        # set fastq input
        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ in ["bamsort", "rmduplicate"]:
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from bamsort.")

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        if chromsizeInput is not None:
            self.setInput("chromsizeInput", chromsizeInput)
        else:
            self.setInput("chromsizeInput", Configure.getConfig("chromSizes"))

        if upstream is None:
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure.getThreads())

        if binlen is not None:
            self.setParam("binlen", binlen)
        else:
            self.setParam("binlen", 100000)

        txtOutput = []
        for x in self.getInput("bamInput"):
            txtOutput.append(
                os.path.join(
                    self.getOutput("outputdir"),
                    self.getMaxFileNamePrefixV2(x).split(".")[0],
                )
                + "_read.txt"
            )
        self.setOutput("txtOutput", txtOutput)

        self.setOutput(
            "bedOutput",
            os.path.join(
                self.getOutput("outputdir"),
                "windows.bed",
            ),
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            multi_run_len = len(self.getInput("bamInput"))
            divide_bin_2(self.getInput("chromsizeInput"), self.getOutput("bedOutput"), self.getParam("binlen"))
            if verbose:
                for i in range(multi_run_len):
                    count_bam(
                        bamInput=self.getInput("bamInput")[i],
                        chromsize=self.getInput("chromsizeInput"),
                        bedOutput=self.getOutput("bedOutput"),
                        txtOutput=self.getOutput("txtOutput")[i],
                        binlen=self.getParam("binlen"),
                    )
            else:
                args = [
                    [
                        self.getInput("bamInput")[i],
                        self.getInput("chromsizeInput"),
                        self.getOutput("bedOutput"),
                        self.getOutput("txtOutput")[i],
                        self.getParam("binlen"),
                    ]
                    for i in range(multi_run_len)
                ]
                self.multiRun(
                    args=args,
                    func=count_bam,
                    nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
                )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
