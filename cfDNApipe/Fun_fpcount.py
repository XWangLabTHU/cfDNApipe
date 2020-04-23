# -*- coding: utf-8 -*-
"""
Created on Wed Apr 8 12:51:24 2020

@author: Jiaqi Huang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, count_fragprof
import os
from .Configure import Configure

__metaclass__ = type


class fpCounter(StepBase):
    def __init__(
        self,
        bedgzInput=None,
        chromsizeInput=None,
        blacklistInput=None,
        gapInput=None,
        domains=None,
        outputdir=None,
        stepNum=None,
        upstream=None,
        **kwargs
    ):
        """
        This function is used for counting short and long reads of the input bedgz files.

        fpCounter(bedgzInput=None, chromsizeInput=None, blacklistInput=None, gapInput=None, domains=None, outputdir=None, stepNum=None, upstream=None)
        {P}arameters:
            bedgzInput: list, paths of input bedgz files waiting to be processed.
            chromsizeInput: str, path of chromsize file.
            blacklistInput: str, path of blacklist file.
            gapInput: str, path of gap file.
            domains: list, [minimum_length_of_short_fragments, maximum_length_of_short_fragments, minimum_length_of_long_fragments, maximum_length_of_long_fragments]; default is [100, 150, 151, 220].
            outputdir: str, output result folder, None means the same folder as input files.
            stepNum: Step number for folder name.
            upstream: Not used parameter, do not set this parameter.
        """

        super(fpCounter, self).__init__(stepNum, upstream)

        # set fastq input
        if (upstream is None) or (upstream is True):
            self.setInput("bedgzInput", bedgzInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "bam2bed":
                self.setInput("bedgzInput", upstream.getOutput("bedgzOutput"))
            else:
                raise commonError("Parameter upstream must from bam2bed.")

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput("outputdir", os.path.dirname(os.path.abspath(self.getInput("bedgzInput")[0])))
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        if chromsizeInput is not None:
            self.setInput("chromsizeInput", chromsizeInput)
        else:
            self.setInput("chromsizeInput", Configure.getConfig("chromSizes"))

        if blacklistInput is not None:
            self.setInput("blacklistInput", blacklistInput)
        else:
            self.setInput("blacklistInput", Configure.getConfig("Blacklist"))

        if gapInput is not None:
            self.setInput("gapInput", gapInput)
        else:
            self.setInput("gapInput", Configure.getConfig("Gaps"))

        if domains is not None:
            self.setParam("domain", domains)
        else:
            self.setParam("domain", [100, 150, 151, 220])

        txtOutput = []
        for x in self.getInput("bedgzInput"):
            txtOutput.append(
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x).split(".")[0]) + "_short.txt"
            )
            txtOutput.append(
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x).split(".")[0]) + "_long.txt"
            )

        self.setOutput("txtOutput", txtOutput)

        self.setOutput("bedOutput", os.path.join(self.getOutput("outputdir"), "windows.bed",))

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            count_fragprof(
                bedgzInput=self.getInput("bedgzInput"),
                chromsize=self.getInput("chromsizeInput"),
                blacklist=self.getInput("blacklistInput"),
                gap=self.getInput("gapInput"),
                bedOutput=self.getOutput("bedOutput"),
                txtOutput=self.getOutput("txtOutput"),
                domain=self.getParam("domain"),
                binlen=5000000,
            )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
