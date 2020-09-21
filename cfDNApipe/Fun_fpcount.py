# -*- coding: utf-8 -*-
"""
Created on Wed Apr 8 12:51:24 2020

@author: Jiaqi Huang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, count_fragprof, divide_bin_1, divide_bin_2
import os
import math
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
        binlen=None,
        outputdir=None,
        threads=1,
        processtype=None,
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for counting total reads or short-and-long-reads of the input bedgz files.

        fpCounter(bedgzInput=None, chromsizeInput=None, blacklistInput=None, gapInput=None, domains=None,
                  binlen=None, outputdir=None, threads=1, processtype=None, stepNum=None, upstream=None, verbose=True,)
        {P}arameters:
            bedgzInput: list, paths of input bedgz files waiting to be processed.
            chromsizeInput: str, path of chromsize file.
            blacklistInput: str, used in fragmentation profile, path of blacklist file.
            gapInput: str, used in fragmentation profile, path of gap file.
            domains: list, used in fragmentation profile, [minimum_length_of_short_fragments, maximum_length_of_short_fragments,
                     minimum_length_of_long_fragments, maximum_length_of_long_fragments]; default is [100, 150, 151, 220].
            binlen: int, length of each bin; default is 5000000(5Mb) for fragmentation profile, or 100000(100kb) for CNV.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            processtype: int, 1 for fragmentation profile, 2 for CNV.
            stepNum: Step number for folder name.
            upstream: Not used parameter, do not set this parameter.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
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
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bedgzInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        if processtype == 1 or processtype == 2 :
            self.setParam("processtype", processtype)
        else:
            raise commonError("Please refer processtype with 1 or 2.")

        if chromsizeInput is not None:
            self.setInput("chromsizeInput", chromsizeInput)
        else:
            self.setInput("chromsizeInput", Configure.getConfig("chromSizes"))

        if self.getParam("processtype") == 1:
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

        if upstream is None:
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure.getThreads())

        if binlen is not None:
            self.setParam("binlen", binlen)
        elif self.getParam("processtype") == 1:
            self.setParam("binlen", 5000000)
        elif self.getParam("processtype") == 2:
            self.setParam("binlen", 100000)

        txtOutput = []
        if self.getParam("processtype") == 1:
            for x in self.getInput("bedgzInput"):
                txtOutput.append(
                    os.path.join(
                        self.getOutput("outputdir"),
                        self.getMaxFileNamePrefixV2(x).split(".")[0],
                    )
                    + "_short.txt"
                )
                txtOutput.append(
                    os.path.join(
                        self.getOutput("outputdir"),
                        self.getMaxFileNamePrefixV2(x).split(".")[0],
                    )
                    + "_long.txt"
                )
        elif self.getParam("processtype") == 2:
            for x in self.getInput("bedgzInput"):
                txtOutput.append(
                    os.path.join(
                        self.getOutput("outputdir"),
                        self.getMaxFileNamePrefixV2(x).split(".")[0],
                    )
                    + "_read.txt"
                )
        self.setOutput("txtOutput", txtOutput)

        self.setOutput(
            "bedOutput", os.path.join(self.getOutput("outputdir"), "windows.bed",)
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            multi_run_len = len(self.getInput("bedgzInput"))
            if not os.path.exists(self.getOutput("bedOutput")):
                if self.getParam("processtype") == 1:
                    divide_bin_1(
                        chromsize=self.getInput("chromsizeInput"),
                        blacklist=self.getInput("blacklistInput"),
                        gap=self.getInput("gapInput"),
                        windows=self.getOutput("bedOutput"),
                        binlen=self.getParam("binlen")
                    )
                elif self.getParam("processtype") == 2:
                    divide_bin_2(
                        chromsize=self.getInput("chromsizeInput"),
                        windows=self.getOutput("bedOutput"),
                        binlen=self.getParam("binlen")
                    )
            if verbose:
                if self.getParam("processtype") == 1:
                    for i in range(multi_run_len):
                        count_fragprof(
                            bedgzInput=self.getInput("bedgzInput")[i],
                            bedOutput=self.getOutput("bedOutput"),
                            txtOutput=self.getOutput("txtOutput")[2 * i : 2 * i + 2],
                            domain=self.getParam("domain"),
                            binlen=self.getParam("binlen"),
                            type=1,
                        )
                elif self.getParam("processtype") == 2:
                    for i in range(multi_run_len):
                        count_fragprof(
                            bedgzInput=self.getInput("bedgzInput")[i],
                            bedOutput=self.getOutput("bedOutput"),
                            txtOutput=self.getOutput("txtOutput")[i],
                            binlen=self.getParam("binlen"),
                            type=2,
                        )
            else:
                if self.getParam("processtype") == 1:
                    args = [
                        [
                            self.getInput("bedgzInput")[i],
                            self.getOutput("bedOutput"),
                            self.getOutput("txtOutput")[2 * i : 2 * i + 2],
                            self.getParam("domain"),
                            self.getParam("binlen"),
                            1,
                        ]
                        for i in range(multi_run_len)
                    ]
                elif self.getParam("processtype") == 2:
                    args = [
                        [
                            self.getInput("bedgzInput")[i],
                            self.getOutput("bedOutput"),
                            self.getOutput("txtOutput")[i],
                            None,
                            self.getParam("binlen"),
                            2,
                        ]
                        for i in range(multi_run_len)
                    ]
                self.multiRun(
                    args=args,
                    func=count_fragprof,
                    nCore=math.ceil(self.getParam("threads") / 4),
                )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
