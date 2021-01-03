# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 17:03:17 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError, maxCore
import os
from .Configure import Configure
import math


__metaclass__ = type


class rmduplicate(StepBase):
    def __init__(
        self,
        bamInput=None,
        outputdir=None,
        Xmx="4G",
        threads=1,
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for removing duplicates in WGS data.
        Note: this function is calling picard.

        rmduplicate(bamInput=None, outputdir=None, threads=1, stepNum=None, upstream=None, verbose=True)
        {P}arameters:
            bamInput: list, bam file input.
            outputdir: str, output result folder, None means the same folder as input files.
            Xmx: How many memory will be used for every thread, default: 4G.
            threads: int, how many thread to use.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """

        super(rmduplicate, self).__init__(stepNum, upstream)

        # set fastq input
        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "bamsort":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from bamsort.")

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            self.setParam("Xmx", Xmx)
            self.setParam("threads", threads)
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[1])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setParam("Xmx", Configure.getJavaMem())
            self.setParam("threads", Configure.getThreads())
            self.setOutput("outputdir", self.getStepFolderPath())

        self.setOutput(
            "bamOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "-rmdup.bam"
                for x in self.getInput("bamInput")
            ],
        )
        self.setOutput("baiOutput", [x + ".bai" for x in self.getOutput("bamOutput")])
        self.setOutput(
            "metricsOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "-rmdup.txt"
                for x in self.getInput("bamInput")
            ],
        )

        multi_run_len = len(self.getInput("bamInput"))

        cmd_step1 = []
        cmd_step2 = []

        for i in range(multi_run_len):
            cmd1 = self.cmdCreate(
                [
                    "gatk",
                    "--java-options",
                    '"-Xmx%s"' % self.getParam("Xmx"),
                    "MarkDuplicates",
                    "--REMOVE_DUPLICATES",
                    "true",
                    "--INPUT",
                    self.getInput("bamInput")[i],
                    "--METRICS_FILE",
                    self.getOutput("metricsOutput")[i],
                    "--OUTPUT",
                    self.getOutput("bamOutput")[i],
                ]
            )

            cmd2 = self.cmdCreate(
                [
                    "samtools index",
                    "-@",
                    self.getParam("threads"),
                    self.getOutput("bamOutput")[i],
                ]
            )

            cmd_step1.append(cmd1)
            cmd_step2.append(cmd2)

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            if verbose:
                self.run(cmd_step1)
                self.run(cmd_step2)
            else:
                self.multiRun(
                    args=cmd_step1,
                    func=None,
                    nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
                )
                self.multiRun(
                    args=cmd_step2,
                    func=None,
                    nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
                )

        self.stepInfoRec(cmds=[cmd_step1, cmd_step2], finishFlag=finishFlag)
