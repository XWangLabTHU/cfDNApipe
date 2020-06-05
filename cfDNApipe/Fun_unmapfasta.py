# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
@author: LY, Shuying He
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import math
import pkg_resources

__metaclass__ = type


class unmapfasta(StepBase):
    def __init__(
        self,
        plInput=None,
        bamInput=None,
        fq1Input=None,
        fq2Input=None,
        OutputDir=None,
        stepNum=None,
        upstream=None,
        threads=1,
        verbose=False,
        **kwargs
    ):

        super(unmapfasta, self).__init__(stepNum, upstream)

        if (upstream is None) or (upstream is True):
            # In this situation, input file and output path should be checked
            self.setInput("bamInput", bamInput)
            if fq1Input:
                self.setInput("unmapped-1", fq1Input)
            if fq2Input:
                self.setInput("unmapped-2", fq2Input)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "bowtie2":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
                self.setInput("unmapped-1", upstream.getOutput("unmapped-1"))
                self.setInput("unmapped-2", upstream.getOutput("unmapped-2"))
            else:
                raise commonError("Parameter upstream must from bowtie2.")

        self.checkInputFilePath()

        if upstream is None:
            if OutputDir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", OutputDir)
            self.setParam("threads", threads)

        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())

        if plInput:
            self.setInput("plInput", plInput)
        else:
            self.setInput(
                "plInput",
                pkg_resources.resource_filename(
                    "cfDNApipe", "data/VirusFinder2.0Plus/preprocessPlus_V2.pl"
                ),
            )

        self.setParam(
            "outdir",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                for x in self.getInput("bamInput")
            ],
        )

        self.setOutput(
            "unmapped-1", [x + "/unmapped.1.fa" for x in self.getParam("outdir")]
        )
        self.setOutput(
            "unmapped-2", [x + "/unmapped.2.fa" for x in self.getParam("outdir")]
        )

        multi_run_len = len(self.getInput("bamInput"))
        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = [
                "perl",
                self.getInput("plInput"),
                "-bam",
                self.getInput("bamInput")[i],
                "-output",
                self.getParam("outdir")[i],
            ]

            if "unmapped-1" in self.getInputs():
                tmp_cmd.append("-fq1 " + self.getInput("unmapped-1")[i])
            if "unmapped-2" in self.getInputs():
                tmp_cmd.append("-fq2 " + self.getInput("unmapped-2")[i])

            tmp_cmd = self.cmdCreate(tmp_cmd)
            all_cmd.append(tmp_cmd)

        finishFlag = self.stepInit(upstream)

        if finishFlag is False:
            for i in self.getParam("outdir"):
                if not os.path.exists(i):
                    os.mkdir(i)
            if verbose:
                self.run(all_cmd)
            else:
                self.multiRun(
                    args=all_cmd,
                    func=None,
                    nCore=math.ceil(self.getParam("threads") / 4),
                )

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)
