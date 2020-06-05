# -*- coding: utf-8 -*-
"""
Created on Tue Mar 3 18:27:32 2020
@author: He Shuying
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import math


__metaclass__ = type


class samtofastq(StepBase):
    def __init__(
        self,
        bamInput=None,
        OutputDir=None,
        stepNum=None,
        upstream=None,
        threads=1,
        verbose=False,
        **kwargs
    ):
        super(samtofastq, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ in [
                "rmduplicate",
                "bismark",
                "bismark_deduplicate",
                "bowtie2",
                "mapQfilter",
            ]:
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from mapQfilter.")
        self.checkInputFilePath()

        if upstream is True:
            self.setParam("threads", threads)
            if OutputDir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", OutputDir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())

        self.setParam(
            "prefix",
            [self.getMaxFileNamePrefixV2(x) for x in self.getInput("bamInput")],
        )

        self.setOutput(
            "fq1Output",
            [
                os.path.join(self.getOutput("outputdir"), x + ".R1.fastq")
                for x in self.getParam("prefix")
            ],
        )

        self.setOutput(
            "fq2Output",
            [
                os.path.join(self.getOutput("outputdir"), x + ".R2.fastq")
                for x in self.getParam("prefix")
            ],
        )

        multi_run_len = len(self.getInput("bamInput"))

        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(
                [
                    "picard",
                    "SamToFastq",
                    "I=" + self.getInput("bamInput")[i],
                    "FASTQ=" + self.getOutput("fq1Output")[i],
                    "SECOND_END_FASTQ=" + self.getOutput("fq2Output")[i],
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
