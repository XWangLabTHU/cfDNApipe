# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
@author: Shuying He
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure


__metaclass__ = type


class mapQ_filter(StepBase):
    def __init__(
        self,
        bamInput=None,
        mapQ=20,
        threads=1,
        OutputDir=None,
        stepNum=None,
        upstream=None,
        **kwargs,
    ):
        super(mapQ_filter, self).__init__(stepNum, upstream)
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
                "sort",
            ]:
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from bismark_deduplicate.")

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

        self.setParam("mapQ", mapQ)

        self.setOutput(
            "bamOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + f".mq{mapQ}.bam"
                for x in self.getInput("bamInput")
            ],
        )

        self.setOutput("baiOutput", [x + ".bai" for x in self.getOutput("bamOutput")])

        multi_run_len = len(self.getInput("bamInput"))
        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(
                [
                    "samtools",
                    "view",
                    "-b -h -q",
                    self.getParam("mapQ"),
                    self.getInput("bamInput")[i],
                    "-@",
                    self.getParam("threads"),
                    "|",
                    "samtools",
                    "sort",
                    "-@",
                    self.getParam("threads"),
                    "-",
                    "-o",
                    self.getOutput("bamOutput")[i],
                    ";",
                    "samtools",
                    "index",
                    "-@",
                    self.getParam("threads"),
                    self.getOutput("bamOutput")[i],
                ]
            )
            all_cmd.append(tmp_cmd)

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            self.run(all_cmd)

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)
