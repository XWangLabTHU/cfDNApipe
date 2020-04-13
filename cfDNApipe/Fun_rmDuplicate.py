# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 17:03:17 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure


__metaclass__ = type


class rmduplicate(StepBase):
    def __init__(self, bamInput=None, outputdir=None, threads=1, stepNum=None, upstream=None, **kwargs):
        """
        This function is used for removing duplicates in WGS data.
        Note: this function is calling picard.

        rmduplicate(bamInput=None, outputdir=None, threads=1, stepNum=None, upstream=None)
        {P}arameters:
            bamInput: list, bam file input.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline.
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
                raise commonError("Parameter upstream must from inputprocess or adapterremoval.")

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(self.getInput("bamInput")[1])),
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
            "bamOutput",
            [
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + "-rmdup.bam"
                for x in self.getInput("bamInput")
            ],
        )
        self.setOutput("baiOutput", [x + ".bai" for x in self.getOutput("bamOutput")])
        self.setOutput(
            "metricsOutput",
            [
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + "-rmdup.txt"
                for x in self.getInput("bamInput")
            ],
        )

        multi_run_len = len(self.getInput("bamInput"))

        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(
                [
                    "picard MarkDuplicates",
                    "REMOVE_DUPLICATES=true",
                    "INPUT=" + self.getInput("bamInput")[i],
                    "OUTPUT=" + self.getOutput("bamOutput")[i],
                    "METRICS_FILE=" + self.getOutput("metricsOutput")[i],
                    ";",
                    "samtools index",
                    "-@",
                    self.getParam("threads"),
                    self.getOutput("bamOutput")[i],
                ]
            )
            all_cmd.append(tmp_cmd)

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            self.run(all_cmd)

        self.stepInfoRec(cmds=[all_cmd], finishFlag=finishFlag)
