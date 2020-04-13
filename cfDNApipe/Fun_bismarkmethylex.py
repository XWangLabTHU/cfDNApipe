# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 18:35:24 2019

@author: Jiaqi Huang

E-mail: huangjq16@mails.tsinghua.edu.cn
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
import math
from .Configure import Configure

__metaclass__ = type


class bismark_methylation_extractor(StepBase):
    def __init__(
        self,
        bamInput=None,
        outputdir=None,
        threads=1,
        other_params={
            "--no_overlap": True,
            "--report": True,
            "--no_header": True,
            "--gzip": True,
            "--bedGraph": True,
            "--zero_based": True,
        },
        paired=True,
        stepNum=None,
        upstream=None,
        **kwargs
    ):
        """
        This function is used for extracting methylation information from bismark output.
        Note: this function is calling bismark.

        bismark_methylation_extractor(bamInput=None, outputdir=None, threads=1,
            other_params={
                "--no_overlap": True,
                "--report": True,
                "--no_header": True,
                "--gzip": True,
                "--bedGraph": True,
                "--zero_based": True,
            }, paired=True, stepNum=None, upstream=None,)
        {P}arameters:
            bamInput: list, input bam files.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            paired: True for paired data, False for single end data.
            other_params: dict, other parameters passing to FASTQC.
                          "-parameter": True means "-parameter" in command line.
                          "-parameter": 1 means "-parameter 1" in command line.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline.
        """

        super(bismark_methylation_extractor, self).__init__(stepNum, upstream)

        # set fastq input
        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "bismark" or "bismark_deduplicate":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from bismark or bismark_deduplicate.")

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

        # set threads, paired,
        if upstream is None:
            self.setParam("threads", threads)
            if paired:
                self.setParam("type", "paired")
            else:
                self.setParam("type", "single")

        else:
            self.setParam("type", Configure.getType())
            self.setParam("threads", Configure.getThreads())

        # update other_params
        if self.getParam("type") == "paired":
            other_params.update({"--paired-end": True})
        elif self.getParam("type") == "single":
            other_params.update({"--single-end": True})
        else:
            commonError("Wrong data type, must be 'single' or 'paired'!")

        other_params.update({"--multicore": math.ceil(self.getParam("threads") / 3)})
        other_params.update({"--output": self.getOutput("outputdir")})

        # set other_params
        if other_params is None:
            self.setParam("other_params", "")
        else:
            self.setParam("other_params", other_params)

        self.setOutput(
            "covOutput",
            [
                os.path.join(self.getOutput("outputdir"), os.path.splitext(os.path.basename(x))[0],)
                + ".bedGraph.gz.bismark.zero.cov"
                for x in self.getInput("bamInput")
            ],
        )
        self.setOutput(
            "bedGraphOutput",
            [
                os.path.join(self.getOutput("outputdir"), os.path.splitext(os.path.basename(x))[0],) + ".bedGraph.gz"
                for x in self.getInput("bamInput")
            ],
        )
        self.setOutput(
            "covgzOutput",
            [
                os.path.join(self.getOutput("outputdir"), os.path.splitext(os.path.basename(x))[0],) + ".bismark.cov.gz"
                for x in self.getInput("bamInput")
            ],
        )
        self.setOutput(
            "reportOutput",
            [
                os.path.join(self.getOutput("outputdir"), os.path.splitext(os.path.basename(x))[0],)
                + "_splitting_report.txt"
                for x in self.getInput("bamInput")
            ],
        )

        all_cmd = []

        tmp_cmd = self.cmdCreate(
            ["bismark_methylation_extractor", self.getParam("other_params"), self.getInput("bamInput"),]
        )
        all_cmd.append(tmp_cmd)

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            self.run(all_cmd)

        self.stepInfoRec(cmds=[all_cmd], finishFlag=finishFlag)
