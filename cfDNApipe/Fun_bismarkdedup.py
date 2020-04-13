# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 11:17:26 2019

@author: Jiaqi Huang

E-mail: huangjq16@mails.tsinghua.edu.cn
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure

__metaclass__ = type


class bismark_deduplicate(StepBase):
    def __init__(
        self,
        bamInput=None,
        outputdir=None,
        threads=1,
        paired=True,
        other_params={},
        stepNum=None,
        upstream=None,
        **kwargs
    ):
        """
        This function is used for removing duplicates from bismark output.
        Note: this function is calling bismark.

        bismark_deduplicate(bamInput=None, outputdir=None, threads=1, paired=True,
                            other_params={}, stepNum=None, upstream=None,)
        {P}arameters:
            bamInput: list, input bam files.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            paired: True for paired data, False for single end data.
            other_params: dict, other parameters passing to Bismark.
                          "-parameter": True means "-parameter" in command line.
                          "-parameter": 1 means "-parameter 1" in command line.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline. This parameter can be True, which means a new pipeline start.
        """

        super(bismark_deduplicate, self).__init__(stepNum, upstream)

        # set bamInput
        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "bismark":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from bismark.")

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

        # set threads and paired
        if upstream is None:
            self.setParam("threads", threads)
            if paired:
                self.setParam("type", "paired")
            else:
                self.setParam("type", "single")
        else:
            self.setParam("threads", Configure.getThreads())
            self.setParam("type", Configure.getType())

        # set other_params
        if self.getParam("type") == "paired":
            other_params.update({"--paired": True})
        elif self.getParam("type") == "single":
            other_params.update({"--single": True})
        else:
            commonError("Wrong data type, must be 'single' or 'paired'!")

        if other_params is None:
            self.setParam("other_params", "")
        else:
            self.setParam("other_params", other_params)

        # set output
        self.setOutput(
            "bamOutput",
            [
                os.path.join(self.getOutput("outputdir"), os.path.splitext(os.path.basename(x))[0],)
                + ".deduplicated.bam"
                for x in self.getInput("bamInput")
            ],
        )
        self.setOutput(
            "reportOutput",
            [
                os.path.join(self.getOutput("outputdir"), os.path.splitext(os.path.basename(x))[0],)
                + ".deduplication_report.txt"
                for x in self.getInput("bamInput")
            ],
        )

        # create cmd
        all_cmd = []
        tmp_cmd = self.cmdCreate(
            [
                "deduplicate_bismark",
                "--output_dir",
                self.getOutput("outputdir"),
                self.getParam("other_params"),
                "--bam",
                self.getInput("bamInput"),
            ]
        )
        all_cmd.append(tmp_cmd)

        # run functions
        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            self.run(all_cmd)

        self.stepInfoRec(cmds=[all_cmd], finishFlag=finishFlag)
