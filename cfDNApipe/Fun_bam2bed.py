# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 08:35:29 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError, bamTobed, bamTobedForSingle
import os
from .Configure import Configure


__metaclass__ = type


class bam2bed(StepBase):
    def __init__(self, bamInput=None, outputdir=None, paired=True, stepNum=None, upstream=None, **kwargs):
        """
        This function is used for converting bam file to bed file.

        bam2bed(bamInput=None, outputdir=None, paired=True, stepNum=None, upstream=None)
        {P}arameters:
            bamInput: list, input bam files.
            outputdir: str, output result folder, None means the same folder as input files.
            paired: boolean, paired end or single end.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline.
        """

        super(bam2bed, self).__init__(stepNum, upstream)

        # set bamInput
        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "bamsort":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            elif upstream.__class__.__name__ == "rmduplicate":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from bamsort or rmduplicate.")

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set paired
        if upstream is None:
            if paired:
                self.setParam("type", "paired")
            else:
                self.setParam("type", "single")
        else:
            self.setParam("type", Configure.getType())

        self.setOutput(
            "bedOutput",
            [
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + ".bed"
                for x in self.getInput("bamInput")
            ],
        )
        self.setOutput(
            "bedgzOutput",
            [
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + ".bed.gz"
                for x in self.getInput("bamInput")
            ],
        )
        self.setOutput(
            "tbiOutput",
            [
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + ".bed.gz.tbi"
                for x in self.getInput("bamInput")
            ],
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            multi_run_len = len(self.getInput("bamInput"))

            if self.getParam("type") == "paired":
                for i in range(multi_run_len):
                    print("Now, converting file: " + self.getInput("bamInput")[i])
                    bamTobed(
                        bamInput=self.getInput("bamInput")[i], bedOutput=self.getOutput("bedOutput")[i],
                    )
            elif self.getParam("type") == "single":
                for i in range(multi_run_len):
                    print("Now, converting file: " + self.getInput("bamInput")[i])
                    bamTobedForSingle(
                        bamInput=self.getInput("bamInput")[i], bedOutput=self.getOutput("bedOutput")[i],
                    )
            else:
                commonError("Wrong data type, must be 'single' or 'paired'!")

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
