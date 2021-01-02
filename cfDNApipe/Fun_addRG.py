# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
Modify on Tues Feb 18 16:27:32 2020
@author: LY, He Shuying
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import math

__metaclass__ = type


class addRG(StepBase):
    def __init__(
        self, bamInput=None, outputdir=None, Xmx="4G", upstream=None, stepNum=None, threads=1, verbose=False, **kwargs
    ):
        """
        This function is used for adding read group info for BAM file.
        Note: this function is calling gatk's AddOrReplaceReadGroups function.

        addRG(bamInput=None, outputdir=None, Xmx="4G", 
              stepNum=None, threads=1, upstream=None, verbose=False)
              
        {P}arameters:
            bamInput: list, input bam files.
            outputdir: str, output result folder, None means the same folder as input files.
            Xmx: str, Xmx mem, default is "4g".
            threads: int, how many thread to use, using multiRun to excute.
            stepNum: int, step number for folder name.
            upstream: upstream output results(rmduplicate / deduplicate_bismark / bamsort / bismark), used for pipeline. This parameter can be True, which means a new pipeline start.
        """

        super(addRG, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
            # In this situation, input file and output path should be checked
            self.setInput("bamInput", bamInput)
        else:
            # check Configure for running pipeline
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "rmduplicate" or "deduplicate_bismark":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from rmduplicate or deduplicate_bismark.")

        self.checkInputFilePath()

        if upstream is None:
            self.setParam("Xmx", Xmx)
            self.setParam("threads", threads)
            if outputdir is None:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setParam("threads", Configure.getThreads())
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("Xmx", Configure.getJavaMem())

        self.setOutput(
            "bamOutput",
            [
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + "-RG.bam"
                for x in self.getInput("bamInput")
            ],
        )
        self.setParam("RGID", [self.getMaxFileNamePrefixV2(x) for x in self.getInput("bamInput")])
        self.setParam("RGLB", [self.getMaxFileNamePrefixV2(x) for x in self.getInput("bamInput")])
        self.setParam("RGPL", "ILLUMINA")
        self.setParam("RGSM", [self.getMaxFileNamePrefixV2(x) for x in self.getInput("bamInput")])
        self.setParam("RGPU", "Null")

        multi_run_len = len(self.getInput("bamInput"))

        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(
                [
                    "gatk",
                    "--java-options",
                    '"-Xmx%s"' % self.getParam("Xmx"),
                    "AddOrReplaceReadGroups",
                    "--INPUT",
                    self.getInput("bamInput")[i],
                    "--OUTPUT",
                    self.getOutput("bamOutput")[i],
                    "--SORT_ORDER coordinate",
                    "--RGID",
                    self.getParam("RGID")[i],
                    "--RGLB",
                    self.getParam("RGLB")[i],
                    "--RGPL",
                    self.getParam("RGPL"),
                    "--RGSM",
                    self.getParam("RGSM")[i],
                    "--RGPU",
                    self.getParam("RGPU"),
                    "--CREATE_INDEX True",
                ]
            )
            all_cmd.append(tmp_cmd)

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            if verbose:
                self.run(all_cmd)
            else:
                # this step is forced to run multiple thread less than 8
                nCore = math.ceil(self.getParam("threads") / 4)
                if nCore > 8:
                    nCore = 8
                    print("The thread number is forced to 8!")

                self.multiRun(
                    args=all_cmd, func=None, nCore=nCore,
                )

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)
