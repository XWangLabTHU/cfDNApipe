# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
Modify on Sun Apr 26 11:27:32 2020
@author: LY, Shuying He
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, maxCore
import os
from .Configure import Configure
import math

__metaclass__ = type


class BQSR(StepBase):
    def __init__(
        self,
        bamInput=None,
        recalInput=None,
        outputdir=None,
        genome=None,
        ref=None,
        stepNum=None,
        upstream=None,
        threads=1,
        verbose=False,
        **kwargs
    ):
        """
        This function is for generating new bam with recaled base.
        Note: this function is calling gatk ApplyBQSR, please install GATK before using.

        BQSR(bamInput=None, recalInput=None, outputdir=None,
        genome=None, ref=None, stepNum=None, upstream=None, threads=1,
        verbose=False, **kwargs)

        {P}arameters:
            bamInput: list, Input bam files.
            recalInput: list, Input recal file from BaseRecalibrator.
            outputdir: str, output result folder, None means the same folder as input files.
            genome: str, genome version, just be hg19 or hg38.
            ref: str, reference genome file path.
            stepNum: int or str, step flag for folder name.
            upstream: upstream output results, used for pipeline, just from BaseRecalibrator. This parameter can be True, which means a new pipeline start.
            threads: int, how many threads used?
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """

        super(BQSR, self).__init__(stepNum, upstream)

        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
            self.setInput("recalInput", recalInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "BaseRecalibrator":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
                self.setInput("recalInput", upstream.getOutput("recalOutput"))
            else:
                raise commonError("Parameter upstream must from rmduplicate.")

        self.checkInputFilePath()

        if upstream is None:
            self.setParam("ref", ref)
            self.setParam("genome", genome)
            self.setParam("threads", threads)

            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)

        else:
            # check Configure for running pipeline
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("ref", Configure.getRefDir())
            self.setParam("genome", Configure.getGenome())
            self.setParam("threads", Configure.getThreads())

        self.BQSRcheck()
        self.setOutput(
            "bamOutput",
            [
                os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + "-BQSR.bam"
                for x in self.getInput("bamInput")
            ],
        )

        multi_run_len = len(self.getInput("bamInput"))

        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(
                [
                    "gatk",
                    "ApplyBQSR",
                    "-I",
                    self.getInput("bamInput")[i],
                    "-R",
                    self.getParam("ref") + "/" + self.getParam("genome") + ".fa",
                    "--bqsr-recal-file",
                    self.getInput("recalInput")[i],
                    "-O",
                    self.getOutput("bamOutput")[i],
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
                    nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
                )

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)

    def BQSRcheck(
        self,
    ):
        fafile = os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa")

        if not os.path.exists(fafile):
            raise commonError("file " + fafile + " don not exist!")
