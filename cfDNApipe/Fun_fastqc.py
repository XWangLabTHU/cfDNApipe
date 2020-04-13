# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 18:27:32 2019

@author: zhang
"""

from .StepBase import StepBase
from .cfDNA_utils import flatten
import os
from .Configure import Configure

__metaclass__ = type


class fastqc(StepBase):
    def __init__(
        self,
        fastqInput=None,
        fastqcOutputDir=None,
        threads=1,
        other_params=None,
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for fastq file quality control.
        Note: this function is calling FASTQC.

        fastqc(fastqInput=None, fastqcOutputDir=None, threads=1, other_params=None, stepNum=None, upstream=None, verbose=True)
        {P}arameters:
            fastqInput: list, fastq files.
            fastqcOutputDir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            other_params: dict, other parameters passing to FASTQC.
                          "-parameter": True means "-parameter" in command line.
                          "-parameter": 1 means "-parameter 1" in command line.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """
        super(fastqc, self).__init__(stepNum, upstream)

        # set fastqInput
        if (upstream is None) or (upstream is True):
            self.setInput("fastqInputs", fastqInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            self.setInput(
                "fastqInputs", list(flatten([upstream.getOutput("fq1"), upstream.getOutput("fq2")])),
            )

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if fastqcOutputDir is None:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(self.getInput("fastqInputs")[0])),
                )
            else:
                self.setOutput("outputdir", fastqcOutputDir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set threads
        if upstream is None:
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure.getThreads())

        # set other_params
        if other_params is None:
            self.setParam("other_params", "")
        else:
            self.setParam("other_params", other_params)

        # create cmd
        cmd = self.cmdCreate(
            [
                "fastqc",
                "--outdir",
                self.getOutput("outputdir"),
                "--threads",
                self.getParam("threads"),
                self.getParam("other_params"),
                self.inputs["fastqInputs"],
            ]
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            if verbose:
                self.run(cmd)
            else:
                self.multiRun(args=[cmd], func=None, nCore=1)

        self.stepInfoRec(cmds=[cmd], finishFlag=finishFlag)
