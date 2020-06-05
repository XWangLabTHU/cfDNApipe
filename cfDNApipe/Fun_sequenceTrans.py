# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 18:27:32 2019

@author: Jiaqi Huang
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError
from .Configure import Configure
import os
import math


__metaclass__ = type


class sequencetransfer(StepBase):
    def __init__(
        self,
        bamInput=None,
        bedInput=None,
        outputdir=None,
        threads=1,
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for transferring sequences.

        bam2bed(bamInput=None, bedInput=None, outputdir=None, threads=1, stepNum=None, upstream=None, verbose=True)
        {P}arameters:
            bamInput: list, input bam files.
            bedInput: list, input bed files.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """

        super(sequencetransfer, self).__init__(stepNum, upstream)

        # set bamInput
        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "rmduplicate":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError(
                    "Parameter upstream must from inputprocess or adapterremoval."
                )

        self.checkInputFilePath()

        bedflag = False

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set bedInput, threads
        if upstream is None:
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure.getThreads())

        if bedInput is not None:
            bedflag = True
            self.setInput("bedInput", bedInput)
        else:
            commonError("Parameter bedInput is None!")

        self.setOutput(
            "txtOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "-seq.txt"
                for x in self.getInput("bamInput")
            ],
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            multi_run_len = len(self.getInput("bamInput"))
            if verbose:
                for i in range(multi_run_len):
                    if bedflag:
                        self.seqTrans(
                            bamInput=self.getInput("bamInput")[i],
                            txtOutput=self.getOutput("txtOutput")[i],
                            bedInput=self.getInput("bedInput"),
                        )
                    else:
                        self.seqTrans(
                            bamInput=self.getInput("bamInput")[i],
                            txtOutput=self.getOutput("txtOutput")[i],
                        )
            else:
                if bedflag:
                    args = [
                        [
                            self.getInput("bamInput")[i],
                            self.getOutput("txtOutput")[i],
                            self.getInput("bedInput"),
                        ]
                        for i in range(multi_run_len)
                    ]
                else:
                    args = [
                        [self.getInput("bamInput")[i], self.getOutput("txtOutput")[i], ]
                        for i in range(multi_run_len)
                    ]
                self.multiRun(
                    args=args,
                    func=self.seqTrans,
                    nCore=math.ceil(self.getParam("threads") / 4),
                )
        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
