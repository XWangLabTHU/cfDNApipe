# -*- coding: utf-8 -*-
"""
Created on Wed Apr 8 14:25:33 2020

@author: Jiaqi Huang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, correctReadCount, maxCore
import os
from .Configure import Configure
import math

__metaclass__ = type


class GCCorrect(StepBase):
    def __init__(
        self,
        readInput=None,
        gcwigInput=None,
        readtype=None,
        corrkey=None,
        outputdir=None,
        threads=1,
        stepNum=None,
        readupstream=None,
        gcupstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for processing GC correction on read count data in wig or csv/txt files.

        GCcorrect(readInput=None, gcwigInput=None, readtype=None, corrkey=None, outputdir=None, threads=1, stepNum=None, readupstream=None, gcupstream=None, verbose=True,)
        {P}arameters:
            readInput: list, paths of input files of read counts.
            gcwigInput: list, paths of wig files of gc contents.
            readtype: int, file type of readInput, 1 for .wig, 2 for .txt/.csv.; 1 is set by default.
            corrkey: char, type of GC correction, "-" for minus, "/" for divide, "0" for process without GC correction; "/" is set by default
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            stepNum: Step number for folder name.
            readupstream: upstream output results, used for pipeline.
            gcupstream: upstream output results, used for pipeline.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """

        super(GCCorrect, self).__init__(stepNum, readupstream)

        # set some input
        if (
            ((readupstream is None) and (gcupstream is None))
            or (readupstream is True)
            or (gcupstream is True)
        ):
            self.setInput("readInput", readInput)
            self.setInput("gcwigInput", gcwigInput)
            if readtype is not None:
                self.setParam("readtype", readtype)
            else:
                self.setParam("readtype", 1)
            if corrkey is not None:
                self.setParam("corrkey", corrkey)
            else:
                self.setParam("corrkey", "/")
        else:
            Configure.configureCheck()
            readupstream.checkFilePath()
            gcupstream.checkFilePath()
            if readupstream.__class__.__name__ == "runCounter":
                self.setInput("readInput", readupstream.getOutput("wigOutput"))
                self.setParam("readtype", 1)
                if corrkey is not None:
                    self.setParam("corrkey", corrkey)
                else:
                    self.setParam("corrkey", "/")
            elif readupstream.__class__.__name__ == "fpCounter":
                self.setInput("readInput", readupstream.getOutput("txtOutput"))
                self.setParam("readtype", 2)
                if corrkey is not None:
                    self.setParam("corrkey", corrkey)
                else:
                    self.setParam("corrkey", "-")
            elif readupstream.__class__.__name__ == "bamCounter":
                self.setInput("readInput", readupstream.getOutput("txtOutput"))
                self.setParam("readtype", 2)
                if corrkey is not None:
                    self.setParam("corrkey", corrkey)
                else:
                    self.setParam("corrkey", "/")
            else:
                raise commonError(
                    "Parameter upstream must from runCounter, fpCounter or bamCounter."
                )
            if gcupstream.__class__.__name__ == "runCounter":
                self.setInput("gcwigInput", gcupstream.getOutput("wigOutput"))
            else:
                raise commonError("Parameter upstream must from runCounter.")

        self.checkInputFilePath()

        # set threads
        if (readupstream is None) and (gcupstream is None):
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure.getThreads())

        if (readupstream is None) and (gcupstream is None):
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("readInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.checkInputFilePath()
            self.setOutput("outputdir", self.getStepFolderPath())

        self.setOutput(
            "txtOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "_gc_cor.txt"
                for x in self.getInput("readInput")
            ],
        )

        self.setOutput(
            "plotOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "_gc_cor.png"
                for x in self.getInput("readInput")
            ],
        )

        finishFlag = self.stepInit(readupstream)

        multi_run_len = len(self.getInput("readInput"))

        if not finishFlag:
            if verbose:
                for i in range(multi_run_len):
                    print(
                        "Now, processing",
                        self.getMaxFileNamePrefixV2(self.getInput("readInput")[i]),
                        "...",
                    )
                    correctReadCount(
                        self.getInput("readInput")[i],
                        self.getInput("gcwigInput")[0],
                        self.getOutput("txtOutput")[i],
                        self.getOutput("plotOutput")[i],
                        self.getParam("corrkey"),
                        self.getParam("readtype"),
                    )
            else:
                args = [
                    [
                        self.getInput("readInput")[i],
                        self.getInput("gcwigInput")[0],
                        self.getOutput("txtOutput")[i],
                        self.getOutput("plotOutput")[i],
                        self.getParam("corrkey"),
                        self.getParam("readtype"),
                    ]
                    for i in range(multi_run_len)
                ]
                self.multiRun(
                    args=args,
                    func=correctReadCount,
                    nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
                )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
