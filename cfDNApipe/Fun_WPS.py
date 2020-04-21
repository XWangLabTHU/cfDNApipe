# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 19:52:09 2020

@author: Jiaqi Huang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, processWPS
import os
from .Configure import Configure
import math

__metaclass__ = type


class runWPS(StepBase):
    def __init__(
        self,
        bedgzInput=None,
        tsvInput=None,
        outputdir=None,
        protect=None,
        empty=None,
        insertsize=None,
        threads=None,
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for running WPS.

        runWPS(bedgzInput=None, tsvInput=None, outputdir=None, protect=None, empty=None, insertsize=None, threads=None, stepNum=None, upstream=None, verbose=True)
        {P}arameters:
            bedgzInput: list, input bed.gz files.
            tsvInput: str, regions of transcript file.
            outputdir: str, output result folder, None means the same folder as input files.
            protect: int, base pair protection assumed for elements (default 120).
            empty: bool, keep files of empty blocks (default False).
            insertsize: list, minimum and maximum read length threshold to consider, shaped like [mininslen, maxinslen] (default None).
            threads: int, how many thread to use.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline. This parameter can be True, which means a new pipeline start.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """
        super(runWPS, self).__init__(stepNum, upstream)

        # set bedgz input
        if (upstream is None) or (upstream is True):
            self.setInput("bedgzInput", bedgzInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "bam2bed":
                self.setInput("bedgzInput", upstream.getOutput("bedzOutput"))
            else:
                raise commonError("Parameter upstream must from bam2bed.")

        self.checkInputFilePath()
        
        # set tsv input
        if tsvInput is not None:
            self.setInput("tsvInput", tsvInput)
        else:
            self.setInput("tsvInput", Configure.getConfig("dummy")) # need to be checked!
            
        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(self.getInput("bedgzInput")[1])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set threads
        if threads is not None:
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure.getThreads())
            
        # set protect
        if protect is not None:
            self.setParam("protect", protect)
        else:
            self.setParam("protect", 120)
            
        # set empty
        if empty is not None:
            self.setParam("empty", empty)
        else:
            self.setParam("empty", False)
            
        # set insertsize
        if insertsize is not None:
            self.setParam("insertsize", insertsize)
        else:
            self.setParam("insertsize", [-1, -1])
        
        dirs = []
        for x in self.getInput("bedgzInput"):
            newdir = os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x))
            if not os.path.exists(newdir):
                os.mkdir(newdir)
            dirs.append(newdir)
        
        self.setOutput(
            "sampleOutputdir",
            dirs,
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            multi_run_len = len(self.getInput("bedgzInput"))
            if verbose:
                for i in range(multi_run_len):
                    processWPS(
                        bedgzInput=self.getInput("bedgzInput")[i],
                        tsvInput=self.getInput("tsvInput"),
                        protectInput=self.getParam("protect"),
                        outputfile=os.path.join(
                            self.getOutput("sampleOutputdir")[i], 
                            self.getMaxFileNamePrefixV2(self.getInput("bedgzInput")[i])
                        ) + "_%s.tsv.gz",
                        empty=self.getParam("empty"),
                        minInsSize=self.getParam("insertsize")[0],
                        maxInsSize=self.getParam("insertsize")[1],
                    )
            else:
                args = [
                    [
                        self.getInput("bedgzInput")[i],
                        self.getInput("tsvInput"),
                        self.getParam("protect"),
                        os.path.join(
                            self.getOutput("sampleOutputdir")[i], 
                            self.getMaxFileNamePrefixV2(self.getInput("bedgzInput")[i])
                        ) + "_%s.tsv.gz",
                        self.getParam("empty"),
                        self.getParam("insertsize")[0],
                        self.getParam("insertsize")[1],
                    ]
                    for i in range(multi_run_len)
                ]
                self.multiRun(args=args, func=processWPS, nCore=math.ceil(self.getParam("threads") / 4))

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
