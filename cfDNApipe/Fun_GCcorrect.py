# -*- coding: utf-8 -*-
"""
Created on Wed Apr 8 14:25:33 2020

@author: Jiaqi Huang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, wig2df, correctReadCount
import pandas as pd
import numpy as np
import os
from .Configure import Configure

__metaclass__ = type


class GCCorrect(StepBase):
    def __init__(
            self,
            readwigInput=None, # list
            gcwigInput=None, # list
            outputdir=None,  # str
            stepNum=None,
            readupstream=None,
            gcupstream=None,
            **kwargs
    ):
        """
        This function is used for processing GC correction on read count data in wig files.

        GCcorrect(readwigInput=None, gcwigInput=None, outputdir=None, stepNum=None, readupstream=None, gcupstream=None)
        {P}arameters:
            readwigInput: list, paths of wig files of read counts.
            gcwigInput: list, paths of wig files of gc contents.
            outputdir: str, output result folder, None means the same folder as input files.
            stepNum: Step number for folder name.
            readupstream: Not used parameter, do not set this parameter.
            gcupstream: Not used parameter, do not set this parameter.
        """
            
        super(GCCorrect, self).__init__(stepNum, readupstream)

        if readupstream is None or gcupstream is None:
            self.setInput("readwigInput", readwigInput)
            self.setInput("gcwigInput", gcwigInput)
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput("outputdir", os.path.dirname(
                    os.path.abspath(self.getInput("readwigInput")[0])))
            else:
                self.setOutput("outputdir", outputdir)
        
        else:
            Configure.configureCheck()
            readupstream.checkFilePath()
            gcupstream.checkFilePath()

            if readupstream.__class__.__name__ == "runCounter":
                self.setInput("readwigInput", readupstream.getOutput("wigOutput"))
            else:
                raise commonError("Parameter upstream must from runCounter.")
            if gcupstream.__class__.__name__ == "runCounter":
                self.setInput("gcwigInput", gcupstream.getOutput("wigOutput"))
            else:
                raise commonError("Parameter upstream must from runCounter.")
            self.checkInputFilePath()

            self.setOutput("outputdir", self.getStepFolderPath())

        self.setOutput("txtOutput", [os.path.join(self.getOutput(
            "outputdir"), self.getMaxFileNamePrefixV2(x)) + ".txt" for x in self.getInput("readwigInput")])
            
        self.setOutput("plotOutput", [os.path.join(self.getOutput(
            "outputdir"), self.getMaxFileNamePrefixV2(x)) + ".png" for x in self.getInput("readwigInput")])

        finishFlag = self.stepInit(readupstream)
        
        multi_run_len = len(self.getInput("readwigInput"))
        
        if not finishFlag:
            gc_df = wig2df(self.getInput("gcwigInput")[0])
            for i in range(multi_run_len):
                print("Now, processing", self.getMaxFileNamePrefixV2(
                    self.getInput("readwigInput")[i]), "...")
                read_df = wig2df(self.getInput("readwigInput")[i])
                correctReadCount(
                    read_df, gc_df, self.getOutput("txtOutput")[i], self.getOutput("plotOutput")[i])

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
