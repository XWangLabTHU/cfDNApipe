# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:26:39 2019

@author: zhang
"""

from .Configure import Configure
from .StepBase import StepBase
from .cfDNA_utils import commonError
import os

__metaclass__ = type


class inputprocess(StepBase):
    def __init__(
        self, fqInput1=None, fqInput2=None, inputFolder=None, stepNum=1, upstream=None,
    ):
        """
        This function is used for auto input files arrangement.
        Note: this function is designed for pipeline, so user must set Configure before using.

        inputprocess(fqInput1=None, fqInput2=None, inputFolder=None, stepNum=1, upstream=None)
        {P}arameters:
            fqInput1: list, fastq files for single end data;  _1 files for paired end data.
            fqInput2: list, [] for single end data;  _2 files for paired end data.
            inputFolder: input folder contains all your input fastq data, the program will detect inputs automatically.
            stepNum: Step number for folder name.
            upstream: Not used parameter, do not set this parameter.
        """
        super(inputprocess, self).__init__(stepNum, upstream)

        # check Configure for running pipeline
        Configure.configureCheck()

        self.setParam("type", Configure.getType())

        # using folder first, ignore "fqInput1" and "fqInput2"
        if inputFolder is not None:
            all_files = os.listdir(inputFolder)
            if len(all_files) == 0:
                raise commonError("The input folder is empty.")
            all_files.sort()
            all_files = list(map(lambda x: os.path.join(inputFolder, x), all_files))
            if self.getParam("type") == "paired":
                print("paired data is processing......")
                fqInput1 = []
                fqInput2 = []
                for i in range(len(all_files)):
                    if i % 2:
                        fqInput2.append(all_files[i])
                    else:
                        fqInput1.append(all_files[i])
                self.setInput("fq1", fqInput1)
                self.setInput("fq2", fqInput2)
            elif self.getParam("type") == "single":
                print("single data is processing......")
                fqInput1 = all_files
                self.setInput("fq1", fqInput1)
                self.setInput("fq2", [])
            else:
                commonError("Wrong data tpye, must be 'single' or 'paired'!")
        else:
            if self.getParam("type") == "paired":
                print("paired data is processing......")
                self.setInput("fq1", fqInput1)
                self.setInput("fq2", fqInput2)
            elif self.getParam("type") == "single":
                print("single data is processing......")
                self.setInput("fq1", fqInput1)
                self.setInput("fq2", [])
            else:
                commonError("Wrong data tpye, must be 'single' or 'paired'!")

        self.checkInputFilePath()

        self.setOutput("fq1", self.getInput("fq1"))
        self.setOutput("fq2", self.getInput("fq2"))

        self.setOutput("outputdir", self.getStepFolderPath())

        finishFlag = self.stepInit(upstream=True)

        if not finishFlag:
            pass

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
