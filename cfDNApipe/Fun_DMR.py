# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:51:54 2019

@author: zhang
"""

from .StepBase2 import StepBase2
from .cfDNA_utils import commonError, computeCUE
import pandas as pd
import os
import math
from .Configure2 import Configure2

__metaclass__ = type


class DMR(StepBase2):
    def __init__(
        self,
        casebedInput=None,
        ctrlbedInput=None,
        refRegInput=None,
        outputdir=None,
        threads=1,
        caseupstream=None,
        ctrlupstream=None,
        stepNum=None, 
        verbose=True,
        **kwargs
    ):
        """
        This function is used for compute OCF values of input bed files.

        computeOCF(casebedInput=None, ctrlbedInput=None, refRegInput=None, outputdir=None, threads=1, caseupstream=None, ctrlupstream=None, labelInput=None, stepNum=None,  verbose=True)
        {P}arameters:
            casebedInput: list, input bed files of case samples.
            ctrlbedInput: list, input bed files of control samples.
            refRegInput: str, reference file for OCF computing.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            caseupstream: upstream output results, used for pipeline.
            ctrlupstream: upstream output results, used for pipeline.
            stepNum: int, step number for folder name.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """