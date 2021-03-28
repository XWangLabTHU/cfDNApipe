# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 17:27:13 2021

@author: Hanwen Xu, Jiaqi Huang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, predecipher, decipher, deconplot
from .Configure import Configure
import os
import pandas as pd
import collections
import numpy as np

__metaclass__ = type


class deconvolution(StepBase):
    def __init__(
        self,
        mixInput=None,
        refInput=None,
        outputdir=None,
        threads=1,
        stepNum=None,
        upstream=None,
        marker_path="",
        scale=0.1,
        delcol_factor=10,
        iter_num=10,
        confidence=0.75,
        w_thresh=10,
        unknown=False,
        is_markers=False,
        **kwargs
    ):
        """
        This function is used for methylation signal deconvolution.

        deconvolution(mixInput=None, refInput=None, outputdir=None, threads=1, stepNum=None, upstream=None, marker_path='', scale=0.1, delcol_factor=10, iter_num=10, confidence=0.75, w_thresh=10, unknown=False, is_markers=False, is_methylation=True)
        {P}arameters:
            mixInput: Input samples need to be deconvoluted.
            refInput: reference files. Default from https://www.pnas.org/content/112/40/E5503.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use. In this function, this number is set to 1.
            upstream: upstream output results, used for pipeline, must from calculate_methyl.
            stepNum: int or str, step flag for folder name.
            marker_path: str, path to markers, if users select to specify certain markers
            scale: float, control the convergence of SVR
            delcol_factor: int, control the extent of removing collinearity
            iter_num: int, iterative numbers of outliers detection
            confidence: float, ratio of remained markers in each outlier detection loop
            w_thresh: int, threshold to cut the weights designer
            unknown: bool, if there is unknown content
            is_markers: bool, if users choose to specify their own markers
        """
        super(deconvolution, self).__init__(stepNum, upstream)

        # set input
        if (upstream is None) or (upstream is True):
            self.setInput("mixInput", mixInput)
            if refInput is None:
                self.setInput("refInput", Configure.getConfig("PlasmaMarker"))
            else:
                self.setInput("refInput", refInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "calculate_methyl":
                self.setInput("mixInput", upstream.getOutput("txtOutput"))
            else:
                raise commonError("Parameter upstream must from calculate_methyl.")

            if refInput is None:
                self.setInput("refInput", Configure.getConfig("PlasmaMarker"))
            else:
                self.setInput("refInput", refInput)

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(self.getInput("mixInput"))),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set threads
        if upstream is None:
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure.getThreads())

        self.setParam("marker_path", marker_path)
        self.setParam("scale", scale)
        self.setParam("delcol_factor", delcol_factor)
        self.setParam("iter_num", iter_num)
        self.setParam("confidence", confidence)
        self.setParam("w_thresh", w_thresh)
        self.setParam("unknown", unknown)
        self.setParam("is_markers", is_markers)
        self.setParam("is_methylation", True)

        self.setOutput("txtOutput", os.path.join(self.getOutput("outputdir"), "result.txt"))
        self.setOutput("plotOutput", os.path.join(self.getOutput("outputdir"), "bar_chart.png"))

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            mix, ref = predecipher(self.getInput("mixInput"), self.getInput("refInput"))
            res_df = decipher(
                ref=ref,
                mix=mix,
                save_path=self.getOutput("txtOutput"),
                marker_path=self.getParam("marker_path"),
                scale=self.getParam("scale"),
                delcol_factor=self.getParam("delcol_factor"),
                iter_num=self.getParam("iter_num"),
                confidence=self.getParam("confidence"),
                w_thresh=self.getParam("w_thresh"),
                unknown=self.getParam("unknown"),
                is_markers=self.getParam("is_markers"),
                is_methylation=self.getParam("is_methylation"),
            )
            deconplot(res_df, self.getOutput("plotOutput"))

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
