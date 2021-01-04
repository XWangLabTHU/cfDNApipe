# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 11:45:21 2020

@author: Shuying He

E-mail: heshuying@fzidt.com
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure


__metaclass__ = type


class cnvHeatmap(StepBase):
    def __init__(
        self,
        cnrInput=None,
        outputdir=None,
        orderfile=None,
        other_params={"-d": True, "-y": True},
        stepNum=None,
        upstream=None,
        **kwargs,
    ):
        """
        This function is used for drawing heatmap plot including all samples.
        Note: This function is calling cnvkit.py diagram / cnvkit.py heatmap, please install cnvkit before using.

        cnvHeatmap(cnrInput=None,
            outputdir=None, orderfile=None,
            other_params={"-d": True, "-y": True},
            stepNum=None, upstream=None, **kwargs)

        {P}arameters:
            cnrInput: list, cnr files( a table of copy number ratios), generating from cnvkit.py batch.
            outputdir: str, output result folder, None means the same folder as input files.
            orderfile: str, a file contain the cnr file name to reorder the outfile, you can miss this parameter which will ordered as the cnrInput.
            diagram: True, drawing diagram plot? default is True.
            other_params: dict, parameter for cnvkit.py heatmap, default is {{"-d": True, "-y": True}. other parameter could setting as this {"-c": "chr8"} as for designate chromosome. {"-g": "EXT1,PXDNL"} is for designate gene analysis.
            stepNum: int or str, step flag for folder name.
            upstream: upstream output results, used for cnv pipeline, just can be cnvbatch. This parameter can be True, which means a new pipeline start.
        """

        super(cnvHeatmap, self).__init__(stepNum, upstream)

        if orderfile:
            self.setInput("orderfile", orderfile)
            self.checkInputFilePath()
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("orderfile"))),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            if (upstream is None) or (upstream is True):
                self.setInput("cnrInput", cnrInput)
                self.checkInputFilePath()
            else:
                Configure.configureCheck()
                upstream.checkFilePath()

                if upstream.__class__.__name__ == "cnvbatch":
                    self.setInput("cnrInput", upstream.getOutput("cnrOutput"))
                else:
                    raise commonError("Parameter upstream must from cnvbatch.")

            if upstream is None:
                if outputdir is None:
                    self.setOutput(
                        "outputdir",
                        os.path.dirname(os.path.abspath(self.getInput("cnrInput")[0])),
                    )
                else:
                    self.setOutput("outputdir", outputdir)

            else:
                self.setOutput("outputdir", self.getStepFolderPath())

        self.setOutput("heatmap", os.path.join(self.getOutput("outputdir"), "heatmap.pdf"))
        self.setParam("other_params", other_params)

        # create cmd
        if "orderfile" in self.getInputs():
            cmd = self.cmdCreate(
                [
                    "cnvkit.py",
                    "heatmap",
                    "`cat",
                    self.getInput("orderfile"),
                    "`",
                    self.getParam("other_params"),
                    "-o",
                    self.getOutput("heatmap"),
                ]
            )
        else:
            cmd = self.cmdCreate(
                [
                    "cnvkit.py",
                    "heatmap",
                    " ".join(self.getInput("cnrInput")),
                    self.getParam("other_params"),
                    "-o",
                    self.getOutput("heatmap"),
                ]
            )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            self.run(cmd)

        self.stepInfoRec(cmds=[cmd], finishFlag=finishFlag)

    def checkorderfile(
        self,
    ):
        """for checking the file in orderfile."""
        with open(self.getInput("orderfile")) as INPUT:
            for line in INPUT:
                tmp = line.strips()
                if not os.path.exists(tmp):
                    raise commonError("%s in orderfile is not exists." % tmp)
