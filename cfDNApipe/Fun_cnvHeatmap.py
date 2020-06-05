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
    # create a heatmap including all samples
    def __init__(
        self,
        cnrInput=None,
        outputdir=None,
        orderfile=None,  # a file contain the cnr file name to reorder the outfile
        other_params={"-d": True, "-y": True, },
        # other_params setting
        # '-c chr8' designate chromosome
        # '-g EXT1,PXDNL'
        stepNum=None,
        upstream=None,
        **kwargs,
    ):

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

        self.setOutput(
            "heatmap", os.path.join(self.getOutput("outputdir"), "heatmap.pdf")
        )
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

    def checkorderfile(self,):
        """for checking the file in orderfile."""
        with open(self.getInput("orderfile")) as INPUT:
            for line in INPUT:
                tmp = line.strips()
                if not os.path.exists(tmp):
                    raise commonError("%s in orderfile is not exists." % tmp)
