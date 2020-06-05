# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 11:45:21 2020

@author: Shuying He

E-mail: heshuying@fzidt.com
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import math

__metaclass__ = type


class cnvPlot(StepBase):
    # draw diagram and scatter plot for one
    def __init__(
        self,
        cnsInput=None,
        cnrInput=None,
        outputdir=None,
        diagram=True,
        diagram_params={"--threshold": 0.5, "--min-probes": 3, "-y": True, },
        scatter=True,
        scatter_params={"--y-max": 2, "--y-min": -2, "--segment-color": "'red'", },
        threads=1,
        verbose=False,
        stepNum=None,
        upstream=None,
        **kwargs
    ):

        super(cnvPlot, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):

            self.setInput("cnsInput", cnsInput)
            self.setInput("cnrInput", cnrInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "cnvbatch":
                self.setInput("cnsInput", upstream.getOutput("cnsOutput"))
                self.setInput("cnrInput", upstream.getOutput("cnrOutput"))
            else:
                raise commonError("Parameter upstream must from cnvbatch.")
        self.checkInputFilePath()

        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("cnsInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
            self.setParam("threads", threads)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())
        # cns and cnr checking...
        figure_number = len(self.getInput("cnsInput"))
        prefix = []
        for i in range(figure_number):
            cns = self.getInput("cnsInput")[i]
            cns_prefix = self.getMaxFileNamePrefixV2(cns)

            cnr = self.getInput("cnrInput")[i]
            cnr_prefix = self.getMaxFileNamePrefixV2(cnr)
            if cns_prefix == cnr_prefix:
                prefix.append(cns_prefix)
            else:
                raise commonError("Error: File %s and %s is not match." % (cns, cnr))
        self.setParam("prefix", prefix)

        # cmd create
        all_cmd = []
        if diagram:
            self.setParam("diagram_params", diagram_params)
            self.setOutput(
                "diagram_pdf",
                [
                    os.path.join(self.getOutput("outputdir"), x + "_diagram.pdf")
                    for x in self.getParam("prefix")
                ],
            )
            for i in range(figure_number):
                cmd = self.cmdCreate(
                    [
                        "cnvkit.py",
                        "diagram",
                        self.getInput("cnrInput")[i],
                        "-s",
                        self.getInput("cnsInput")[i],
                        "-o",
                        self.getOutput("diagram_pdf")[i],
                        "--title",
                        self.getParam("prefix")[i],
                        self.getParam("diagram_params"),
                    ]
                )
                all_cmd.append(cmd)

        if scatter:
            self.setParam("scatter_params", scatter_params)
            self.setOutput(
                "scatter_pdf",
                [
                    os.path.join(self.getOutput("outputdir"), x + "_scatter.pdf")
                    for x in self.getParam("prefix")
                ],
            )
            for i in range(figure_number):
                cmd = self.cmdCreate(
                    [
                        "cnvkit.py",
                        "scatter",
                        self.getInput("cnrInput")[i],
                        "-s",
                        self.getInput("cnsInput")[i],
                        "-o",
                        self.getOutput("scatter_pdf")[i],
                        "--title",
                        self.getParam("prefix")[i],
                        self.getParam("scatter_params"),
                    ]
                )
                all_cmd.append(cmd)

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            if verbose:
                self.run(all_cmd)
            else:
                self.multiRun(
                    args=all_cmd,
                    func=None,
                    nCore=math.ceil(self.getParam("threads") / 4),
                )

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)
