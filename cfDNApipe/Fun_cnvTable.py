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


class cnvTable(StepBase):
    # caculate break, genemetrics and gene for one
    def __init__(
        self,
        cnsInput=None,
        cnrInput=None,
        outputdir=None,
        breaks=True,
        breaks_params={"--min-probes": 1, },
        genemetrics=True,
        genemetrics_params={"--threshold": 0.1, "--min-probes": 3, "-y": True, },
        threads=1,
        verbose=False,
        stepNum=None,
        upstream=None,
        **kwargs,
    ):

        super(cnvTable, self).__init__(stepNum, upstream)

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
                    os.path.dirname(os.path.abspath(self.getInput("cnrInput")[0])),
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
        if breaks:
            self.setParam("breaks_params", breaks_params)
            self.setOutput(
                "breaks_txt",
                [
                    os.path.join(self.getOutput("outputdir"), x + "_breaks.txt")
                    for x in self.getParam("prefix")
                ],
            )
            for i in range(figure_number):
                cmd = self.cmdCreate(
                    [
                        "cnvkit.py",
                        "breaks",
                        self.getInput("cnrInput")[i],
                        self.getInput("cnsInput")[i],
                        self.getParam("breaks_params"),
                        "-o",
                        self.getOutput("breaks_txt")[i],
                    ]
                )
                all_cmd.append(cmd)

        if genemetrics:
            self.setParam("genemetrics_params", genemetrics_params)
            self.setOutput(
                "genemetrics_cnrs",
                [
                    os.path.join(
                        self.getOutput("outputdir"), x + "_genemetrics_cnrs.txt"
                    )
                    for x in self.getParam("prefix")
                ],
            )
            self.setOutput(
                "genemetrics_cnr",
                [
                    os.path.join(
                        self.getOutput("outputdir"), x + "_genemetrics_cnr.txt"
                    )
                    for x in self.getParam("prefix")
                ],
            )
            self.setOutput(
                "cnrs_gene",
                [
                    os.path.join(self.getOutput("outputdir"), x + "_cnrs_gene.txt")
                    for x in self.getParam("prefix")
                ],
            )
            self.setOutput(
                "cnr_gene",
                [
                    os.path.join(self.getOutput("outputdir"), x + "_cnr_gene.txt")
                    for x in self.getParam("prefix")
                ],
            )
            self.setOutput(
                "genemetrics_gene",
                [
                    os.path.join(
                        self.getOutput("outputdir"), x + "_genemetrics_gene.txt"
                    )
                    for x in self.getParam("prefix")
                ],
            )

            for i in range(figure_number):
                cmd = self.cmdCreate(
                    [
                        "cnvkit.py",
                        "genemetrics",
                        self.getInput("cnrInput")[i],
                        "-s",
                        self.getInput("cnsInput")[i],
                        self.getParam("genemetrics_params"),
                        "-o",
                        self.getOutput("genemetrics_cnrs")[i],
                        ";",
                        "cnvkit.py",
                        "genemetrics",
                        self.getInput("cnrInput")[i],
                        self.getParam("genemetrics_params"),
                        "-o",
                        self.getOutput("genemetrics_cnr")[i],
                        ";",
                        "tail -n +2",
                        self.getOutput("genemetrics_cnrs")[i],
                        "| cut -f 1 | sort >",
                        self.getOutput("cnrs_gene")[i],
                        ";",
                        "tail -n +2",
                        self.getOutput("genemetrics_cnr")[i],
                        "| cut -f 1 | sort >",
                        self.getOutput("cnr_gene")[i],
                        ";",
                        "comm",
                        "-12",
                        self.getOutput("cnrs_gene")[i],
                        self.getOutput("cnr_gene")[i],
                        ">",
                        self.getOutput("genemetrics_gene")[i],
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
