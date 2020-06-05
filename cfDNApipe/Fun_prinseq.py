# -*- coding: utf-8 -*-
"""
Created on Tue Mar 3 18:27:32 2020
@author: He Shuying
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure


__metaclass__ = type


class lowcomplexityfilter(StepBase):
    def __init__(
        self,
        seq1Input=None,
        seq2Input=None,
        threads=1,
        type="paired",
        other_params={"-y": True, "-Y": 30},
        OutputDir=None,
        stepNum=None,
        upstream=None,
        **kwargs
    ):
        super(lowcomplexityfilter, self).__init__(stepNum, upstream)

        if (upstream is None) or (upstream is True):
            self.setParam("type", type)
            self.setInput("seq1Input", seq1Input)
            if self.getParam("type") == "paired":
                self.setInput("seq2Input", seq2Input)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            self.setParam("type", Configure.getType())

            if upstream.__class__.__name__ in ["bismark", "bowtie2"]:
                if self.getParam("type") == "paired":
                    self.setInput("seq1Input", upstream.getOutput("unmapped-1"))
                    self.setInput("seq2Input", upstream.getOutput("unmapped-2"))
                elif self.getParam("type") == "single":
                    self.setInput("seq1Input", upstream.getOutput("unmapped"))
                else:
                    raise commonError("Analysis date type should be paired or single.")
            else:
                raise commonError("Parameter upstream must from bismark.")

        self.checkInputFilePath()

        # because the warning:
        # fastp uses up to 16 threads although you specified 25
        if upstream is None:
            self.setParam("threads", min(threads, 16))
            if OutputDir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("seq1Input")[0])),
                )
            else:
                self.setOutput("outputdir", OutputDir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", min(Configure.getThreads(), 16))

        self.setParam(
            "prefix",
            [self.getMaxFileNamePrefixV2(x) for x in self.getInput("seq1Input")],
        )

        self.setOutput(
            "jsonOutput",
            [
                os.path.join(self.getOutput("outputdir"), x + ".json")
                for x in self.getParam("prefix")
            ],
        )

        self.setOutput(
            "htmlOutput",
            [
                os.path.join(self.getOutput("outputdir"), x + ".html")
                for x in self.getParam("prefix")
            ],
        )

        self.setOutput(
            "seq1Output",
            [
                os.path.join(self.getOutput("outputdir"), x + ".lc_filter_R1.fq.gz")
                for x in self.getParam("prefix")
            ],
        )

        if self.getParam("type") == "paired":
            self.setOutput(
                "seq2Output",
                [
                    os.path.join(self.getOutput("outputdir"), x + ".lc_filter_R2.fq.gz")
                    for x in self.getParam("prefix")
                ],
            )

        if other_params is None:
            self.setParam("other_params", "")
        else:
            self.setParam("other_params", other_params)

        multi_run_len = len(self.getInput("seq1Input"))

        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = [
                "fastp",
                "--thread",
                self.getParam("threads"),
                self.getParam("other_params"),
                "--json",
                self.getOutput("jsonOutput")[i],
                "--html",
                self.getOutput("htmlOutput")[i],
                "--report_title '" + self.getParam("prefix")[i] + "'",
                "-i",
                self.getInput("seq1Input")[i],
                "-o",
                self.getOutput("seq1Output")[i],
            ]
            if self.getParam("type") == "paired":
                tmp_cmd.extend(
                    [
                        "-I",
                        self.getInput("seq2Input")[i],
                        "-O",
                        self.getOutput("seq2Output")[i],
                    ]
                )

            all_cmd.append(self.cmdCreate(tmp_cmd))

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            self.run(all_cmd)

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)
