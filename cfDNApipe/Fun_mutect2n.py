# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
@author: LY
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import math

__metaclass__ = type


class mutect2n(StepBase):
    def __init__(
        self,
        bamInput=None,
        outputdir=None,
        other_params={"-max-mnp-distance": 0},
        genome=None,
        ref=None,  # str
        stepNum=None,
        upstream=None,
        threads=1,
        verbose=False,
        **kwargs
    ):

        super(mutect2n, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
            # In this situation, input file and output path should be checked
            self.setInput("bamInput", bamInput)
        else:
            # check Configure for running pipeline
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "addRG":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            elif upstream.__class__.__name__ == "BQSR":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from addRG or BQSR.")
        self.checkInputFilePath()

        if upstream is None:
            self.setParam("ref", ref)
            self.setParam("genome", genome)
            self.setParam("threads", threads)

            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("ref", Configure.getRefDir())
            self.setParam("genome", Configure.getGenome())
            self.setParam("threads", Configure.getThreads())

        self.mutect2ncheck()
        self.setOutput(
            "vcfOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + ".unfiltered.vcf.gz"
                for x in self.getInput("bamInput")
            ],
        )
        self.setOutput("tbiOutput", [x + ".tbi" for x in self.getOutput("vcfOutput")])
        self.setOutput("bamOutput", self.getInput("bamInput"))

        self.setParam(
            "prefix",
            [self.getMaxFileNamePrefixV2(x) for x in self.getInput("bamInput")],
        )
        self.setOutput(
            "outdir",
            [
                os.path.join(self.getOutput("outputdir"), z)
                for z in self.getParam("prefix")
            ],
        )

        chromosome = ["chr%i" % x for x in range(1, 23)]
        chromosome.extend(["chrX", "chrY", "chrM"])

        for i in range(len(self.getParam("prefix"))):
            self.setParam(
                "%s_vcfInput" % self.getParam("prefix")[i],
                [
                    os.path.join(
                        self.getOutput("outdir")[i],
                        self.getParam("prefix")[i] + "_" + y + ".unfiltered.vcf.gz",
                    )
                    for y in chromosome
                ],
            )

        if other_params is None:
            self.setParam("other_params", "")
        else:
            self.setParam("other_params", other_params)

        # cmd create
        multi_run_len = len(self.getInput("bamInput"))
        all_cmd = []

        for i in range(multi_run_len):
            for y in range(len(chromosome)):
                tmp_cmd = self.cmdCreate(
                    [
                        "gatk",
                        "Mutect2",
                        "-I",
                        self.getInput("bamInput")[i],
                        "-R",
                        self.getParam("ref") + "/" + self.getParam("genome") + ".fa",
                        self.getParam("other_params"),
                        "-L",
                        chromosome[y],
                        "-O",
                        self.getOutput("outdir")[i]
                        + "/"
                        + self.getParam("prefix")[i]
                        + "_"
                        + chromosome[y]
                        + ".unfiltered.vcf.gz",
                    ]
                )
                all_cmd.append(tmp_cmd)

        gather_cmd = []
        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(
                [
                    "gatk",
                    "GatherVcfs",
                    "-I",
                    " -I ".join(
                        self.getParam("%s_vcfInput" % self.getParam("prefix")[i])
                    ),
                    "-O",
                    self.getOutput("vcfOutput")[i],
                    ";" "gatk",
                    "IndexFeatureFile",
                    "-I",
                    self.getOutput("vcfOutput")[i],
                ]
            )
            gather_cmd.append(tmp_cmd)

        # excute
        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            for i in self.getOutput("outdir"):
                if not os.path.exists(i):
                    os.mkdir(i)
            if verbose:
                self.run(all_cmd)
                self.run(gather_cmd)
            else:
                self.multiRun(
                    args=all_cmd,
                    func=None,
                    nCore=math.ceil(self.getParam("threads") / 4),
                )
                self.multiRun(
                    args=gather_cmd,
                    func=None,
                    nCore=math.ceil(self.getParam("threads") / 4),
                )

        self.stepInfoRec(cmds=all_cmd + gather_cmd, finishFlag=finishFlag)

    def mutect2ncheck(self,):
        """check ref exists."""
        fafile = os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa")

        if not os.path.exists(fafile):
            raise commonError("file " + fafile + " don not exist!")
