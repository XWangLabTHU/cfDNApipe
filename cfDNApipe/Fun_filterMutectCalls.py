# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
Modify on Sun Apr 26 11:27:32 2020
@author: LY, Shuying He
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import math

__metaclass__ = type


class filterMutectCalls(StepBase):
    def __init__(
        self,
        vcfInput=None,
        contaminationInput=None,
        outputdir=None,
        genome=None,
        ref=None,
        other_params=None,
        stepNum=None,
        upstream=None,
        threads=1,
        verbose=False,
        **kwargs,
    ):

        super(filterMutectCalls, self).__init__(stepNum, upstream)
        chromosome = ["chr%i" % x for x in range(1, 23)]
        chromosome.extend(
            ["chrX", "chrY", "chrM", ]
        )

        if upstream is None:
            self.setParam("ref", ref)
            self.setParam("genome", genome)
            self.setParam("threads", threads)
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("vcfInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)

        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("ref", Configure.getRefDir())
            self.setParam("genome", Configure.getGenome())
            self.setParam("threads", Configure.getThreads())

        if (upstream is None) or (upstream is True):
            self.setInput("vcfInput", vcfInput)
            if contaminationInput is not None:
                self.setInput(
                    "contaminationInput", self.convertToList(contaminationInput)
                )

            self.setOutput(
                "vcfOutput",
                [
                    os.path.join(
                        self.getOutput("outputdir"),
                        self.getMaxFileNamePrefixV2(x) + ".filtered.vcf.gz",
                    )
                    for x in self.getInput("vcfInput")
                ],
            )
            self.setOutput(
                "summaryOutput",
                [
                    os.path.join(
                        self.getOutput("outputdir"),
                        self.getMaxFileNamePrefixV2(x)
                        + ".unfiltered.vcf.gz.filteringStats.tsv",
                    )
                    for x in self.getInput("vcfInput")
                ],
            )
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "mutect2t" or "mutect2n":
                self.setInput("indir", upstream.getOutput("outputdir"))
                self.setParam(
                    "prefix",
                    [os.path.basename(x) for x in upstream.getOutput("outdir")],
                )
                if "contaminationOutput" in upstream.getOutputs():
                    self.setInput(
                        "contaminationInput", upstream.getOutput("contaminationOutput")
                    )
                elif contaminationInput is not None:
                    self.setInput(
                        "contaminationInput", self.convertToList(contaminationInput)
                    )
                else:
                    pass
            else:
                raise commonError("Parameter upstream must from mutect2t or mutect2n.")

            vcfInput = []
            vcfOutput = []
            summaryOutput = []

            for x in self.getParam("prefix"):
                if not os.path.exists(os.path.join(self.getOutput("outputdir"), x)):
                    os.makedirs(os.path.join(self.getOutput("outputdir"), x))

                vcfInput.extend(
                    [
                        os.path.join(
                            self.getInput("indir"), f"{x}/{x}_{y}.unfiltered.vcf.gz"
                        )
                        for y in chromosome
                    ]
                )
                vcfOutput.extend(
                    [
                        os.path.join(
                            self.getOutput("outputdir"), f"{x}/{x}_{y}.filtered.vcf.gz"
                        )
                        for y in chromosome
                    ]
                )
                summaryOutput.extend(
                    [
                        os.path.join(
                            self.getOutput("outputdir"),
                            f"{x}/{x}_{y}.filtered.vcf.gz.filteringStats.tsv",
                        )
                        for y in chromosome
                    ]
                )

            self.setInput("vcfInput", vcfInput)
            self.setOutput("vcfOutput", vcfOutput)
            self.setOutput("summaryOutput", summaryOutput)

        self.checkInputFilePath()
        self.filterMutectCallscheck()

        if other_params is None:
            self.setParam("other_params", "")
        else:
            self.setParam("other_params", other_params)

        # cmd create
        vcfnum = len(self.getInput("vcfInput"))
        all_cmd = []
        if "contaminationInput" not in self.getInputs():
            for i in range(vcfnum):
                tmp_cmd = self.cmdCreate(
                    [
                        "gatk",
                        "FilterMutectCalls",
                        "-V",
                        self.getInput("vcfInput")[i],
                        "-R",
                        self.getParam("ref") + "/" + self.getParam("genome") + ".fa",
                        self.getParam("other_params"),
                        "-O",
                        self.getOutput("vcfOutput")[i],
                    ]
                )
                all_cmd.append(tmp_cmd)

        elif len(self.getInput("contaminationInput")) == 1:
            for i in range(vcfnum):
                tmp_cmd = self.cmdCreate(
                    [
                        "gatk",
                        "FilterMutectCalls",
                        "-V",
                        self.getInput("vcfInput")[i],
                        "-R",
                        self.getParam("ref") + "/" + self.getParam("genome") + ".fa",
                        "--contamination-table",
                        self.getInput("contaminationInput")[0],
                        self.getParam("other_params"),
                        "-O",
                        self.getOutput("vcfOutput")[i],
                    ]
                )
                all_cmd.append(tmp_cmd)

        elif (
            len(self.getInput("contaminationInput")) > 1
            and len(self.getInput("vcfInput")) % 25 == 0
        ):
            for i in range(vcfnum):
                tmp_cmd = self.cmdCreate(
                    [
                        "gatk",
                        "FilterMutectCalls",
                        "-V",
                        self.getInput("vcfInput")[i],
                        "-R",
                        self.getParam("ref") + "/" + self.getParam("genome") + ".fa",
                        "--contamination-table",
                        self.getInput("contaminationInput")[int(i / 25)],
                        "-L",
                        chromosome[int(i % 25)],
                        self.getParam("other_params"),
                        "-O",
                        self.getOutput("vcfOutput")[i],
                    ]
                )
                all_cmd.append(tmp_cmd)
        else:
            raise commonError("Unmatch number for contaminationInput and vcfInput.")

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

    def filterMutectCallscheck(self,):
        fafile = os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa")

        if not os.path.exists(fafile):
            raise commonError("file " + fafile + " don not exist!")
