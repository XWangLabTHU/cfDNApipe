# -*- coding: utf-8 -*-
"""
Created on Thur Feb 27 18:27:32 2020
@author: He Shuying
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import math
import pkg_resources

__metaclass__ = type


class VCFpostprocess(StepBase):
    def __init__(
        self,
        BisSNP=None,  # Bissnp.jar
        vcfInput=None,
        memSize="2g",
        OutputDir=None,
        genome=None,
        ref=None,  # str
        stepNum=None,
        upstream=None,
        threads=1,
        verbose=False,
        **kwargs,
    ):

        super(VCFpostprocess, self).__init__(stepNum, upstream)

        if upstream is None:
            self.setParam("ref", ref)
            self.setParam("genome", genome)
            self.setParam("threads", threads)
            self.setParam("memSize", memSize)
            if OutputDir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("vcfInput")[0])),
                )
            else:
                self.setOutput("outputdir", OutputDir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("ref", Configure.getRefDir())
            self.setParam("genome", Configure.getGenome())
            self.setParam("threads", Configure.getThreads())
            self.setParam("memSize", Configure.getJavaMem())

        if (upstream is None) or (upstream is True):
            self.setInput("vcfInput", vcfInput)
            self.setParam(
                "prefix",
                [self.getMaxFileNamePrefixV2(x) for x in self.getInput("vcfInput")],
            )
            self.setOutput(
                "vcfOutput",
                [
                    os.path.join(self.getOutput("outputdir"), x + ".filter.vcf")
                    for x in self.getParam("prefix")
                ],
            )
            self.setOutput(
                "summaryOutput",
                [
                    os.path.join(self.getOutput("outputdir"), x + ".filter.summary.txt")
                    for x in self.getParam("prefix")
                ],
            )
        else:
            # check Configure for running pipeline
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "BisulfiteGenotyper":
                self.setInput("indir", upstream.getOutput("outputdir"))
                self.setParam(
                    "prefix",
                    [os.path.basename(x) for x in upstream.getOutput("outdir")],
                )
            else:
                raise commonError("Parameter upstream must from BisulfiteGenotyper.")

            chromosome = ["chr%i" % x for x in range(1, 23)]
            chromosome.extend(["chrX", "chrY", "chrM"])
            vcfInput = []
            vcfOutput = []
            summaryOutput = []

            for x in self.getParam("prefix"):
                if not os.path.exists(os.path.join(self.getOutput("outputdir"), x)):
                    os.makedirs(os.path.join(self.getOutput("outputdir"), x))

                vcfInput.extend(
                    [
                        os.path.join(self.getInput("indir"), f"{x}/{x}_{y}.snp.raw.vcf")
                        for y in chromosome
                    ]
                )
                vcfOutput.extend(
                    [
                        os.path.join(
                            self.getOutput("outputdir"), f"{x}/{x}_{y}.filter.vcf"
                        )
                        for y in chromosome
                    ]
                )
                summaryOutput.extend(
                    [
                        os.path.join(
                            self.getOutput("outputdir"),
                            f"{x}/{x}_{y}.filter.summary.txt",
                        )
                        for y in chromosome
                    ]
                )

            self.setOutput(
                "outdir",
                [
                    os.path.join(self.getOutput("outputdir"), x)
                    for x in self.getParam("prefix")
                ],
            )
            self.setInput("vcfInput", vcfInput)
            self.setOutput("vcfOutput", vcfOutput)
            self.setOutput("summaryOutput", summaryOutput)

        self.checkInputFilePath()
        self.refcheck()

        if BisSNP:
            self.setInput("jarInput", BisSNP)
        else:
            self.setInput(
                "jarInput",
                pkg_resources.resource_filename("cfDNApipe", "data/BisSNP-0.82.2.jar"),
            )

        all_cmd = []

        for i in range(len(self.getInput("vcfInput"))):
            tmp_cmd = self.cmdCreate(
                [
                    "java",
                    "-Xmx10G",
                    "-jar",
                    self.getInput("jarInput"),
                    "-T",
                    "VCFpostprocess",
                    "-R",
                    self.getParam("ref") + "/" + self.getParam("genome") + ".fa",
                    "-oldVcf",
                    self.getInput("vcfInput")[i],
                    "-newVcf",
                    self.getOutput("vcfOutput")[i],
                    "-snpVcf",
                    self.getInput("vcfInput")[i],
                    "-o",
                    self.getOutput("summaryOutput")[i],
                ]
            )
            all_cmd.append(tmp_cmd)

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

    def refcheck(self,):
        """reference fasta check."""
        fafile = os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa")

        if not os.path.exists(fafile):
            raise commonError("file " + fafile + " do not exist!")
