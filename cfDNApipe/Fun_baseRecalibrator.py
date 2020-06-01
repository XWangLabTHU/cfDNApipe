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


class BaseRecalibrator(StepBase):
    def __init__(
        self,
        bamInput=None,
        knownSitesDir=None,
        outputdir=None,
        stepNum=None,
        upstream=None,
        genome=None,
        ref=None,  # str
        threads=1,
        verbose=True,
        **kwargs
    ):

        super(BaseRecalibrator, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
            # In this situation, input file and output path should be checked
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "addRG":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from addRG.")

        self.checkInputFilePath()

        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)

            self.setParam("ref", ref)
            self.setParam("genome", genome)
            self.setParam("threads", threads)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("ref", Configure.getRefDir())
            self.setParam("genome", Configure.getGenome())
            self.setParam("threads", Configure.getThreads())

        self.setParam("knownSitesDir", knownSitesDir)
        extension = [
            "1000G_omni2.5.hg19.sites.vcf",
            "1000G_phase1.indels.hg19.sites.vcf",
            "1000G_phase1.snps.high_confidence.hg19.sites.vcf",
            "dbsnp_138.hg19.vcf",
            "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf",
            "1000G_omni2.5.hg19.sites.vcf.idx",
            "1000G_phase1.indels.hg19.sites.vcf.idx",
            "1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx",
            "dbsnp_138.hg19.vcf.idx",
            "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx",
        ]

        vcffile = [os.path.join(self.getParam("knownSitesDir"), x) for x in extension]

        self.setParam("vcffile", vcffile)
        self.recalcheck()
        self.setOutput(
            "recalOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + ".recal.table"
                for x in self.getInput("bamInput")
            ],
        )

        self.setOutput("bamOutput", self.getInput("bamInput"))

        multi_run_len = len(self.getInput("bamInput"))

        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(
                [
                    "gatk",
                    "BaseRecalibrator ",
                    "-I",
                    self.getInput("bamInput")[i],
                    "-R",
                    self.getParam("ref") + "/" + self.getParam("genome") + ".fa",
                    "--known-sites",
                    self.getParam("vcffile")[0],
                    "--known-sites",
                    self.getParam("vcffile")[1],
                    "--known-sites",
                    self.getParam("vcffile")[2],
                    "--known-sites",
                    self.getParam("vcffile")[3],
                    "--known-sites",
                    self.getParam("vcffile")[4],
                    "-O",
                    self.getOutput("recalOutput")[i],
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

    def recalcheck(self,):

        fafile = [os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa")]
        vcffile = self.getParam("vcffile")
        filePaths = fafile + vcffile
        for filePath in filePaths:
            if not os.path.exists(filePath):
                raise commonError("file " + filePath + " don not exist!")
