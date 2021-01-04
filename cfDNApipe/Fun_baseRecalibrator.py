# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
Modify on Sun Apr 26 11:27:32 2020
@author: LY, Shuying He
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, maxCore
import os
from .Configure import Configure
import math

__metaclass__ = type


class BaseRecalibrator(StepBase):
    def __init__(
        self,
        bamInput=None,
        knownSitesDir=None,  # necessary
        outputdir=None,
        stepNum=None,
        upstream=None,
        genome=None,
        ref=None,  # str
        threads=1,
        verbose=True,
        other_params=None,
        **kwargs
    ):
        """
        This function is used for generating recalibration table for Base Quality Score Recalibration using gatk.
        Note: this function is calling gatk BaseRecalibrator.

        BaseRecalibrator(bamInput=None, knownSitesDir=None, 
                outputdir=None, threads=1, genome=None, ref=None,
                stepNum=None, upstream=None, verbose=True)
                
        {P}arameters:
            bamInput: list, input bam files.
            knownSitesDir: str, dirname for knownsites file.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            genome: str, human genome version, just support "hg19" and "hg38"
            ref: str, reference folderpath
            stepNum: int or str, step flag for folder name.
            upstream: upstream output results, used for pipeline, just can be addRG. This parameter can be True, which means a new pipeline start.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
            other_params: dict, other parameters passing to gatk, default is None.
        """

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
        if self.getParam("genome") == "hg19":
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
        elif self.getParam("genome") == "hg38":
            extension = [
                "1000G_omni2.5.hg38.vcf",
                "1000G_phase1.snps.high_confidence.hg38.vcf",
                "dbsnp_146.hg38.vcf",
                "hapmap_3.3.hg38.vcf",
                "Mills_and_1000G_gold_standard.indels.hg38.vcf",
                "1000G_omni2.5.hg38.vcf.idx",
                "1000G_phase1.snps.high_confidence.hg38.vcf.idx",
                "dbsnp_146.hg38.vcf.idx",
                "hapmap_3.3.hg38.vcf.idx",
                "Mills_and_1000G_gold_standard.indels.hg38.vcf.idx",
            ]
        else:
            raise commonError("Wrong genome type! Just be hg19 or hg38...")

        vcffile = [os.path.join(self.getParam("knownSitesDir"), x) for x in extension]

        self.setParam("vcffile", vcffile)

        if other_params:
            self.setParam("other_params", other_params)
        else:
            self.setParam("other_params", "")

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
                    self.getParam("other_params"),
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
                    nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
                )

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)

    def recalcheck(self,):
        fafile = [os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa")]
        vcffile = self.getParam("vcffile")
        filePaths = fafile + vcffile
        for filePath in filePaths:
            if not os.path.exists(filePath):
                raise commonError("file " + filePath + " don not exist!")
