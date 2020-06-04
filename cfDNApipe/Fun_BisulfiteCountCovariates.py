# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
@author: He Shuying
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import pkg_resources

__metaclass__ = type


class BisulfiteCountCovariates(StepBase):
    def __init__(
        self,
        BisSNP=None,  # necessary, eg. ~/soft/BisSNP-0.82.2.jar
        bamInput=None,
        memSize="2g",
        knownSitesDir=None,  # necessary
        OutputDir=None,
        threads=1,
        Other_Params="-cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate",
        genome=None,
        ref=None,  # str
        stepNum=None,
        upstream=None,
        **kwargs
    ):

        super(BisulfiteCountCovariates, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
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
            # In this situation, input file and output path should be checked
            self.setParam("ref", ref)
            self.setParam("genome", genome)
            self.setParam("threads", threads)
            self.setParam("memSize", memSize)
            if OutputDir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", OutputDir)

        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("ref", Configure.getRefDir())
            self.setParam("genome", Configure.getGenome())
            self.setParam("threads", Configure.getThreads())
            self.setParam("memSize", Configure.getJavaMem())

        self.setParam("knownSitesDir", knownSitesDir)
        extension = [
            "/dbsnp_135.hg19.sort.vcf",
            "/1000G_phase1.indels.hg19.sort.vcf",
            "/Mills_and_1000G_gold_standard.indels.hg19.sites.sort.vcf",
        ]
        vcffile = [self.getParam("knownSitesDir") + x for x in extension]
        self.setInput("vcffile", vcffile)
        self.recalcheck()

        if BisSNP:
            self.setInput("jarInput", BisSNP)
        else:
            self.setInput(
                "jarInput",
                pkg_resources.resource_filename("cfDNApipe", "data/BisSNP-0.82.2.jar"),
            )

        self.setParam("other_params", Other_Params)
        self.setOutput(
            "CovOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + ".recalFile_before.csv"
                for x in self.getInput("bamInput")
            ],
        )

        multi_run_len = len(self.getInput("bamInput"))

        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(
                [
                    "java",
                    "-Xmx%s" % self.getParam("memSize"),
                    "-jar",
                    self.getInput("jarInput"),
                    "-R",
                    self.getParam("ref") + "/" + self.getParam("genome") + ".fa",
                    "-T BisulfiteCountCovariates",
                    "-I",
                    self.getInput("bamInput")[i],
                    "--knownSites",
                    self.getInput("vcffile")[0],
                    "--knownSites",
                    self.getInput("vcffile")[1],
                    "--knownSites",
                    self.getInput("vcffile")[2],
                    self.getParam("other_params"),
                    "-recalFile",
                    self.getOutput("CovOutput")[i],
                    "-nt",
                    self.getParam("threads"),
                ]
            )

            all_cmd.append(tmp_cmd)

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            self.run(all_cmd)

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)

    def recalcheck(self,):
        """Check ref file and known-sites vcfs exists"""
        fafile = [os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa")]
        vcffile = self.getInput("vcffile")
        filePaths = fafile + vcffile
        for filePath in filePaths:
            if not os.path.exists(filePath):
                raise commonError("file " + filePath + " do not exist!")
