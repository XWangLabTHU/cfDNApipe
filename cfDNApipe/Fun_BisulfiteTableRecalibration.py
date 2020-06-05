# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 18:27:32 2020
@author: He Shuying
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import math
import pkg_resources

__metaclass__ = type


class BisulfiteTableRecalibration(StepBase):
    def __init__(
        self,
        BisSNP=None,  # necessary BisSNP.jar
        bamInput=None,
        CovFile=None,
        memSize="2g",
        OutputDir=None,
        genome=None,
        ref=None,  # str
        stepNum=None,
        upstream=None,
        threads=1,
        verbose=False,
        **kwargs
    ):

        super(BisulfiteTableRecalibration, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
            self.setInput("CovFile", CovFile)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "BisulfiteCountCovariates":
                self.setInput("bamInput", upstream.getInput("bamInput"))
                self.setInput("CovFile", upstream.getOutput("CovOutput"))
            else:
                raise commonError(
                    "Parameter upstream must from BisulfiteCountCovariates."
                )
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

        if BisSNP:
            self.setInput("jarInput", BisSNP)
        else:
            self.setInput(
                "jarInput",
                pkg_resources.resource_filename("cfDNApipe", "data/BisSNP-0.82.2.jar"),
            )

        self.refcheck()
        self.setOutput(
            "bamOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"),
                    self.getMaxFileNamePrefixV2(x) + "-BQSR.bam",
                )
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
                    "-T",
                    "BisulfiteTableRecalibration",
                    "-R",
                    self.getParam("ref") + "/" + self.getParam("genome") + ".fa",
                    "-I",
                    self.getInput("bamInput")[i],
                    "-recalFile",
                    self.getInput("CovFile")[i],
                    "-o",
                    self.getOutput("bamOutput")[i],
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
        """check ref file exists."""
        fafile = os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa")

        if not os.path.exists(fafile):
            raise commonError("file " + fafile + " do not exist!")
