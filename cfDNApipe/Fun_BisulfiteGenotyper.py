# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 18:27:32 2019
@author: He Shuying
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import math
import pkg_resources

__metaclass__ = type


class BisulfiteGenotyper(StepBase):
    def __init__(
        self,
        BisSNP=None,  # BisSNP
        bamInput=None,
        knownsite=None,
        memSize="2g",
        OutputDir=None,
        genome=None,
        ref=None,  # str
        threads=1,
        Other_Params="-C CG,1 -C CH,1 -out_modes EMIT_VARIANT_AND_CYTOSINES -stand_call_conf 10 -stand_emit_conf 0 -nt 1 -mmq 20 -mbq 5",
        stepNum=None,
        upstream=None,
        verbose=False,
        **kwargs
    ):

        super(BisulfiteGenotyper, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
        else:
            # check Configure for running pipeline
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "BisulfiteTableRecalibration" or "addRG":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError(
                    "Parameter upstream must from BisulfiteTableRecalibration or addRG."
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

        self.refcheck()
        if BisSNP:
            self.setInput("jarInput", BisSNP)
        else:
            self.setInput(
                "jarInput",
                pkg_resources.resource_filename("cfDNApipe", "data/BisSNP-0.82.2.jar"),
            )

        self.setInput("knownsite", knownsite)
        self.setParam(
            "prefix",
            [self.getMaxFileNamePrefixV2(x) for x in self.getInput("bamInput")],
        )
        self.setParam("other_params", Other_Params)
        self.setOutput(
            "outdir",
            [
                os.path.join(self.getOutput("outputdir"), z)
                for z in self.getParam("prefix")
            ],
        )

        # chromosome set
        chromosome = ["chr%i" % x for x in range(1, 23)]
        chromosome.extend(["chrX", "chrY", "chrM"])

        bamnum = len(self.getInput("bamInput"))
        all_cmd = []
        for i in range(bamnum):
            for y in chromosome:
                tmp_cmd = self.cmdCreate(
                    [
                        "java",
                        "-Xmx%s" % self.getParam("memSize"),
                        "-jar",
                        self.getInput("jarInput"),
                        "-T",
                        "BisulfiteGenotyper",
                        "-R",
                        self.getParam("ref") + "/" + self.getParam("genome") + ".fa",
                        "-I",
                        self.getInput("bamInput")[i],
                        "-D",
                        self.getInput("knownsite"),
                        "-L",
                        y,
                        "-vfn1",
                        self.getOutput("outdir")[i]
                        + "/"
                        + self.getParam("prefix")[i]
                        + "_"
                        + y
                        + ".cpg.raw.vcf",
                        "-vfn2",
                        self.getOutput("outdir")[i]
                        + "/"
                        + self.getParam("prefix")[i]
                        + "_"
                        + y
                        + ".snp.raw.vcf",
                        self.getParam("other_params"),
                    ]
                )
                all_cmd.append(tmp_cmd)

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            for i in range(bamnum):
                if not os.path.exists(self.getOutput("outdir")[i]):
                    os.mkdir(self.getOutput("outdir")[i])
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
        """check ref file exists? """
        fafile = os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa")

        if not os.path.exists(fafile):
            raise commonError("file " + fafile + " do not exist!")
