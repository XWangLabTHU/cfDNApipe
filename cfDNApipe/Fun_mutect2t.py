# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
Modify on Wed Feb 26 18:27:32 2020
@author: LY, He Shuying
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, maxCore
import os
from .Configure import Configure
import math

__metaclass__ = type


class mutect2t(StepBase):
    def __init__(
        self,
        bamInput=None,
        ponbedInput=None,
        vcfInput=None,
        outputdir=None,
        genome=None,
        ref=None,
        threads=1,
        stepNum=None,
        caseupstream=None,
        ctrlupstream=None,
        verbose=False,
        **kwargs
    ):
        """
        This function is used for Call somatic SNVs and indels via local assembly of haplotypes using gatk.
        Note: This function is calling gatk Mutect2, please install gatk before using.

        mutect2t(bamInput=None, ponbedInput=None, vcfInput=None,
            outputdir=None, genome=None, ref=None, threads=1,
            stepNum=None, other_params=None, caseupstream=None,
            ctrlupstream=None, verbose=False, **kwargs)

        {P}arameters:
            bamInput: list, bam files.
            ponbedInput: str, panel-of-normals file, you can generating this file from mutect2n or using download   file (like somatic-hg38_1000g_pon.hg38.vcf.gz...).
            vcfInput: str, like af-only-gnomad.raw.sites.hg19.vcf.gz. this is for gatk parameter "--germline-resource".
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            genome: str, human genome version, just support "hg19" and "hg38"
            ref: str, reference folderpath.
            stepNum: int or str, step flag for folder name.
            other_params: str or dict. other parameters for gatk mutect, default is None.
            caseupstream: case upstream output results, used for call snp pipeline, just can be contamination / BQSR / addRG. This parameter can be True, which means a new pipeline start.
            ctrlupstream: ctrl upstream output results, used for pon bed file, just can be createPON. This parameter can be None, which means you need to give me ponbedInput.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """

        super(mutect2t, self).__init__(stepNum, caseupstream)

        if (caseupstream is None) or (caseupstream is True):
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            caseupstream.checkFilePath()

            if caseupstream.__class__.__name__ == "contamination":
                self.setInput("bamInput", caseupstream.getOutput("bamOutput"))
                self.setOutput("contaminationOutput", caseupstream.getOutput("contaminationOutput"))
            elif caseupstream.__class__.__name__ in ["BQSR", "addRG"]:
                self.setInput("bamInput", caseupstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter case upstream must from contamination or BQSR or addRG.")

        # PON file
        if ctrlupstream and ctrlupstream.__class__.__name__ == "createPON":
            ctrlupstream.checkFilePath()
            self.setInput("ponbedInput", ctrlupstream.getOutput("createPONOutput"))
        else:
            ponbed = self.convertToList(ponbedInput)
            if len(ponbed) == 1:
                ponbed = ponbed * 25
            self.setInput("ponbedInput", ponbed)

        if self.getInput("ponbedInput") is None:
            raise commonError("ponbedInput cannot be none.")

        self.checkInputFilePath()

        # set outputdir / ref / threads
        if caseupstream is None:
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

        self.mutect2tcheck()
        self.setInput("vcfInput", vcfInput)
        self.setParam(
            "prefix",
            [self.getMaxFileNamePrefixV2(x) for x in self.getInput("bamInput")],
        )

        self.setOutput(
            "outdir",
            [os.path.join(self.getOutput("outputdir"), z) for z in self.getParam("prefix")],
        )

        chromosome = ["chr%i" % x for x in range(1, 23)]
        chromosome.extend(["chrX", "chrY", "chrM"])

        bamnum = len(self.getParam("prefix"))
        all_cmd = []
        for i in range(bamnum):
            for y in range(len(chromosome)):
                tmp_cmd = self.cmdCreate(
                    [
                        "gatk",
                        "Mutect2",
                        "-R",
                        self.getParam("ref") + "/" + self.getParam("genome") + ".fa",
                        "-I",
                        self.getInput("bamInput")[i],
                        "--germline-resource",
                        self.getInput("vcfInput"),
                        "--panel-of-normals",
                        self.getInput("ponbedInput")[y],
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

        finishFlag = self.stepInit(caseupstream)

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
                    nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
                )

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)

    def mutect2tcheck(
        self,
    ):

        fafile = os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa")

        if not os.path.exists(fafile):
            raise commonError("file " + fafile + " don not exist!")
