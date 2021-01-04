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


class dbimport(StepBase):
    def __init__(
        self,
        vcfInput=None,  # vcffiles
        outputdir=None,
        genome=None,
        ref=None,  # str
        stepNum=None,
        upstream=None,
        threads=1,
        verbose=False,
        **kwargs
    ):
        """
        This function is to import the normal samples' mutation VCF into GenomicsDB using gatk.
        Note: This function is calling gatk GenomicsDBImport, please install gatk before using.

        dbimport(vcfInput=None, outputdir=None,
            genome=None, ref=None, threads=1,
            stepNum=None, upstream=None, verbose=False, **kwargs)

        {P}arameters:
            vcfInput: list, vcf files.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            genome: str, human genome version, just support "hg19" and "hg38"
            ref: str, reference folderpath.
            stepNum: int or str, step flag for folder name.
            upstream: upstream output results, used for call snp pipeline, just can be mutect2n. This parameter can be True, which means a new pipeline start.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """

        super(dbimport, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
            self.setInput("vcfInput", vcfInput)
        else:
            # check Configure for running pipeline
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "mutect2n":
                self.setInput("vcfInput", upstream.getOutput("vcfOutput"))
            else:
                raise commonError("Parameter upstream must from mutect2n.")

        self.checkInputFilePath()

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

        self.dbimportcheck()

        all_cmd = []

        chromosome = [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
            "chrX",
            "chrY",
            "chrM",
        ]

        # all vcf input
        input_vcfs = " -V ".join(self.getInput("vcfInput"))

        self.setOutput(
            "dbimportOutput",
            [os.path.join(self.getOutput("outputdir"), "pon_db_" + x) for x in chromosome],
        )

        for x in range(len(chromosome)):
            tmp_cmd = self.cmdCreate(
                [
                    "gatk",
                    "GenomicsDBImport",
                    "-V",
                    input_vcfs,
                    "-R",
                    self.getParam("ref") + "/" + self.getParam("genome") + ".fa",
                    "-L",
                    chromosome[x],
                    "--genomicsdb-workspace-path",
                    self.getOutput("dbimportOutput")[x],
                    "--java-options",
                    "'-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'",
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

    def dbimportcheck(
        self,
    ):
        fafile = os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa")

        if not os.path.exists(fafile):
            raise commonError("file " + fafile + " don not exist!")
