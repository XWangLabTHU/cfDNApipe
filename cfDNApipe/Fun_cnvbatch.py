# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 11:45:21 2020

@author: Shuying He

E-mail: heshuying@fzidt.com
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure

__metaclass__ = type


class cnvbatch(StepBase):
    def __init__(
        self,
        casebamInput=None,
        ctrlbamInput=None,
        outputdir=None,
        threads=1,
        access=None,
        annotate=None,
        ref=None,
        genome=None,
        reference_cnn=None,
        other_params={"-m": "wgs", "-y": True},
        stepNum=None,
        caseupstream=None,
        ctrlupstream=None,
        **kwargs
    ):
        """
        This function is used for cnv calling using cnvkit.
        Note: This function is calling cnvkit.py batch, please install cnvkit before using.

        cnvbatch(casebamInput=None, ctrlbamInput=None,
            outputdir=None, genome=None, ref=None, threads=1,
            access=None, annotate=None, reference_cnn=None,
            other_params={"-m": "wgs", "-y": True},
            stepNum=None, caseupstream=None, ctrlupstream=None,
            verbose=False, **kwargs)

        {P}arameters:
            casebamInput: list, case bam files, only this part of bam files will call cnv.
            casebamInput: list, ctrl bam files, this bam files is using for generating reference cnn, will not call cnv.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            genome: str, human genome version, just support "hg19" and "hg38"
            ref: str, reference folderpath.
            access: str, an "access" file precomputed for the UCSC reference human genome, with some know low-mappability regions excluded. eg: access-5kb-mappable.hg19.bed.
            annotate: str, Gene annotation databases, eg. refFlat_hg19.txt.
            reference_cnn: str, reference cnn file, you can using the ready-made cnn file as baseline.
            stepNum: int or str, step flag for folder name.
            other_params: str or dict. other parameters for cnvkit batch. default is {"-m": "wgs", "-y": True}, which means WGS data, and the ref is male.
            caseupstream: case upstream output results, used for call cnv pipeline. This parameter can be True, which means a new pipeline start.
            ctrlupstream: ctrl upstream output results.
        """

        super(cnvbatch, self).__init__(stepNum, caseupstream)

        if (caseupstream is None) or (caseupstream is True):
            self.setInput("casebamInput", casebamInput)
        else:
            Configure.configureCheck()
            caseupstream.checkFilePath()

            if caseupstream.__class__.__name__ in [
                "rmduplicate",
                "bamsort",
            ]:
                self.setInput("casebamInput", caseupstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from rmduplicate, bamsort.")

        if ctrlupstream is None:
            if ctrlbamInput is not None:
                self.setInput("ctrlbamInput", ctrlbamInput)
                self.checkInputFilePath()
        else:
            ctrlupstream.checkFilePath()

            if ctrlupstream.__class__.__name__ in [
                "rmduplicate",
                "bamsort",
            ]:
                self.setInput("ctrlbamInput", ctrlupstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from rmduplicate, bamsort.")

        self.checkInputFilePath()

        if caseupstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("casebamInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        if ref is None:
            self.setParam("ref", Configure.getRefDir())
        else:
            self.setParam("ref", ref)

        if genome is None:
            self.setParam("genome", Configure.getGenome())
        else:
            self.setParam("genome", genome)

        if threads is None:
            self.setParam("threads", Configure.getThreads())
        else:
            self.setParam("threads", threads)

        if reference_cnn is not None:
            self.setInput("reference_cnn", reference_cnn)
        else:
            self.setOutput(
                "reference_cnn",
                os.path.join(self.getOutput("outputdir"), "reference.cnn"),
            )

            self.setInput("access", access)
            self.setInput("annotate", annotate)
            self.refcheck()

        self.setParam("other_params", other_params)

        # Output cnr and cns
        self.setOutput(
            "cnsOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"),
                    os.path.basename(x).replace(".bam", "") + ".cns",
                )
                for x in self.getInput("casebamInput")
            ],
        )

        self.setOutput(
            "cnrOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"),
                    os.path.basename(x).replace(".bam", "") + ".cnr",
                )
                for x in self.getInput("casebamInput")
            ],
        )

        if "reference_cnn" in self.getInputs():
            # Reusing a reference for additional samples
            cmd = self.cmdCreate(
                [
                    "cnvkit.py",
                    "batch",
                    " ".join(self.getInput("casebamInput")),
                    "-r",
                    self.getInput("reference_cnn"),
                    "-d",
                    self.getOutput("outputdir"),
                    "-p",
                    self.getParam("threads"),
                ]
            )
        else:
            if "ctrlbamInput" in self.getInputs():
                normal = " ".join(self.getInput("ctrlbamInput"))
            else:
                normal = ""

            # classic cnv pipeline for with or without ctrlbamInput
            cmd = self.cmdCreate(
                [
                    "cnvkit.py",
                    "batch",
                    " ".join(self.getInput("casebamInput")),
                    "--normal",
                    normal,
                    "--annotate",
                    self.getInput("annotate"),
                    "--access",
                    self.getInput("access"),
                    "--fasta",
                    os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa"),
                    self.getParam("other_params"),
                    "--output-reference",
                    self.getOutput("reference_cnn"),
                    "--output-dir",
                    self.getOutput("outputdir"),
                    "-p",
                    self.getParam("threads"),
                ]
            )

        finishFlag = self.stepInit(caseupstream)

        if not finishFlag:
            self.run(cmd)

        self.stepInfoRec(cmds=[cmd], finishFlag=finishFlag)

    def refcheck(
        self,
    ):
        """check reference FilePath."""
        fafile = os.path.join(self.getParam("ref"), self.getParam("genome") + ".fa")

        if not os.path.exists(fafile):
            raise commonError("file " + fafile + " don not exist!")
