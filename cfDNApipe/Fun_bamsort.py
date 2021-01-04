# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 11:45:21 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure

__metaclass__ = type


class bamsort(StepBase):
    def __init__(
        self,
        bamInput=None,
        outputdir=None,
        threads=1,
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for sorting bam files.
        Note: this function is calling samtools.

        bamsort(bamInput=None, outputdir=None, threads=1, stepNum=None, upstream=None)
        {P}arameters:
            bamInput: list, bam file input.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            stepNum: int or str, step flag for folder name.
            upstream: upstream output results, used for pipeline.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """

        super(bamsort, self).__init__(stepNum, upstream)

        # set fastq input
        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ in [
                "bowtie2",
                "bismark",
                "bismark_deduplicate",
            ]:
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from bowtie2 or bismark.")

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set threads
        if upstream is None:
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure.getThreads())

        self.setOutput(
            "bamOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "_sorted.bam"
                for x in self.getInput("bamInput")
            ],
        )

        self.setOutput("baiOutput", [x + ".bai" for x in self.getOutput("bamOutput")])

        multi_run_len = len(self.getInput("bamInput"))

        cmd1 = []
        cmd2 = []

        for i in range(multi_run_len):
            cmd1.append(
                self.cmdCreate(
                    [
                        "samtools sort",
                        "-@",
                        self.getParam("threads"),
                        "-o",
                        self.getOutput("bamOutput")[i],
                        self.getInput("bamInput")[i],
                    ]
                )
            )
            cmd2.append(
                self.cmdCreate(
                    [
                        "samtools index",
                        "-@",
                        self.getParam("threads"),
                        self.getOutput("bamOutput")[i],
                    ]
                )
            )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            if verbose:
                print("Sorting BAM files......")
                self.run(cmd1)
                print("Indexing BAM files......")
                self.run(cmd2)
            else:
                print("Sorting BAM files......")
                self.multiRun(args=cmd1, func=None, nCore=1)
                print("Indexing BAM files......")
                self.multiRun(args=cmd2, func=None, nCore=1)

        self.stepInfoRec(cmds=[cmd1, cmd2], finishFlag=finishFlag)
