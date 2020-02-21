# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 09:28:19 2019

@author: zhang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure


__metaclass__ = type


class bowtie2(StepBase):
    def __init__(
        self,
        seqInput1=None,  # list
        seqInput2=None,  # list
        ref=None,  # str
        outputdir=None,  # str
        threads=1,
        paired=True,
        other_params={"-q": True, "-N": 1, "--time": True},
        stepNum=None,
        upstream=None,
        **kwargs
    ):
        super(bowtie2, self).__init__(stepNum, upstream)
        if upstream is None:
            self.setInput("seq1", seqInput1)
            self.setInput("seq2", seqInput2)
            self.checkInputFilePath()

            if paired:
                self.setParam("type", "paired")
            else:
                self.setParam("type", "single")

            self.setParam("ref", ref)

            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("fq1")[1])),
                )
            else:
                self.setOutput("outputdir", outputdir)

            self.setParam("threads", threads)

        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            self.setParam("type", Configure.getType())

            if upstream.__class__.__name__ == "inputprocess":
                self.setInput("seq1", upstream.getOutput("fq1"))
                self.setInput("seq2", upstream.getOutput("fq2"))
            elif upstream.__class__.__name__ == "adapterremoval":
                self.setInput("seq1", upstream.getOutput("pair1"))
                self.setInput("seq2", upstream.getOutput("pair2"))
            else:
                raise commonError(
                    "Parameter upstream must from inputprocess or adapterremoval."
                )

            self.setParam(
                "ref", os.path.join(Configure.getRefDir(), Configure.getGenome())
            )
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())

        # check reference for bowtie2
        self.bt2refcheck()

        if self.getParam("type") == "paired":
            # generate base name
            prefix = []
            for seq1, seq2 in zip(self.getInput("seq1"), self.getInput("seq2")):
                prefix.append(self.getMaxFileNamePrefix(seq1, seq2))
            self.setParam("prefix", prefix)

            self.setParam(
                "outPrefix",
                [
                    os.path.join(self.getOutput("outputdir"), x)
                    for x in self.getParam("prefix")
                ],
            )

            if other_params is None:
                self.setParam("other_params", "")
            else:
                self.setParam("other_params", other_params)

            self.setParam(
                "unmapped", [x + ".unmapped.gz" for x in self.getParam("outPrefix")]
            )
            self.setOutput(
                "unmapped-1", [x + ".unmapped.1.gz" for x in self.getParam("outPrefix")]
            )
            self.setOutput(
                "unmapped-2", [x + ".unmapped.2.gz" for x in self.getParam("outPrefix")]
            )
            self.setOutput(
                "bamOutput", [x + ".bam" for x in self.getParam("outPrefix")]
            )

            if len(self.getInput("seq1")) == len(self.getInput("seq1")):
                multi_run_len = len(self.getInput("seq1"))
            else:
                raise commonError("Paired end Input files are not consistent.")

            all_cmd = []

            for i in range(multi_run_len):
                tmp_cmd = self.cmdCreate(
                    [
                        "bowtie2",
                        "-x",
                        self.getParam("ref"),
                        "-1",
                        self.getInput("seq1")[i],
                        "-2",
                        self.getInput("seq2")[i],
                        self.getParam("other_params"),
                        "--un-conc-gz",
                        self.getParam("unmapped")[i],
                        "-p",
                        self.getParam("threads"),
                        "|",
                        "samtools view -b -S -@",
                        self.getParam("threads"),
                        "-o",
                        self.getOutput("bamOutput")[i],
                        "-",
                    ]
                )
                all_cmd.append(tmp_cmd)

        elif self.getParam("type") == "single":
            self.setParam(
                "prefix",
                [self.getMaxFileNamePrefixV2(x) for x in self.getInput("seq1")],
            )
            self.setParam(
                "outPrefix",
                [
                    os.path.join(self.getOutput("outputdir"), x)
                    for x in self.getParam("prefix")
                ],
            )

            if other_params is None:
                self.setParam("other_params", "")
            else:
                self.setParam("other_params", other_params)

            if other_params is None:
                self.setParam("other_params", "")
            else:
                self.setParam("other_params", other_params)

            self.setOutput(
                "bamOutput", [x + ".bam" for x in self.getParam("outPrefix")]
            )

            multi_run_len = len(self.getInput("seq1"))

            all_cmd = []

            for i in range(multi_run_len):
                tmp_cmd = self.cmdCreate(
                    [
                        "bowtie2",
                        "-x",
                        self.getParam("ref"),
                        "-U",
                        self.getInput("seq1")[i],
                        self.getParam("other_params"),
                        "-p",
                        self.getParam("threads"),
                        "|",
                        "samtools view -b -S -@",
                        self.getParam("threads"),
                        "-o",
                        self.getOutput("bamOutput")[i],
                        "-",
                    ]
                )
                all_cmd.append(tmp_cmd)

        else:
            commonError("Wrong data tpye, must be 'single' or 'paired'!")

        self.setParam("cmd", all_cmd)

        finishFlag = self.stepInit(upstream)

        self.excute(finishFlag)

    # ref check

    def bt2refcheck(self,):
        extension = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
        bt2Ref = [self.getParam("ref") + x for x in extension]
        for filePath in bt2Ref:
            if not os.path.exists(filePath):
                raise commonError("Bowtie2 index file " + filePath + " don not exist!")
