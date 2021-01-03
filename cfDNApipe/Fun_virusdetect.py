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


class virusdetect(StepBase):
    def __init__(
        self,
        seqInput1=None,
        seqInput2=None,
        ref=None,
        outputdir=None,
        threads=1,
        paired=True,
        other_params={"--min-hitlen": 25},
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for detect virus in sequencing data.
        Note: this function is calling centrifuge, only used for WGS data.

        virusdetect(seqInput1=None, seqInput2=None, ref=None, outputdir=None, threads=1, paired=True,
                    other_params={"-q": True, "-N": 1, "--time": True}, stepNum=None, upstream=None,)
        {P}arameters:
            seqInput1: list, input _1 fastq files.
            seqInput2: list, input _2 fastq files, None for single end.
            ref: bowtie2 reference path.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            paired: True for paired data, False for single end data.
            other_params: dict, other parameters passing to Bismark.
                          "-parameter": True means "-parameter" in command line.
                          "-parameter": 1 means "-parameter 1" in command line.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline. This parameter can be True, which means a new pipeline start.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """

        super(virusdetect, self).__init__(stepNum, upstream)

        # set fastq input
        # set sequencing input
        if (upstream is None) or (upstream is True):
            self.setInput("seq1", seqInput1)
            self.setInput("seq2", seqInput2)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "inputprocess":
                self.setInput("seq1", upstream.getOutput("fq1"))
                self.setInput("seq2", upstream.getOutput("fq2"))
            elif upstream.__class__.__name__ == "adapterremoval":
                self.setInput("seq1", upstream.getOutput("pair1"))
                self.setInput("seq2", upstream.getOutput("pair2"))
            elif upstream.__class__.__name__ == "bowtie2":
                self.setInput("seq1", upstream.getOutput("unmapped-1"))
                self.setInput("seq2", upstream.getOutput("unmapped-2"))
            else:
                raise commonError("Parameter upstream must from inputprocess, adapterremoval or bowtie2.")

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("seq1")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set ref, threads, paired, ref should be set by user
        if ref is not None:
            self.setParam("ref", ref)
        else:
            self.setParam("ref", Configure.getConfig("virus.ref"))
        if upstream is None:

            self.setParam("threads", threads)
            if paired:
                self.setParam("type", "paired")
            else:
                self.setParam("type", "single")
        else:
            self.setParam("threads", Configure.getThreads())
            self.setParam("type", Configure.getType())

        # check reference for centrifuge
        self.ctfrefcheck()

        # paired or single
        if self.getParam("type") == "paired":
            # generate base name
            prefix = []
            for seq1, seq2 in zip(self.getInput("seq1"), self.getInput("seq2")):
                prefix.append(self.getMaxFileNamePrefix(seq1, seq2))
            self.setParam("prefix", prefix)

            self.setParam(
                "outPrefix",
                [os.path.join(self.getOutput("outputdir"), x) for x in self.getParam("prefix")],
            )

            if other_params is None:
                self.setParam("other_params", "")
            else:
                self.setParam("other_params", other_params)

            self.setOutput("output", [x + ".output" for x in self.getParam("outPrefix")])
            self.setOutput("report", [x + ".report" for x in self.getParam("outPrefix")])

            if len(self.getInput("seq1")) == len(self.getInput("seq2")):
                multi_run_len = len(self.getInput("seq1"))
            else:
                raise commonError("Paired end Input files are not consistent.")

            all_cmd = []

            for i in range(multi_run_len):
                tmp_cmd = self.cmdCreate(
                    [
                        "centrifuge",
                        "--threads",
                        self.getParam("threads"),
                        "-x",
                        self.getParam("ref"),
                        self.getParam("other_params"),
                        "-1",
                        self.getInput("seq1")[i],
                        "-2",
                        self.getInput("seq2")[i],
                        "-S",
                        self.getOutput("output")[i],
                        "--report-file",
                        self.getOutput("report")[i],
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
                [os.path.join(self.getOutput("outputdir"), x) for x in self.getParam("prefix")],
            )

            if other_params is None:
                self.setParam("other_params", "")
            else:
                self.setParam("other_params", other_params)

            self.setOutput("output", [x + ".output" for x in self.getParam("outPrefix")])
            self.setOutput("report", [x + ".report" for x in self.getParam("outPrefix")])

            multi_run_len = len(self.getInput("seq1"))

            all_cmd = []

            for i in range(multi_run_len):
                tmp_cmd = self.cmdCreate(
                    [
                        "centrifuge",
                        "--threads",
                        self.getParam("threads"),
                        "-x",
                        self.getParam("ref"),
                        self.getParam("other_params"),
                        "-U",
                        self.getInput("seq1")[i],
                        "-S",
                        self.getOutput("output")[i],
                        "--report-file",
                        self.getOutput("report")[i],
                    ]
                )
                all_cmd.append(tmp_cmd)

        else:
            commonError("Wrong data tpye, must be 'single' or 'paired'!")

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            if verbose:
                self.run(all_cmd)
            else:
                self.multiRun(args=all_cmd, func=None, nCore=1)

        self.stepInfoRec(cmds=[all_cmd], finishFlag=finishFlag)

    # ref check
    def ctfrefcheck(
        self,
    ):
        extension = [".1.cf", ".2.cf", ".3.cf"]
        ctf2Ref = [self.getParam("ref") + x for x in extension]
        for filePath in ctf2Ref:
            if not os.path.exists(filePath):
                raise commonError("Centrifuge index file " + filePath + " don not exist!")
