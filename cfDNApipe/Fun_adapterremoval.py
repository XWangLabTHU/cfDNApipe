# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:58:34 2019

@author: zhang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
import re
from .Configure import Configure

__metaclass__ = type


class adapterremoval(StepBase):
    def __init__(
        self,
        fqInput1=None,
        fqInput2=None,
        outputdir=None,
        threads=1,
        paired=True,
        adapter1=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
        adapter2=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"],
        other_params={"--qualitybase": 33, "--gzip": True},
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for removing adapters in paired end fastq files.
        Note: this function is calling AdapterRemoval.

        adapterremoval(fqInput1=None, fqInput2=None, outputdir=None, threads=1, paired=True,
                       adapter1=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
                       adapter2=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"],
                       other_params={"--qualitybase": 33, "--gzip": True},
                       stepNum=None, upstream=None, verbose=True)
        {P}arameters:
            fqInput1: list, fastq 1 files or single end files.
            fqInput2: list, fastq 2 files or None for single end data.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            paired: boolean, paired end or single end.
            adapter1: list, adapters for corresponding fqInput1.
            adapter2: list, adapters for corresponding fqInput2, None for single end data.
            other_params: dict, other parameters passing to command "AdapterRemoval".
                          "-parameter": True means "-parameter" in command line.
                          "-parameter": 1 means "-parameter 1" in command line.
            stepNum: int, step number for folder name.
            upstream: upstream output results, used for pipeline.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """

        super(adapterremoval, self).__init__(stepNum, upstream)

        # set fastq input
        if (upstream is None) or (upstream is True):
            self.setInput("fq1", fqInput1)
            self.setInput("fq2", fqInput2)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            self.setInput("fq1", upstream.getOutput("fq1"))
            self.setInput("fq2", upstream.getOutput("fq2"))

        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(self.getInput("fq1")[1])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set threads, paired and adapters
        if upstream is None:
            if paired:
                self.setParam("type", "paired")
            else:
                self.setParam("type", "single")
            self.setParam("threads", threads)
            self.setParam("adapter1", adapter1)
            self.setParam("adapter2", adapter2)

        else:
            self.setParam("type", Configure.getType())
            self.setParam("threads", Configure.getThreads())
            if upstream.__class__.__name__ == "identifyAdapter":
                # using identified adapters
                adapter1 = []
                adapter2 = []
                for fq1, fq2 in zip(self.getInput("fq1"), self.getInput("fq2")):
                    basename = self.getMaxFileNamePrefix(fq1, fq2)
                    adapter_file = upstream.getOutput(basename + "-adapterFile")
                    adapters = self.getAdapetrFromFile(adapter_file)
                    adapter1.append(adapters[0])
                    adapter2.append(adapters[1])
            elif upstream.__class__.__name__ == "inputprocess":
                if len(adapter1) != len(self.getInput("fq1")) and len(adapter1) == 1:
                    adapter1 = adapter1 * len(self.getInput("fq1"))
                if len(adapter2) != len(self.getInput("fq2")) and len(adapter2) == 1:
                    adapter2 = adapter2 * len(self.getInput("fq2"))
            else:
                raise commonError("Parameter upstream must from inputprocess or identifyAdapter.")

            self.setParam("adapter1", adapter1)
            self.setParam("adapter2", adapter2)

        if self.getParam("type") == "paired":
            basenames = []
            for fq1, fq2 in zip(self.getInput("fq1"), self.getInput("fq2")):
                basenames.append(self.getMaxFileNamePrefix(fq1, fq2))
            self.setParam("basename", basenames)

            self.setParam(
                "outPrefix", [os.path.join(self.getOutput("outputdir"), x) for x in self.getParam("basename")],
            )

            if other_params is None:
                self.setParam("other_params", "")
            else:
                self.setParam("other_params", other_params)

            if len(self.getInput("fq1")) == len(self.getInput("fq2")):
                multi_run_len = len(self.getInput("fq1"))
            else:
                raise commonError("Paired end Input files are not consistent.")

            self.setOutput(
                "discarded",
                [(os.path.join(self.getOutput("outputdir"), x) + ".discarded.gz") for x in self.getParam("basename")],
            )
            self.setOutput(
                "settings",
                [(os.path.join(self.getOutput("outputdir"), x) + ".settings") for x in self.getParam("basename")],
            )
            self.setOutput(
                "pair1",
                [
                    (os.path.join(self.getOutput("outputdir"), x) + ".pair1.truncated.gz")
                    for x in self.getParam("basename")
                ],
            )
            self.setOutput(
                "pair2",
                [
                    (os.path.join(self.getOutput("outputdir"), x) + ".pair2.truncated.gz")
                    for x in self.getParam("basename")
                ],
            )
            self.setOutput(
                "singleton",
                [
                    (os.path.join(self.getOutput("outputdir"), x) + ".singleton.truncated.gz")
                    for x in self.getParam("basename")
                ],
            )

            all_cmd = []

            for i in range(multi_run_len):
                tmp_cmd = self.cmdCreate(
                    [
                        "AdapterRemoval",
                        "--threads",
                        self.getParam("threads"),
                        "--file1",
                        self.getInput("fq1")[i],
                        "--file2",
                        self.getInput("fq2")[i],
                        "--adapter1",
                        self.getParam("adapter1")[i],
                        "--adapter2",
                        self.getParam("adapter2")[i],
                        "--basename",
                        self.getParam("outPrefix")[i],
                        self.getParam("other_params"),
                    ]
                )
                all_cmd.append(tmp_cmd)

        elif self.getParam("type") == "single":
            self.setParam(
                "basename", [self.getMaxFileNamePrefixV2(x) for x in self.getInput("fq1")],
            )
            self.setParam(
                "outPrefix", [os.path.join(self.getOutput("outputdir"), x) for x in self.getParam("basename")],
            )

            if other_params is None:
                self.setParam("other_params", "")
            else:
                self.setParam("other_params", other_params)

            multi_run_len = len(self.getInput("fq1"))

            self.setOutput(
                "discarded",
                [(os.path.join(self.getOutput("outputdir"), x) + ".discarded.gz") for x in self.getParam("basename")],
            )
            self.setOutput(
                "settings",
                [(os.path.join(self.getOutput("outputdir"), x) + ".settings") for x in self.getParam("basename")],
            )
            self.setOutput(
                "pair1",
                [(os.path.join(self.getOutput("outputdir"), x) + ".truncated.gz") for x in self.getParam("basename")],
            )
            self.setOutput("pair2", [])
            self.setOutput("singleton", [])

            all_cmd = []

            for i in range(multi_run_len):
                tmp_cmd = self.cmdCreate(
                    [
                        "AdapterRemoval",
                        "--threads",
                        self.getParam("threads"),
                        "--file1",
                        self.getInput("fq1")[i],
                        "--adapter1",
                        self.getParam("adapter1")[i],
                        "--basename",
                        self.getParam("outPrefix")[i],
                        self.getParam("other_params"),
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

    # get adapter from file
    def getAdapetrFromFile(self, file):
        adapter = []
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if re.search(r"Consensus:", line):
                    tmp_line = line.split()
                    adapter.append(tmp_line[1])

        return adapter
