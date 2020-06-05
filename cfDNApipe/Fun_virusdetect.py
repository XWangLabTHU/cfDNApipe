# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
Modify on Mon Apr 17 16:07:09 2020
@author: LY, HeShuying
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import pkg_resources
import re

__metaclass__ = type


class virusdetect(StepBase):
    def __init__(
        self,
        plInput=None,  # option
        InputDir=None,  # list
        OutputDir=None,
        virusDB=None,
        blastnIdxH=None,
        blastnIdxV=None,
        virusIDfile=None,  # option, for adding the virus note info
        stepNum=None,
        upstream=None,
        threads=1,
        verbose=False,
        **kwargs
    ):

        super(virusdetect, self).__init__(stepNum, upstream)

        if (upstream is None) or (upstream is True):
            # In this situation, input file and output path should be checked
            self.setInput("unmapped-1", [x + "/unmapped.1.fa" for x in InputDir])
            self.setInput("unmapped-2", [x + "/unmapped.2.fa" for x in InputDir])
            self.setParam(
                "prefix", [os.path.basename(os.path.abspath(x)) for x in InputDir]
            )
        else:
            # check Configure for running pipeline
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "unmapfasta":
                self.setInput("unmapped-1", upstream.getOutput("unmapped-1"))
                self.setInput("unmapped-2", upstream.getOutput("unmapped-2"))
                self.setParam(
                    "prefix",
                    [
                        os.path.basename(os.path.dirname(x))
                        for x in upstream.getOutput("unmapped-1")
                    ],
                )
            else:
                raise commonError("Parameter upstream must from unmapfasta.")
        self.checkInputFilePath()

        if upstream is None:
            self.setParam("threads", threads)

            if not OutputDir:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(InputDir[0]))
                )
            else:
                self.setOutput("outputdir", OutputDir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())

        self.setOutput(
            "outdir",
            [
                os.path.abspath(os.path.join(self.getOutput("outputdir"), x))
                for x in self.getParam("prefix")
            ],
        )

        if plInput:
            self.setInput("plInput", plInput)
        else:
            self.setInput(
                "plInput",
                pkg_resources.resource_filename(
                    "cfDNApipe", "data/VirusFinder2.0Plus/detect_virusPlus.pl"
                ),
            )

        self.setInput("virusDB", virusDB)
        self.setParam("blastnIdxH", blastnIdxH)
        self.setParam("blastnIdxV", blastnIdxV)

        # self.setInput('virusIDfile', virusIDfile)

        self.setOutput("virusDB", self.getInput("virusDB"))
        self.setOutput(
            "faOutput", [x + "/results-virus-top1.fa" for x in self.getOutput("outdir")]
        )
        self.setOutput(
            "viruslistOutput", [x + "/virus-list.txt" for x in self.getOutput("outdir")]
        )
        self.setOutput(
            "virusOutput", [x + "/virus.txt" for x in self.getOutput("outdir")]
        )

        all_cmd = []

        for i in range(len(self.getOutput("outdir"))):
            tmp_cmd = self.cmdCreate(
                [
                    "perl",
                    self.getInput("plInput"),
                    "--fa1",
                    self.getInput("unmapped-1")[i],
                    "--fa2",
                    self.getInput("unmapped-2")[i],
                    "--output",
                    self.getOutput("outdir")[i],
                    "--virusDB",
                    self.getInput("virusDB"),
                    "--blastnIdxH",
                    self.getParam("blastnIdxH"),
                    "--blastnIdxV",
                    self.getParam("blastnIdxV"),
                    "--thread_no",
                    self.getParam("threads"),
                ]
            )
            all_cmd.append(tmp_cmd)

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            for i in self.getOutput("outdir"):
                if not os.path.exists(i):
                    os.mkdir(i)
            if verbose:
                self.run(all_cmd)
            else:
                self.multiRun(args=all_cmd, func=None, nCore=1)

            if virusIDfile:
                self.setInput("virusIDfile", virusIDfile)
                for i in self.getOutput("viruslistOutput"):
                    add_virus_note(self.getInput("virusIDfile"), i, True)

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)


def add_virus_note(virus_name, piko, head=False):
    """
    This function is for adding the virus note info into the viruslistOutput.
    piko means the viruslistOutput.
    """
    virus_dict = {"*": "*"}
    with open(virus_name) as INPUT:
        for line in INPUT:
            tmp = re.split(r"\S+", line.strip(), 1)
            try:
                virus_dict.setdefault(tmp[0], tmp[1])
            except IndexError:
                virus_dict.setdefault(tmp[0], tmp[0])

    with open(piko) as INPUT:
        content = INPUT.readlines()

    with open(piko, "w") as OUTPUT:
        if head:
            head = content.pop(0).split("\t")
            head.insert(1, "Annotation")
            OUTPUT.write("\t".join(head))
        else:
            OUTPUT.write(
                "virusID\tAnnotation\tvirus length\tmapping reads\tunmapping reads\n"
            )

        for line in content:
            tmp = line.strip().split("\t")
            note_info = virus_dict.get(tmp[0], tmp[0])
            tmp.insert(1, note_info)
            OUTPUT.write("\t".join(tmp) + "\n")
