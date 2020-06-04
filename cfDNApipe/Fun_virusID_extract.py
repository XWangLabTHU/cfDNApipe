# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
@author: Shuying He
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import re


__metaclass__ = type


class virusID_extract(StepBase):
    def __init__(
        self,
        bamInput=None,
        virusDB=None,
        virusIDfile=None,
        threads=1,
        OutputDir=None,
        stepNum=None,
        upstream=None,
        **kwargs
    ):
        super(virusID_extract, self).__init__(stepNum, upstream)

        if (upstream is None) or (upstream is True):
            self.setInput("bamInput", bamInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "mapQ_filter":
                self.setInput("bamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from mapQ_filter.")
        self.checkInputFilePath()

        if upstream is None:
            if OutputDir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
                )
            else:
                self.setOutput("outputdir", OutputDir)
            self.setParam("threads", threads)

        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())

        self.setOutput("virusDB", virusDB)

        self.setOutput(
            "viruslistOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + ".stat.txt"
                for x in self.getInput("bamInput")
            ],
        )

        multi_run_len = len(self.getInput("bamInput"))

        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(
                [
                    "samtools",
                    "idxstats",
                    "-@",
                    self.getParam("threads"),
                    self.getInput("bamInput")[i],
                    "|",
                    "sort -k 3 -rn - >",
                    self.getOutput("viruslistOutput")[i],
                ]
            )
            all_cmd.append(tmp_cmd)

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            self.run(all_cmd)
            if virusIDfile:
                self.setInput("virusIDfile", virusIDfile)
                for i in self.getOutput("viruslistOutput"):
                    add_virus_note(
                        self.getInput("virusIDfile"), i,
                    )

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
