# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
@author: LY, Shuying He
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, flatten
import os
from .Configure import Configure

__metaclass__ = type


class seqtk(StepBase):
    """for extracting the virus fasta sequence."""

    def __init__(
        self,
        virusID=None,
        virusDB=None,
        OutputDir=None,
        stepNum=None,
        upstream=None,
        **kwargs
    ):

        super(seqtk, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
            self.setInput("virusDB", virusDB)
        else:
            # check Configure for running pipeline
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ in ["virusdetect", "virusID_extract"]:
                self.setInput("virusDB", upstream.getOutput("virusDB"))
            else:
                raise commonError(
                    "Parameter upstream must from virusdetect or virusID_extract."
                )
        self.checkInputFilePath()

        if upstream is None:
            if not OutputDir:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("virusDB"))),
                )
            else:
                self.setOutput("outputdir", OutputDir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set virusID
        if virusID:
            self.setParam("virusID", treat_virusID(virusID))
        elif upstream.__class__.__name__ in ["virusdetect", "virusID_extract"]:
            virusID = [
                grep_top3_virus(x) for x in upstream.getOutput("viruslistOutput")
            ]
            virusID = set(flatten(virusID))
            self.setParam("virusID", virusID)
        else:
            raise commonError("Please give me virusID.")

        self.setOutput(
            "txtOutput", os.path.join(self.getOutput("outputdir"), "virus-id.txt")
        )
        self.setOutput(
            "faOutput", os.path.join(self.getOutput("outputdir"), "virus.fa")
        )

        cmd = self.cmdCreate(
            [
                "seqtk",
                "subseq",
                self.getInput("virusDB"),
                self.getOutput("txtOutput"),
                ">",
                self.getOutput("faOutput"),
            ]
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            with open(self.getOutput("txtOutput"), "w") as f:
                f.write("\n".join(self.getParam("virusID")) + "\n")
            self.run(cmd)

        self.stepInfoRec(cmds=[cmd], finishFlag=finishFlag)


def treat_virusID(virusID):
    "virus_ID should convert to list"
    if isinstance(virusID, str):
        virusID = [virusID]
    return virusID


def grep_top3_virus(viruslistFile):
    """grep the first column virus id from virus-list.txt."""
    with open(viruslistFile) as INPUT:
        tmp = INPUT.readlines()
        top3_line = tmp[1:3]
    viruslist = list(map(lambda x: x.strip().split("\t")[0], top3_line))
    return viruslist
