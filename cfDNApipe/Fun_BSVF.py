# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:27:32 2019
Modify on Sun Apr 26 11:27:32 2020
@author: LY, Shuying He
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import math
import pkg_resources

__metaclass__ = type


class BSVF(StepBase):
    def __init__(
        self,
        fqInput1=None,
        fqInput2=None,
        plInput=None,
        bsuit=None,
        hostRef=None,
        virusRef=None,
        OutputDir=None,
        aln="bwa",
        MinVirusLength="20",
        stepNum=None,
        upstream1=None,
        upstream2=None,
        threads=1,
        verbose=False,
        **kwargs,
    ):

        if upstream2:
            super(BSVF, self).__init__(stepNum, upstream2)
        else:
            super(BSVF, self).__init__(stepNum, upstream1)

        if (upstream1 is None) or (upstream1 is True):
            # In this situation, input file and output path should be checked
            self.setInput("fqInput1", fqInput1)
            self.setInput("fqInput2", fqInput2)
        else:
            Configure.configureCheck()
            if upstream1.__class__.__name__ == "adapterremoval":
                upstream1.checkFilePath()
                self.setInput("fqInput1", upstream1.getOutput("pair1"))
                self.setInput("fqInput2", upstream1.getOutput("pair2"))
            else:
                raise commonError("Parameter upstream1 must from adapterremoval.")
        self.checkInputFilePath()

        if upstream1 is None:
            self.setParam("threads", threads)
            if OutputDir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("fqInput1")[0])),
                )
            else:
                self.setOutput("outputdir", OutputDir)
        else:
            # check Configure for running pipeline
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())

        if virusRef:
            self.setInput("virusRef", virusRef)
        elif upstream2.__class__.__name__ == "seqtk":
            upstream2.checkFilePath()
            self.setInput("virusRef", upstream2.getOutput("faOutput"))
        else:
            raise commonError("Please set virusRef or set seqtk as upstream2.")

        self.setParam(
            "sampleName",
            [self.getMaxFileNamePrefixV2(x) for x in self.getInput("fqInput1")],
        )

        if bsuit:
            self.setInput("bsuit", bsuit)
        else:
            self.setInput(
                "bsuit", pkg_resources.resource_filename("cfDNApipe", "data/BSVF/bsuit")
            )

        if plInput:
            self.setInput("plInput", plInput)
        else:
            self.setInput(
                "plInput",
                pkg_resources.resource_filename(
                    "cfDNApipe", "data/BSVF/BSVF_prepare_configFile.pl"
                ),
            )

        self.setInput("hostRef", hostRef)

        self.setParam("aln", aln)
        self.setParam("MinVirusLength", MinVirusLength)

        self.setOutput(
            "breakpointFile",
            [
                os.path.join(self.getOutput("outputdir"), f"{x}_analyse/{x}.analyse")
                for x in self.getParam("sampleName")
            ],
        )

        index_cmd = []

        # prepare ini file
        for i in range(len(self.getParam("sampleName"))):
            tmp_cmd = self.cmdCreate(
                [
                    "perl",
                    self.getInput("plInput"),
                    self.getInput("hostRef"),
                    self.getInput("virusRef"),
                    self.getInput("fqInput1")[i],
                    self.getInput("fqInput2")[i],
                    self.getOutput("outputdir"),
                    self.getParam("sampleName")[i],
                    self.getParam("aln"),
                    self.getParam("MinVirusLength"),
                ]
            )
            index_cmd.append(tmp_cmd)

        # create index
        tmp_cmd = self.cmdCreate(
            [
                self.getInput("bsuit"),
                "prepare",
                self.getOutput("outputdir")
                + "/"
                + self.getParam("sampleName")[0]
                + ".ini",
            ]
        )
        index_cmd.append(tmp_cmd)

        tmp_cmd = self.cmdCreate(
            [
                "sh",
                self.getOutput("outputdir")
                + "/"
                + self.getParam("sampleName")[0]
                + "_index.sh",
            ]
        )
        index_cmd.append(tmp_cmd)

        # bsuit aln
        aln_cmd = []
        for i in range(len(self.getParam("sampleName"))):
            tmp_cmd = self.cmdCreate(
                [
                    self.getInput("bsuit"),
                    "aln",
                    self.getOutput("outputdir")
                    + "/"
                    + self.getParam("sampleName")[i]
                    + ".ini",
                    ";",
                    "sh",
                    self.getOutput("outputdir")
                    + "/"
                    + self.getParam("sampleName")[i]
                    + "_aln.sh",
                ]
            )
            aln_cmd.append(tmp_cmd)

        finishFlag = self.stepInit(upstream1)

        if not finishFlag:
            self.run(index_cmd)
            if verbose:
                self.run(aln_cmd)
            else:
                self.multiRun(
                    args=aln_cmd,
                    func=None,
                    nCore=math.ceil(self.getParam("threads") / 4),
                )

        self.stepInfoRec(cmds=index_cmd + aln_cmd, finishFlag=finishFlag)
