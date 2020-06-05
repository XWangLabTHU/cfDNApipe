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

__metaclass__ = type


class gatherVCF(StepBase):
    def __init__(
        self,
        vcfInput=None,  # list [['samp1-1.vcf', 'samp1-2.vcf'...],[]...]
        outputdir=None,
        stepNum=None,
        upstream=None,
        threads=1,
        verbose=False,
        **kwargs
    ):

        super(gatherVCF, self).__init__(stepNum, upstream)
        if upstream is None:
            self.setParam("threads", threads)
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("vcfInput")[0][0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())

        if (upstream is None) or (upstream is True):
            # In this situation, input file and output path should be checked
            self.setInput("vcfInput", vcfInput)
            self.setParam("prefix", [self.getMaxFileNamePrefixV2(vcfInput[0])])
            self.setInput(
                "%s_vcfInput" % self.getParam("prefix")[0], self.getInput("vcfInput")
            )
        else:
            # check Configure for running pipeline
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "filterMutectCalls" or "VCFpostprocess":
                self.setInput("vcfInput", upstream.getOutput("vcfOutput"))
                vcfnum = len(self.getInput("vcfInput"))
                prefix = []
                if vcfnum % 25 == 0:
                    for i in range(0, vcfnum, 25):
                        pf = self.getMaxFileNamePrefixV2(
                            self.getInput("vcfInput")[i]
                        ).replace("_chr1", "")
                        prefix.append(pf)
                        self.setInput(
                            "%s_vcfInput" % pf, self.getInput("vcfInput")[i : i + 25]
                        )
                    self.setParam("prefix", prefix)
                else:
                    raise commonError("VCFInputs length is not a multiple of 25.")
            else:
                raise commonError(
                    "Parameter upstream must from filterMutectCalls or VCFpostprocess."
                )
        self.checkInputFilePath()

        self.setOutput(
            "vcfOutput",
            [
                os.path.join(self.getOutput("outputdir"), x + ".vcf.gz")
                for x in self.getParam("prefix")
            ],
        )

        all_cmd = []

        for i in range(len(self.getParam("prefix"))):
            tmp_cmd = self.cmdCreate(
                [
                    "gatk",
                    "GatherVcfs",
                    "-I",
                    " -I ".join(
                        self.getInput("%s_vcfInput" % self.getParam("prefix")[i])
                    ),
                    "-O",
                    self.getOutput("vcfOutput")[i],
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
                    nCore=math.ceil(self.getParam("threads") / 4),
                )

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)
