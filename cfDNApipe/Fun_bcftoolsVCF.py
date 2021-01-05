# -*- coding: utf-8 -*-
"""
Created on Wed May 27 11:27:32 2020
@author: Shuying He
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, maxCore
import os
from .Configure import Configure
import math

__metaclass__ = type


class bcftoolsVCF(StepBase):
    def __init__(
        self,
        vcfInput=None,  # list ['samp1-1.vcf', 'samp1-2.vcf's...]
        outputdir=None,
        other_params={"-f": "'PASS'"},
        suffix="somatic",
        stepNum=None,
        upstream=None,
        threads=1,
        verbose=False,
        **kwargs
    ):
        '''
       This funtion is for picking out designate variants from vcf by bcftools view.
       You can set the selected params in other_params and set out suffix name.
       Default setting is for selecting the somatic mutation.

        bcftoolsVCF(vcfInput=None, outputdir=None,
            other_params={"-f": "'PASS'"},
            suffix="somatic",
            stepNum=None,
            upstream=None,
            threads=1,
            verbose=False,)

        {P}arameters:
        -  vcfInput: list, input vcf files.
        -  outputdir: str, output result folder, None means the same folder as input files.
        -  suffix: str, file name suffix.
        -  stepNum: int or str, step flag for folder name.
        -  threads: int, how many thread to use.
        -  upstream: upstream output results, used for pipeline, just can be addRG. This parameter can be True, which means a new pipeline start.
        -  verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        -  other_params: filter germline using {"-f": "'germline'"}, filter somatic using {"-f": "'PASS'"}.

        '''
        super(bcftoolsVCF, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
            # In this situation, input file and output path should be checked
            self.setInput("vcfInput", vcfInput)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ in ["gatherVCF", "filterMutectCalls"]:
                self.setInput("vcfInput", upstream.getOutput("vcfOutput"))
            else:
                raise commonError(
                    "Parameter upstream must from gatherVCF / filterMutectCalls."
                )
        self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            self.setParam("threads", threads)
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("vcfInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())

        self.setOutput(
            "vcfOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "."
                + suffix
                + ".vcf.gz"
                for x in self.getInput("vcfInput")
            ],
        )

        if other_params is not None:
            self.setParam("other_params", other_params)
        else:
            self.setParam("other_params", "")

        all_cmd = []

        for i in range(len(self.getInput("vcfInput"))):
            tmp_cmd = self.cmdCreate(
                [
                    "bcftools",
                    "view",
                    self.getParam("other_params"),
                    self.getInput("vcfInput")[i],
                    "-o",
                    self.getOutput("vcfOutput")[i],
                    "-O",
                    "z",
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
                    nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
                )

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)
