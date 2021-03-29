# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 10:51:10 2019

@author: Jiaqi Huang

"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, fraglendistribution, fraglencompplot, maxCore
import os
import math
from .Configure2 import Configure2

__metaclass__ = type


class fraglenplot_comp(StepBase):
    def __init__(
        self,
        casebedInput=None,
        ctrlbedInput=None,
        outputdir=None,
        maxLimit=500,
        ratio1=150,
        ratio2=400,
        labelInput=None,
        threads=1,
        stepNum=None,
        caseupstream=None,
        ctrlupstream=None,
        verbose=True,
        **kwargs
    ):
        """
        This function is used for compare the fragments' lengths in case and control samples.

        fraglenplot_comp(casebedInput=None, ctrlbedInput=None, outputdir=None, maxLimit=500, ratio1=150, ratio2=400, labelInput=None,
                         threads=1, stepNum=None, caseupstream=None, ctrlupstream=None, verbose=True,)
        {P}arameters:
            casebedInput: list, input bed files of case samples.
            ctrlbedInput: list, input bed files of control samples.
            outputdir: str, output result folder, None means the same folder as input files.
            maxLimit: int, maximum length to be considered.
            ratio1, ratio2: proportion statistics break point, default: 150, 400
            labelInput: list, [name_of_case, name_of_control](e.g. ["HCC", "CTR"])
            threads: int, how many thread to use.
            stepNum: int or str, step flag for folder name.
            caseupstream: upstream output results, used for pipeline.
            ctrlupstream: upstream output results, used for pipeline.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """

        if (stepNum is None) and (caseupstream is not None) and (ctrlupstream is None):
            super(fraglenplot_comp, self).__init__(stepNum, caseupstream)
        elif (
            (stepNum is None) and (caseupstream is None) and (ctrlupstream is not None)
        ):
            super(fraglenplot_comp, self).__init__(stepNum, ctrlupstream)
        elif (
            (stepNum is None)
            and (caseupstream is not None)
            and (ctrlupstream is not None)
        ):
            if caseupstream.getStepID() >= ctrlupstream.getStepID():
                super(fraglenplot_comp, self).__init__(stepNum, caseupstream)
            else:
                super(fraglenplot_comp, self).__init__(stepNum, ctrlupstream)
        else:
            super(fraglenplot_comp, self).__init__(stepNum)

        labelflag = False

        # set casetxtInput and ctrltxtInput
        if (
            ((caseupstream is None) and (ctrlupstream is None))
            or (caseupstream is True)
            or (ctrlupstream is True)
        ):
            self.setInput("casebedInput", casebedInput)
            self.setInput("ctrlbedInput", ctrlbedInput)
        else:
            Configure2.configureCheck()
            caseupstream.checkFilePath()
            ctrlupstream.checkFilePath()
            if caseupstream.__class__.__name__ == "bam2bed":
                self.setInput("casebedInput", caseupstream.getOutput("bedOutput"))
            else:
                raise commonError("Parameter upstream must from bam2bed.")

            if ctrlupstream.__class__.__name__ == "bam2bed":
                self.setInput("ctrlbedInput", ctrlupstream.getOutput("bedOutput"))
            else:
                raise commonError("Parameter upstream must from bam2bed.")

        self.checkInputFilePath()

        # set outputdir
        if (caseupstream is None) and (ctrlupstream is None):
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("casebedInput")[1])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set labelInput
        if labelInput is not None:
            self.setParam("label", labelInput)
            labelflag = True

        # set threads
        if (caseupstream is None) and (ctrlupstream is None):
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure2.getThreads())

        # set maxLimit
        self.setParam("maxLimit", maxLimit)
        
        self.setParam("ratio", [ratio1, ratio2])

        self.setOutput(
            "caseplotOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "_fraglen.png"
                for x in self.getInput("casebedInput")
            ],
        )

        self.setOutput(
            "ctrlplotOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "_fraglen.png"
                for x in self.getInput("ctrlbedInput")
            ],
        )

        self.setOutput(
            "plotOutput",
            [
                self.getOutput("outputdir") + "/" + "length_distribution.png",
                self.getOutput("outputdir") + "/" + "propotion.png",
            ],
        )

        self.setOutput(
            "casepickleOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "_fraglen.pickle"
                for x in self.getInput("casebedInput")
            ],
        )

        self.setOutput(
            "ctrlpickleOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "_fraglen.pickle"
                for x in self.getInput("ctrlbedInput")
            ],
        )
        
        self.setOutput(
            "txtOutput",
            os.path.join(self.getOutput("outputdir"), "statistic.txt"),
        )

        finishFlag = self.stepInit(caseupstream)

        if not finishFlag:
            case_multi_run_len = len(self.getInput("casebedInput"))
            ctrl_multi_run_len = len(self.getInput("ctrlbedInput"))
            if verbose:
                for i in range(case_multi_run_len):
                    print(
                        "Now, ploting fragment length distribution for "
                        + self.getInput("casebedInput")[i]
                    )
                    fraglendistribution(
                        bedInput=self.getInput("casebedInput")[i],
                        plotOutput=self.getOutput("caseplotOutput")[i],
                        pickleOutput=self.getOutput("casepickleOutput")[i],
                        maxLimit=self.getParam("maxLimit"),
                    )
                for i in range(ctrl_multi_run_len):
                    print(
                        "Now, ploting fragment length distribution for "
                        + self.getInput("ctrlbedInput")[i]
                    )
                    fraglendistribution(
                        bedInput=self.getInput("ctrlbedInput")[i],
                        plotOutput=self.getOutput("ctrlplotOutput")[i],
                        pickleOutput=self.getOutput("ctrlpickleOutput")[i],
                        maxLimit=self.getParam("maxLimit"),
                    )
                if labelflag:
                    fraglencompplot(
                        caseInput=self.getOutput("casepickleOutput"),
                        ctrlInput=self.getOutput("ctrlpickleOutput"),
                        plotOutput=self.getOutput("plotOutput"),
                        txtOutput=self.getOutput("txtOutput"),
                        ratio=self.getParam("ratio"),
                        labelInput=self.getParam("label"),
                    )
                else:
                    fraglencompplot(
                        caseInput=self.getOutput("casepickleOutput"),
                        ctrlInput=self.getOutput("ctrlpickleOutput"),
                        plotOutput=self.getOutput("plotOutput"),
                        txtOutput=self.getOutput("txtOutput"),
                        ratio=self.getParam("ratio"),
                    )
            else:
                case_args = [
                    [
                        self.getInput("casebedInput")[i],
                        self.getOutput("caseplotOutput")[i],
                        self.getOutput("casepickleOutput")[i],
                        self.getParam("maxLimit"),
                    ]
                    for i in range(case_multi_run_len)
                ]
                self.multiRun(
                    args=case_args,
                    func=fraglendistribution,
                    nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
                )
                ctrl_args = [
                    [
                        self.getInput("ctrlbedInput")[i],
                        self.getOutput("ctrlplotOutput")[i],
                        self.getOutput("ctrlpickleOutput")[i],
                        self.getParam("maxLimit"),
                    ]
                    for i in range(ctrl_multi_run_len)
                ]
                self.multiRun(
                    args=ctrl_args,
                    func=fraglendistribution,
                    nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
                )
                if labelflag:
                    fraglencompplot(
                        caseInput=self.getOutput("casepickleOutput"),
                        ctrlInput=self.getOutput("ctrlpickleOutput"),
                        plotOutput=self.getOutput("plotOutput"),
                        txtOutput=self.getOutput("txtOutput"),
                        ratio=self.getParam("ratio"),
                        labelInput=self.getParam("label"),
                    )
                else:
                    fraglencompplot(
                        caseInput=self.getOutput("casepickleOutput"),
                        ctrlInput=self.getOutput("ctrlpickleOutput"),
                        plotOutput=self.getOutput("plotOutput"),
                        txtOutput=self.getOutput("txtOutput"),
                        ratio=self.getParam("ratio"),
                    )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
