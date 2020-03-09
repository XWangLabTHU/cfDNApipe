# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 10:51:10 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""


from .StepBase2 import StepBase2
from .cfDNA_utils import commonError, fraglendistribution, fraglencompplot
import os
from .Configure2 import Configure2


__metaclass__ = type


class fraglenplot_comp(StepBase2):
    def __init__(
        self,
        casebedInput=None,  # list
        ctrlbedInput=None,  # list
        outputdir=None,  # str
        maxLimit=500,
        stepNum=None,
        caseupstream=None,
        ctrlupstream=None,
        **kwargs
    ):
        if (stepNum is None) and (caseupstream is not None) and (ctrlupstream
                                                                 is None):
            super(fraglenplot_comp, self).__init__(stepNum, caseupstream)
        elif ((stepNum is None) and (caseupstream is None)
              and (ctrlupstream is not None)):
            super(fraglenplot_comp, self).__init__(stepNum, ctrlupstream)
        elif ((stepNum is None) and (caseupstream is not None)
              and (ctrlupstream is not None)):
            if caseupstream.getStepID() >= ctrlupstream.getStepID():
                super(fraglenplot_comp, self).__init__(stepNum, caseupstream)
            else:
                super(fraglenplot_comp, self).__init__(stepNum, ctrlupstream)
        else:
            super(fraglenplot_comp, self).__init__(stepNum)
            
        if caseupstream is None and ctrlupstream is None:
            self.setInput("casebedInput", casebedInput)
            self.setInput("ctrlbedInput", ctrlbedInput)
            self.checkInputFilePath()

            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("casebedInput")[1])),
                )
            else:
                self.setOutput("outputdir", outputdir)

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

            self.setOutput("outputdir", self.getStepFolderPath())

        self.setParam("maxLimit", maxLimit)
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
            "compplotOutput",
            [
                self.getOutput("outputdir") + "/" + "length_distribution.png",
                self.getOutput("outputdir") + "/" + "propotion.png",
            ]
        )
        self.setOutput(
            "casenpyOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "_fraglen.npy"
                for x in self.getInput("casebedInput")
            ],
        )
        self.setOutput(
            "ctrlnpyOutput",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                + "_fraglen.npy"
                for x in self.getInput("ctrlbedInput")
            ],
        )

        finishFlag = self.stepInit(caseupstream)

        if finishFlag:
            self.excute(finishFlag)
        else:
            case_multi_run_len = len(self.getInput("casebedInput"))
            ctrl_multi_run_len = len(self.getInput("ctrlbedInput"))
            case_len_data = []
            ctrl_len_data = []
            for i in range(case_multi_run_len):
                print(
                    "Now, ploting fragment length distribution for "
                    + self.getInput("casebedInput")[i]
                )
                case_len_data.append(
                    fraglendistribution(
                        bedInput=self.getInput("casebedInput")[i],
                        plotOutput=self.getOutput("caseplotOutput")[i],
                        binOutput=self.getOutput("casenpyOutput")[i],
                        maxLimit=self.getParam("maxLimit"),
                    )
                )
            for i in range(ctrl_multi_run_len):
                print(
                    "Now, ploting fragment length distribution for "
                    + self.getInput("ctrlbedInput")[i]
                )
                ctrl_len_data.append(
                    fraglendistribution(
                        bedInput=self.getInput("ctrlbedInput")[i],
                        plotOutput=self.getOutput("ctrlplotOutput")[i],
                        binOutput=self.getOutput("ctrlnpyOutput")[i],
                        maxLimit=self.getParam("maxLimit"),
                    )
                )

            fraglencompplot(
                caseInput = case_len_data,
                ctrlInput = ctrl_len_data,
                plotOutput = self.getOutput("compplotOutput"),
            )

            self.excute(finishFlag, runFlag=False)
