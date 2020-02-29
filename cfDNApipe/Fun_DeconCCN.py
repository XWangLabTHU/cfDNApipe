# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:17:42 2020

@author: Jiaqi Huang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, DeconCCN, DeconCCNplot
from .Configure import Configure
import numpy as np
import os

__metaclass__ = type


class runDeconCCN(StepBase):
    def __init__(
        self,
        mixInput=None,
        refInput=None,
        celltypeInput=None, #list
        outputdir=None,  # str
        threads=1,
        stepNum=None,
        upstream=None,
        **kwargs
    ):
        super(runDeconCCN, self).__init__(stepNum, upstream)

        if upstream is None:
            self.setInput("mixInput", mixInput)
            if refInput is None:
                self.setInput("refInput", Configure.getConfig("methylref")) #mark
            else:
                self.setInput("refInput", refInput)
            self.checkInputFilePath()

            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("mixInput"))),
                )
            else:
                self.setOutput("outputdir", outputdir)

            self.setParam("threads", threads)

        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "calculate_methyl":
                self.setInput("mixInput", upstream.getOutput("txtOutput"))
            else:
                raise commonError("Parameter upstream must from calculate_methyl.")

            if refInput is None:
                self.setInput("refInput", Configure.getConfig("methylref")) #mark
            else:
                self.setInput("refInput", refInput)
                
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())

        self.setInput("celltypeInput", celltypeInput)
        
        self.setOutput("txtOutput", self.getOutput("outputdir") + "/mixture.txt")
        self.setOutput("plotOutput", self.getOutput("outputdir") + "/bar-chart.png")

        finishFlag = self.stepInit(upstream)

        if finishFlag:
            self.excute(finishFlag)
        else:
            multi_run_len = len(self.getInput("mixInput"))
            mix = [[] for i in range(multi_run_len)]
            ref = np.load(self.getInput("refInput"))
            for i in range(multi_run_len):
                data = pd.read_csv(
                    self.getInput("mixInput")[i],
                    sep="\t",
                    header=0,
                    names=[
                        "chr",
                        "start",
                        "end",
                        "unmCpG",
                        "mCpG",
                        "mlCpG",
                    ],
                )
                mix[i] += DeconCCN(ref, data["mlCpG"].tolist())
            mix_df = pd.DataFrame(
                np.transpose(mix), 
                index = self.getInput("celltypeInput"),
                columns = [self.getMaxFileNamePrefixV2(x).split('.')[0] for x in self.getInput("mixInput")]
            )
            print(mix_df)
            mix_df.to_csv(self.getOutput("txtOutput"), sep = "\t", index = True)
            DeconCCNplot(mix_df, self.getOutput("plotOutput"))
            
            self.excute(finishFlag, runFlag=False)
