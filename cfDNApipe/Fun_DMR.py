# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:51:54 2019

@author: zhang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, processDMR
import pandas as pd
import os
from .Configure2 import Configure2

__metaclass__ = type


class computeDMR(StepBase):
    def __init__(
        self,
        casetxtInput=None,
        ctrltxtInput=None,
        outputdir=None,
        threads=1,
        adjmethod=None,
        caseupstream=None,
        ctrlupstream=None,
        stepNum=None,
        **kwargs
    ):
        """
        This function is used for compute DMR between case and aontrol samples.

        computeDMR(casetxtInput=None, ctrltxtInput=None, outputdir=None, threads=1, adjmethod=None, caseupstream=None, ctrlupstream=None, stepNum=None)
        {P}arameters:
            casetxtInput: list, input methylation level files of case samples.
            ctrlbedInput: list, input methylation level files of control samples.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            adjmethod: str, method of p_value correction, must be "bonferroni", "fdr_bh"(default), "fdr_by" or "holm"
            caseupstream: upstream output results, used for pipeline.
            ctrlupstream: upstream output results, used for pipeline.
            stepNum: int or str, step flag for folder name.
        """

        if (stepNum is None) and (caseupstream is not None) and (ctrlupstream is None):
            super(computeDMR, self).__init__(stepNum, caseupstream)
        elif (
            (stepNum is None) and (caseupstream is None) and (ctrlupstream is not None)
        ):
            super(computeDMR, self).__init__(stepNum, ctrlupstream)
        elif (
            (stepNum is None)
            and (caseupstream is not None)
            and (ctrlupstream is not None)
        ):
            if caseupstream.getStepID() >= ctrlupstream.getStepID():
                super(computeDMR, self).__init__(stepNum, caseupstream)
            else:
                super(computeDMR, self).__init__(stepNum, ctrlupstream)
        else:
            super(computeDMR, self).__init__(stepNum)

        # set casetxtInput and ctrltxtInput
        if (
            ((caseupstream is None) and (ctrlupstream is None))
            or (caseupstream is True)
            or (ctrlupstream is True)
        ):
            self.setInput("casetxtInput", casetxtInput)
            self.setInput("ctrltxtInput", ctrltxtInput)
        else:
            Configure2.configureCheck()
            caseupstream.checkFilePath()
            ctrlupstream.checkFilePath()
            if caseupstream.__class__.__name__ == "calculate_methyl":
                self.setInput("casetxtInput", caseupstream.getOutput("txtOutput"))
            else:
                raise commonError("Parameter caseupstream must from calculate_methyl.")
            if ctrlupstream.__class__.__name__ == "calculate_methyl":
                self.setInput("ctrltxtInput", ctrlupstream.getOutput("txtOutput"))
            else:
                raise commonError("Parameter ctrlupstream must from calculate_methyl.")

        self.checkInputFilePath()

        # set outputdir
        if (caseupstream is None) and (ctrlupstream is None):
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("casetxtInput")[1])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set threads
        if (caseupstream is None) and (ctrlupstream is None):
            self.setParam("threads", threads)
        else:
            self.setParam("threads", Configure2.getThreads())

        # set adjmethod
        if adjmethod is not None:
            self.setParam("adjmethod", adjmethod)
        else:
            self.setParam("adjmethod", "fdr_bh")

        self.setOutput(
            "casetxtOutput", os.path.join(self.getOutput("outputdir"), "case_DMR.txt"),
        )

        self.setOutput(
            "ctrltxtOutput",
            os.path.join(self.getOutput("outputdir"), "control_DMR.txt"),
        )

        self.setOutput(
            "txtOutput", os.path.join(self.getOutput("outputdir"), "DMR.txt"),
        )

        finishFlag = self.stepInit(caseupstream)

        if not finishFlag:
            case_multi_run_len = len(self.getInput("casetxtInput"))
            ctrl_multi_run_len = len(self.getInput("ctrltxtInput"))
            for i in range(case_multi_run_len):
                data = pd.read_csv(
                    self.getInput("casetxtInput")[i],
                    sep="\t",
                    header=0,
                    names=["chr", "start", "end", "unmCpG", "mCpG", "mlCpG", ],
                )
                if i == 0:
                    ml_df = pd.DataFrame(
                        {
                            "chr": data["chr"],
                            "start": data["start"],
                            "end": data["end"],
                            os.path.split(self.getInput("casetxtInput")[i])[1]: data[
                                "mlCpG"
                            ],
                        }
                    )
                else:
                    ml_df = pd.concat(
                        [
                            ml_df,
                            pd.DataFrame(
                                {
                                    os.path.split(self.getInput("casetxtInput")[i])[
                                        1
                                    ]: data["mlCpG"]
                                }
                            ),
                        ],
                        axis=1,
                    )
            for i in range(ctrl_multi_run_len):
                data = pd.read_csv(
                    self.getInput("ctrltxtInput")[i],
                    sep="\t",
                    header=0,
                    names=["chr", "start", "end", "unmCpG", "mCpG", "mlCpG", ],
                )
                ml_df = pd.concat(
                    [
                        ml_df,
                        pd.DataFrame(
                            {
                                os.path.split(self.getInput("ctrltxtInput")[i])[
                                    1
                                ]: data["mlCpG"]
                            }
                        ),
                    ],
                    axis=1,
                )
            processDMR(
                ml_df,
                case_multi_run_len,
                self.getParam("adjmethod"),
                self.getOutput("casetxtOutput"),
                self.getOutput("ctrltxtOutput"),
                self.getOutput("txtOutput"),
            )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
