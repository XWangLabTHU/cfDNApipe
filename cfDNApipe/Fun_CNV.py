# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:51:54 2019

@author: zhang, Huang
"""

from .StepBase2 import StepBase2
from .cfDNA_utils import commonError, wig2df, correctReadCount, plotCNVscatter, sumChromarm, plotCNVheatmap
import pandas as pd
import numpy as np
import os
from .Configure2 import Configure2

__metaclass__ = type


class computeCNV(StepBase2):
    def __init__(
            self,
            fastaInput=None,
            casebamInput=None,  # list
            ctrlbamInput=None,  # list
            outputdir=None,  # str
            caseupstream=None,
            ctrlupstream=None,
            cytoBandInput=None,
            stepNum=None,
            **kwargs):
        if (stepNum is None) and (caseupstream is not None) and (ctrlupstream is None):
            super(computeCNV, self).__init__(stepNum, caseupstream)
        elif ((stepNum is None) and (caseupstream is None) and (ctrlupstream is not None)):
            super(computeCNV, self).__init__(stepNum, ctrlupstream)
        elif ((stepNum is None) and (caseupstream is not None) and (ctrlupstream is not None)):
            if caseupstream.getStepID() >= ctrlupstream.getStepID():
                super(computeCNV, self).__init__(stepNum, caseupstream)
            else:
                super(computeCNV, self).__init__(stepNum, ctrlupstream)
        else:
            super(computeCNV, self).__init__(stepNum)

        if caseupstream is None and ctrlupstream is None:
            self.setInput("casebamInput", casebamInput)
            self.setInput("ctrlbamInput", ctrlbamInput)
            self.checkInputFilePath()

            if fastaInput is None:
                self.setInput("fastaInput", Configure2.getConfig("genome.seq"))
            else:
                self.setInput("fastaInput", fastaInput)

            if outputdir is None:
                self.setOutput("outputdir", os.path.dirname(
                    os.path.abspath(self.getInput("casebamInput")[1])))
            else:
                self.setOutput("outputdir", outputdir)

        else:
            Configure2.configureCheck()
            caseupstream.checkFilePath()
            ctrlupstream.checkFilePath()

            if caseupstream.__class__.__name__ == "bamsort":
                self.setInput("casebamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from bamsort.")
            if ctrlupstream.__class__.__name__ == "bamsort":
                self.setInput("ctrlbamInput", upstream.getOutput("bamOutput"))
            else:
                raise commonError("Parameter upstream must from bamsort.")

            if fastaInput is None:
                self.setInput("fastaInput", Configure2.getConfig("genome.seq"))
            else:
                self.setInput("fastaInput", fastaInput)

            self.checkInputFilePath()

            self.setOutput("outputdir", self.getStepFolderPath())
            
        if cytoBandInput is None:
            self.setInput("cytoBandInput", Configure2.getConfig("cytoBand"))
        else:
            self.setInput("cytoBandInput", cytoBandInput)
            
        self.setOutput("gcOutput", os.path.join(self.getOutput(
            "outputdir"), self.getMaxFileNamePrefixV2(self.getInput("fastaInput"))) + ".gc.wig")
            
        self.setOutput("casereadOutput", [os.path.join(self.getOutput(
            "outputdir"), self.getMaxFileNamePrefixV2(x)) + ".read.wig" for x in self.getInput("casebamInput")])   
            
        self.setOutput("ctrlreadOutput", [os.path.join(self.getOutput(
            "outputdir"), self.getMaxFileNamePrefixV2(x)) + ".read.wig" for x in self.getInput("ctrlbamInput")])
            
        self.setOutput(
            "txtOutput",
            self.getOutput("outputdir") + "/Z-score.txt")

        self.setOutput(
            "casereadplotOutput",
            [self.getOutput("outputdir") + "/" + self.getMaxFileNamePrefixV2(x) +
             "-gc_cor.png" for x in self.getInput("casebamInput")]
        )

        self.setOutput(
            "ctrlreadplotOutput",
            [self.getOutput("outputdir") + "/" + self.getMaxFileNamePrefixV2(x) +
             "-gc_cor.png" for x in self.getInput("ctrlbamInput")]
        )

        self.setOutput(
            "casecnvplotOutput",
            [self.getOutput("outputdir") + "/" + self.getMaxFileNamePrefixV2(x) +
             "-Z_score.png" for x in self.getInput("casebamInput")]
        )

        self.setOutput(
            "ctrlcnvplotOutput",
            [self.getOutput("outputdir") + "/" + self.getMaxFileNamePrefixV2(x) +
             "-Z_score.png" for x in self.getInput("ctrlbamInput")]
        )

        self.setOutput(
            "plotOutput",
            self.getOutput("outputdir") + "/CNV.png")

        finishFlag = self.stepInit(caseupstream)
        
        #add cmds into all_cmd
        case_multi_run_len = len(self.getInput("casebamInput"))
        ctrl_multi_run_len = len(self.getInput("ctrlbamInput"))
        all_cmd = []
        if not os.path.exists(self.getOutput("gcOutput")):
            gc_tmp_cmd = self.cmdCreate(["gcCounter",
                                         "-w", 100000,
                                         self.getInput("fastaInput"),
                                         ">", self.getOutput("gcOutput")])
            all_cmd.append(gc_tmp_cmd)
        for i in range(case_multi_run_len):
            if not os.path.exists(self.getOutput("casereadOutput")[i]):
                case_read_tmp_cmd = self.cmdCreate(["readCounter",
                                               "-w", 100000,
                                               self.getInput("casebamInput")[i],
                                               ">", self.getOutput("casereadOutput")[i]])
                all_cmd.append(case_read_tmp_cmd)
        for i in range(ctrl_multi_run_len):
            if not os.path.exists(self.getOutput("ctrlreadOutput")[i]):
                ctrl_read_tmp_cmd = self.cmdCreate(["readCounter",
                                               "-w", 100000,
                                               self.getInput("ctrlbamInput")[i],
                                               ">", self.getOutput("ctrlreadOutput")[i]])
                all_cmd.append(ctrl_read_tmp_cmd)
        
        if not finishFlag:
            #run commands
            self.run(all_cmd)
            
            #process CNV
            case_binpos = [[[], []] for i in range(case_multi_run_len)]
            ctrl_binpos = [[[], []] for i in range(ctrl_multi_run_len)]
            case_bin = [[] for i in range(case_multi_run_len)]
            ctrl_bin = [[] for i in range(ctrl_multi_run_len)]
            case_chrom = [[] for i in range(case_multi_run_len)]
            ctrl_chrom = [[] for i in range(ctrl_multi_run_len)]
            case_tmp_maxlen = -1
            ctrl_tmp_maxlen = -1
            genes = [] #stores the names of the chromosome arms
            #case
            case_df_gc = wig2df(self.getOutput("gcOutput"))
            for i in range(case_multi_run_len):
                print("Now, processing", self.getMaxFileNamePrefixV2(
                    self.getOutput("casereadOutput")[i]), "...")
                case_df_read = wig2df(self.getOutput("casereadOutput")[i])
                
                #GC correction
                case_read_correct, case_read_correct2 = correctReadCount(
                    case_df_read, case_df_gc, self.getOutput("casereadplotOutput")[i])
                case_binpos[i][0] = case_read_correct2["chrom"].tolist()
                case_binpos[i][1] = [int(x.split("-")[0]) for x in case_read_correct2["start-end"].tolist()]
                case_bin[i] = case_read_correct2["value"].tolist()
                case_tmp_maxlen = max(case_tmp_maxlen, len(case_bin[i]))
                
                #sum the read count data into each chromosome arm
                case_chrom[i], genes = sumChromarm(
                    case_read_correct, self.getInput("cytoBandInput"))
            '''        
            #add Nonetype values to the tail of each row in case_bin to fullfill it
            for i in range(case_multi_run_len):
                for j in range(case_tmp_maxlen - len(case_bin[i])):
                    case_bin[i].append(None)
            '''        
            #build the case bin-divided dataframe, columns are case samples, rows are bins
            case_bin_df = pd.DataFrame(
                np.transpose(case_bin),
                columns=[self.getMaxFileNamePrefixV2(x).split(
                    ".")[0] for x in self.getOutput("casereadOutput")],
                index=None
            )
            #build the case chromarms-divided dataframe, columns are case samples, rows are chromosome arms
            case_chrom_df = pd.DataFrame(
                np.transpose(case_chrom),
                columns=[self.getMaxFileNamePrefixV2(x).split(
                    ".")[0] for x in self.getOutput("casereadOutput")],
                index=genes
            )

            #control
            ctrl_df_gc = wig2df(self.getOutput("gcOutput"))
            for i in range(ctrl_multi_run_len):
                print("Now, processing", self.getMaxFileNamePrefixV2(
                    self.getOutput("ctrlreadOutput")[i]), "...")
                ctrl_df_read = wig2df(self.getOutput("ctrlreadOutput")[i])
                ctrl_read_correct, ctrl_read_correct2 = correctReadCount(
                    ctrl_df_read, ctrl_df_gc, self.getOutput("ctrlreadplotOutput")[i])
                ctrl_binpos[i][0] = ctrl_read_correct2["chrom"].tolist()
                ctrl_binpos[i][1] = [int(x.split("-")[0]) for x in ctrl_read_correct2["start-end"].tolist()]
                ctrl_bin[i] = ctrl_read_correct2["value"].tolist()
                ctrl_tmp_maxlen = max(ctrl_tmp_maxlen, len(ctrl_bin[i]))
                ctrl_chrom[i], genes = sumChromarm(
                    ctrl_read_correct, self.getInput("cytoBandInput"))
            for i in range(ctrl_multi_run_len):
                for j in range(ctrl_tmp_maxlen - len(ctrl_bin[i])):
                    ctrl_bin[i].append(None)
            ctrl_bin_df = pd.DataFrame(
                np.transpose(ctrl_bin),
                columns=[self.getMaxFileNamePrefixV2(x).split(
                    ".")[0] for x in self.getOutput("ctrlreadOutput")],
                index=None
            )
            ctrl_chrom_df = pd.DataFrame(
                np.transpose(ctrl_chrom),
                columns=[self.getMaxFileNamePrefixV2(x).split(
                    ".")[0] for x in self.getOutput("ctrlreadOutput")],
                index=genes
            )
            
            #compute z-score and plot scatter-plot
            plotCNVscatter(case_bin_df, ctrl_bin_df, case_binpos, ctrl_binpos, 
                self.getInput("cytoBandInput"), self.getOutput(
                "casecnvplotOutput"), self.getOutput("ctrlcnvplotOutput"))
            #compute z-score and plot heatmap
            plotCNVheatmap(case_chrom_df, ctrl_chrom_df, self.getOutput(
                "txtOutput"), self.getOutput("plotOutput"))

        self.stepInfoRec(cmds=[all_cmd], finishFlag=finishFlag)
