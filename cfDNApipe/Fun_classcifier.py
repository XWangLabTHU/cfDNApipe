# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:51:54 2019

@author: zhang, Huang
"""

from .StepBase2 import StepBase2
from .cfDNA_utils import commonError, get_cross_validation
import pandas as pd
import os
import numpy as np
from .Configure2 import Configure2
from sklearn.svm import SVC
import itertools
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import joblib

__metaclass__ = type


class classifier(StepBase2):
    def __init__(
        self,
        caseInput=None,
        ctrlInput=None,
        model=None,
        LOOCV=False,
        test_size=0.3,
        outputdir=None,
        threads=1,
        caseupstream=None,
        ctrlupstream=None,
        stepNum=None,
        **kwargs
    ):
        """
        This function is used for building classifier for analysis cfDNA data.

        classifier(caseInput=None, ctrlInput=None, model=None, LOOCV=FALSE, test_size=0.3,
                   outputdir=None, threads=1, caseupstream=None, ctrlupstream=None, stepNum=None,)
        {P}arameters:
            caseInput: str, integrated file for case samples. must be a text file with
                       rownames(which indicates feature name) and colnames(which indicates
                       sample ID).
            ctrlInput: str, integrated file for control samples. Same format with caseInput.
            model: Classifier model from sklearn with fit method. None means sklearn.svm.SVC with default parameters.
            LOOCV: Using leave-one-out-cross-validation or not, if True, no model will be saved.
            test_size: if LOOCV is False, please specify train test split proportion.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            caseupstream: upstream output results, used for pipeline.
            ctrlupstream: upstream output results, used for pipeline.
            stepNum: int, step number for folder name.
        """
        if (stepNum is None) and (caseupstream is not None) and (ctrlupstream is None):
            super(classifier, self).__init__(stepNum, caseupstream)
        elif (stepNum is None) and (caseupstream is None) and (ctrlupstream is not None):
            super(classifier, self).__init__(stepNum, ctrlupstream)
        elif (stepNum is None) and (caseupstream is not None) and (ctrlupstream is not None):
            if caseupstream.getStepID() >= ctrlupstream.getStepID():
                super(classifier, self).__init__(stepNum, caseupstream)
            else:
                super(classifier, self).__init__(stepNum, ctrlupstream)
        else:
            super(classifier, self).__init__(stepNum)

        # set caseInput and ctrlInput
        if ((caseupstream is None) and (ctrlupstream is None)) or (caseupstream is True) or (ctrlupstream is True):
            self.setInput("caseInput", caseInput)
            self.setInput("ctrlInput", ctrlInput)
        else:
            Configure2.configureCheck()
            caseupstream.checkFilePath()
            ctrlupstream.checkFilePath()
            if (caseupstream.__class__.__name__ == "computeDMR") and (ctrlupstream.__class__.__name__ == "computeDMR"):
                self.setInput("caseInput", caseupstream.getOutput("casetxtOutput"))
                self.setInput("ctrlInput", ctrlupstream.getOutput("ctrltxtOutput"))
            elif (caseupstream.__class__.__name__ == "runDeconCCN") and (ctrlupstream.__class__.__name__ == "runDeconCCN"):
                self.setInput("caseInput", caseupstream.getOutput("txtOutput"))
                self.setInput("ctrlInput", ctrlupstream.getOutput("txtOutput"))
            else:
                raise commonError("Parameter 'caseupstream' and 'ctrlupstream' must from computeDMR or runDeconCCN.")

        self.checkInputFilePath()

        # set default model
        if model is None:
            model = SVC(gamma="auto")

        # set outputdir
        if (caseupstream is None) and (ctrlupstream is None):
            if outputdir is None:
                self.setOutput(
                    "outputdir", os.path.dirname(os.path.abspath(self.getInput("caseInput")[1])),
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

        self.setOutput("modelOutput", os.path.join(self.getOutput("outputdir"), "prediction.model"))

        # all data
        case_data = pd.read_table(self.getInput("caseInput"), header=0, index_col=0)
        ctrl_data = pd.read_table(self.getInput("ctrlInput"), header=0, index_col=0)
        if not all(case_data.index == ctrl_data.index):
            commonError("Rowname(index) is not the same between case and control.")

        x = np.transpose(pd.concat([case_data, ctrl_data], axis=1).values)
        y = np.asarray(
            list(itertools.repeat(0, len(case_data.columns))) + list(itertools.repeat(1, len(ctrl_data.columns)))
        )

        finishFlag = self.stepInit(caseupstream)

        if not finishFlag:
            if LOOCV:
                print("Start leave-one-out-cross-validation......")
                trueLabel, predLabel = get_cross_validation(model, x, y)
            else:
                print("Start build model......")
                X_train, X_test, y_train, y_test = train_test_split(
                    x, y, test_size=test_size, random_state=np.random.randint(1000000)
                )
                model.fit(X_train, y_train)
                trueLabel = y_test
                predLabel = model.predict(X_test)

            acc = accuracy_score(trueLabel, predLabel)
            print("trueLabel:", trueLabel)
            print("predLabel:", predLabel)
            print("Accuracy:", acc)

            joblib.dump(model, self.getOutput("modelOutput"))

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
