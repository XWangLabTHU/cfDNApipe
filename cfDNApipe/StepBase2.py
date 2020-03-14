# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 09:55:10 2019

@author: zhang-wei
"""


from .StepBase import StepBase
from .Configure2 import Configure2


__metaclass__ = type


class StepBase2(StepBase):
    def __init__(self, stepNum=None, prevStep=None):
        if stepNum is not None:
            self.__stepID = stepNum
        elif (stepNum is None) and (prevStep is not None):
            self.__stepID = prevStep.getStepID() + 1
        else:
            self.__stepID = 0
        self.inputs = {}
        self.outputs = {}
        self.params = {}
        self.logpath = {}
        self.__isFinished = False

    # get this step folder path
    def getStepFolderPath(self,):
        return Configure2.getTmpPath(self.getStepFolderName())

    # get configure value
    def getConfigVal(key):
        return Configure2.getConfig(key)
