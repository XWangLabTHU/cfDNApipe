# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 19:32:30 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""

class a:
    stepCounter = 0
    
    def __init__(self, initStep = False):
        if initStep:
            self.initStepId()
        
        self.__stepID = self.regStepID()
    
    def regStepID(self,):
        a.stepCounter = a.stepCounter + 1
        return a.stepCounter
    
    def getStepId(self,):
        return(self.__stepID)
        
    def initStepId(self,):
        a.stepCounter = 0





class b(a):
    def __init__(self, initStep = False):
        super(b, self).__init__(initStep)

class c(a):
    def __init__(self, initStep = False):
        super(c, self).__init__(initStep)


b1 = b()

c1 = c()

b2 = b()

c2 = c(initStep = True)