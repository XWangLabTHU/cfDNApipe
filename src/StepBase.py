# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 09:55:10 2019

@author: zhang-wei
"""


from Configure import Configure
import os
from hashlib import md5
from cfDNA_utils import commonError, flatten, isAlphaOrDigit, rmEndString
import time
import ast
import subprocess, sys


__metaclass__ = type



class StepBase:
    def __init__(self, stepCounter = 0):
        self.inputs = {}
        self.outputs = {}
        self.params = {}
        self.stepCounter = stepCounter
        self.logpath = {}
        self.__isFinished = False
        self.__stepID = self.regStepID()

    # set input parameters, all the input will be absolute Path
    def setInput(self, inputName, inputValue):
        if isinstance(inputName, list):
            if len(inputName) != len(inputValue):
                raise commonError('Number of input name and value not equal.')
            values = self.absolutePath(inputValue)
            for name, value in zip(inputName, values):
                self.inputs[name] = value
        else:
            self.inputs[inputName] = self.absolutePath(inputValue)

    # get input value by dict name (key)
    def getInput(self, inputName):
        return self.inputs[inputName]

    # get all input keys
    def getInputs(self,):
        return list(self.inputs.keys())

    # set output parameters
    def setOutput(self, outputName, outputValue):
        if isinstance(outputName, list):
            if len(outputName) != len(outputValue):
                raise commonError('Number of output key name and value not equal.')
            values = self.absolutePath(outputValue)
            for name, value in zip(outputName, values):
                self.outputs[name] = value
        else:
            self.outputs[outputName] = self.absolutePath(outputValue)

    # get output value by dict name (key)
    def getOutput(self, outputName):
        return self.outputs[outputName]

    # get all output keys
    def getOutputs(self,):
        return list(self.outputs.keys())

    # check input and output file path
    def checkFilePath(self, checkExist = True):
        self.checkOutputFilePath(checkExist)
        self.checkInputFilePath(checkExist)   

    # check output file path
    def checkOutputFilePath(self, checkExist = True):
        for key in self.outputs.keys():
            outPaths = self.convertToList(self.outputs[key])
            self.checkFilePathList(outPaths, 'output', key, checkExist)

    # check input file path
    def checkInputFilePath(self, checkExist = True):
        for key in self.inputs.keys():
            inPaths = self.convertToList(self.inputs[key])
            self.checkFilePathList(inPaths, 'input', key, checkExist)

    # check path in a list
    def checkFilePathList(self, filePathList, iodict, key = None, checkExist = True):
        for filePath in filePathList:
            if filePath is None:
                raise commonError('File path of ' + iodict + ' ' + key + ' can not be None.')
            if checkExist:
                if not os.path.exists(filePath):
                    raise commonError('File path of ' + iodict + ' ' + key + ' not found: ' + filePath + '.')

    # add 1 to StepID
    def regStepID(self,):
        self.stepCounter = self.stepCounter + 1
        return self.stepCounter

    # get stepID
    def getStepID(self,):
        return self.__stepID

    # get this Step folder name
    def getStepFolderName(self,):
        return 'step_' + str(self.getStepID()).zfill(2) + '_' + self.__class__.__name__

    # get this step folder path
    def getStepFolderPath(self,):
        return Configure.getTmpPath(self.getStepFolderName())

    # create this step folder
    def stepFolderInit(self,):
        if not os.path.exists(self.getStepFolderPath()):
            os.mkdir(self.getStepFolderPath())
    
    # step init, include create step folder, log file and record file
    def stepInit(self, upstream = None):
        if upstream is not None:
            self.stepFolderInit()
            self.setPipeLogPath()
            self.setPipeRecPath()
        else:
            self.setLogPath(os.path.join(self.getOutput('outputdir'), self.getLogName()))
            self.setRecPath(os.path.join(self.getOutput('outputdir'), self.getRecName()))
        
        if os.path.exists(self.getLogPath()):
            finishFlag = self.checkFinish()
        else:
            finishFlag = False
            
        if not finishFlag:
            self.createLog(overwrite = True)
            self.createRec(overwrite = True)
        
        return finishFlag

    # get log file name
    def getLogName(self,):
        return self.__class__.__name__ + '.' + self.getParaMD5code() + '.log'

    # set log path
    def setLogPath(self, path):
        self.logpath = path

    # set pipeline log path
    def setPipeLogPath(self, ):
        self.logpath = os.path.join(self.getStepFolderPath(), self.getLogName())

    # get log file path
    def getLogPath(self,):
        return self.logpath
    
    # create log file
    def createLog(self, overwrite = False):
        if not os.path.exists(self.logpath):
            open(self.logpath, 'a').close()
        else:
            if overwrite:
                open(self.logpath, 'w').close()
            else:
                pass

    # write log
    def writeLogLines(self, strlines):
        if not os.path.exists(self.logpath):
            raise commonError("can not write log when log file is not created!")
        if not isinstance(strlines, list):
            strlines = [strlines]
            
        strlines = list(flatten(strlines))
        logfile = open(self.logpath, 'a')
        mess = '\t'.join([str(x) for x in strlines]) + '\n'
        logfile.writelines(mess)
        logfile.close()
        
    # get log file value
    def getLogValue(self,):
        logdict = {}
        with open(self.logpath) as f:
            for line in f:
                tok = line.split()
                logdict[tok[0]] = tok[1:]
                
        return logdict
    
    # get record file name
    def getRecName(self,):
        return self.__class__.__name__ + '.' + self.getParaMD5code() + '-record.txt'
    
    # set record path
    def setRecPath(self, path):
        self.recpath = path

    # set pipeline record path
    def setPipeRecPath(self, ):
        self.recpath = os.path.join(self.getStepFolderPath(), self.getRecName())
        
    # get log file path
    def getRecPath(self,):
        return self.recpath

    # create recort file
    def createRec(self, overwrite = False):
        if not os.path.exists(self.recpath):
            open(self.recpath, 'a').close()
        else:
            if overwrite:
                open(self.recpath, 'w').close()
            else:
                pass

    # write log
    def writeRec(self, mess):
        recfile = open(self.recpath, 'a')
        recfile.writelines(mess + '\n')
        recfile.close()

    # check whether finished
    def checkFinish(self,):  
         logdict = self.getLogValue()
         try:
             flag = ast.literal_eval(logdict['FinishedOrNot'][0])
         except KeyError:
             return False
         else:
             return flag

    # set input and output parameters
    def setParam(self, paramName, paramValue):
        self.params[paramName] = paramValue

    # get input and output parameters
    def getParam(self, paramName):
        return self.params[paramName]

    # get all parameter keys
    def getParams(self,):
        return list(self.params.keys())

    # get time
    def getCurTime(self,):
        return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))

    # what a value is?
    def getIOtype(self, value):
        if value is None:
            return None
        value = os.path.abspath(value)
        for key in self.inputs.keys():
            for path in self.convertToList(self.inputs[key]):
                apath = os.path.abspath(path)
                if apath.startswith(value):
                    if os.path.isdir(value): 
                        return 'inputDir'
                    elif os.path.isfile(value):
                        return 'inputFile'
                    elif os.path.isdir(os.path.dirname(value)):
                        return 'inputPrefix'
        for key in self.outputs.keys():
            for path in self.convertToList(self.outputs[key]):
                apath = os.path.abspath(path)
                if apath.startswith(value):                    
                    if value == apath:
                        return 'outputFile'
        for key in self.outputs.keys():
            for path in self.convertToList(self.outputs[key]):
                apath = os.path.abspath(path)
                if apath.startswith(value):                    
                    if apath[0:(len(value)+1)] == os.path.sep:
                        return 'outputDir'
                    else:                        
                        return 'outputPrefix'
        return None

    # get file absolute path
    def absolutePath(self, pathOrPathList):
        if pathOrPathList is None:
            return None
        elif isinstance(pathOrPathList, list):
            return [os.path.abspath(s) for s in pathOrPathList]
        else:
            return os.path.abspath(pathOrPathList)

    # convert a odj to a list
    def convertToList(self, obj):
        if not isinstance(obj, list):
            return [obj]
        else:
            return obj

    # get single input file name
    def getFileNamePrefix(self, fileName):
        return os.path.split(fileName)[-1]

    # get paired input file name prefix
    def getMaxFileNamePrefix(self, file1, file2):
        file1 = self.getFileNamePrefix(file1)
        file2 = self.getFileNamePrefix(file2)
        len1 = len(file1)
        len2 = len(file2)
        for i in range(min(len1, len2)):
            if file1[i] != file2[i]:
                break
        if i == 0:
            return ''
        elif i == (min(len1, len2) - 1):
            tmp_str = file1[: i + 1]
        else:
            tmp_str = file1[: i]
        
        for k in reversed(range(len(tmp_str))):
            if isAlphaOrDigit(tmp_str[k]):
                final_name = tmp_str[: k + 1]
                final_name = rmEndString(final_name, ['.pair'])
                return final_name
            else:
                k = k - 1
        
        raise commonError('File names must contain at least one alphabet or number.')

    # single prefix
    def getMaxFileNamePrefixV2(self, file):
        final_name = os.path.splitext(os.path.basename(file))[0]
        final_name = rmEndString(final_name, ['_pe', '-sorted', '-rmdup'])
        return final_name

    # get file name and size
    def getFileNameAndSize(self, filePath, fileSize = True):
        namesizelist = []
        if not isinstance(filePath, list):
            filePath = [filePath]
        for filepath in filePath:
            filename = os.path.split(filepath)[-1]
            if fileSize: 
                filesize = os.path.getsize(filepath)
            else:
                filesize = ''
            namesizelist.append(filename + str(filesize))
        return namesizelist

    # get md5 code
    def getParaMD5code(self,):
        checklist = ['']
        checklist1 = []
        checklist2 = []
        checklist3 = []
        for key in self.inputs.keys():
            checklist1.extend(self.getFileNameAndSize(self.inputs[key]))        
        checklist1.sort()
        checklist.extend(checklist1) 
        for key in self.outputs.keys():
            checklist2.extend(self.getFileNameAndSize(self.outputs[key], fileSize = False))
        checklist2.sort()
        checklist.extend(checklist2)
        keys = list(self.params.keys())
        keys.sort()
        for key in keys:
            checklist3.append(self.params[key])
        checklist.extend(str(checklist3))
        checklist = ''.join(checklist)
        return md5(checklist.encode('utf8')).hexdigest()[0:8]

    # get configure value
    def getConfigVal(key):
        return Configure.getConfig(key)

    # command create
    def cmdCreate(self, cmdlist):
        if not isinstance(cmdlist, list):
            raise commonError("Parameter 'cmdlist' must be a list!")
        else:
            tmp_cmd = []
            for tmp_value in cmdlist:
                if isinstance(tmp_value, dict):
                    for k, v in tmp_value.items():
                        if isinstance(v, bool) and v:
                            tmp_cmd.append(k)
                        elif isinstance(v, bool) and (not v):
                            pass
                        else:
                            tmp_cmd.append(k)
                            tmp_cmd.append(v)
                elif isinstance(tmp_value, list):
                    for l in tmp_value:
                        tmp_cmd.append(l)
                else:
                    tmp_cmd.append(tmp_value)
        
        cmd = ' '.join([str(x) for x in tmp_cmd])
        return cmd

    # run the command line
    def run(self,):
        self.writeRec('#############################################################################################')
        if isinstance(self.getParam('cmd'), list):  # cmd is a list
            for idx, cmd in enumerate(self.getParam('cmd')):
                self.writeLogLines(['Cmd_No.' + str(idx), cmd])
                self.writeRec('Cmd: {}'.format(cmd))
                print('Now, running command: {}'.format(cmd))
                proc = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, universal_newlines = True)
                while True:
                    nextline = proc.stdout.readline()
                    if (nextline == '') and (proc.poll() is not None):
                        break
                    sys.stdout.write(nextline)
                    self.writeRec(nextline.strip())
                    sys.stdout.flush()
                    
                output, error = proc.communicate()
                exitCode = proc.returncode
                # catch error
                if exitCode != 0:
                    self.writeRec('Exit code: {}'.format(exitCode))
                    self.writeRec('Exit cmd: {}'.format(cmd))
                    self.writeLogLines(['FinishedOrNot', 'False'])
                    raise commonError('**********CMD running error**********')
                
                self.writeRec('#############################################################################################')

        else:                                       # cmd is a string
            self.writeLogLines(['Cmd', self.getParam('cmd')])
            self.writeRec('Cmd: {}'.format(self.getParam('cmd')))
            print('Now, running command: {}'.format(self.getParam('cmd')))
            proc = subprocess.Popen(self.getParam('cmd'), shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, universal_newlines = True)
            while True:
                nextline = proc.stdout.readline()
                if (nextline == '') and (proc.poll() is not None):
                    break
                sys.stdout.write(nextline)
                self.writeRec(nextline.strip("\n"))
                sys.stdout.flush()

            output, error = proc.communicate()
            exitCode = proc.returncode
            # catch error
            if exitCode != 0:
                self.writeRec('Exit code: {}'.format(exitCode))
                self.writeRec('Exit cmd: {}'.format(self.getParam('cmd')))
                self.writeLogLines(['FinishedOrNot', 'False'])
                raise commonError('**********CMD running error**********')
            
            self.writeRec('#############################################################################################')


    # excute program
    def excute(self, finishFlag, runFlag = True):
        if finishFlag:
            print("***************************************************************************************")
            print("***************************Program finished before, skip*******************************")
            print("********If you want run it again, please change parameters or delete log file**********")
            print("***************************************************************************************")
        else:
            self.writeLogLines(['Classname', self.__class__.__name__])
            self.writeLogLines(['Start_time', self.getCurTime()])
            self.writeLogLines(['Input_files', self.inputs.values()])
            self.writeLogLines(['Outputs', self.outputs.values()])
            self.writeLogLines(['Log_file', self.getLogPath()])
            self.writeLogLines(['Record_file', self.getRecPath()])
            
            if runFlag:
                self.run()
            else:
                self.writeRec('There is nothing to be recorded in this program.')
            
            self.writeLogLines(['End_time', self.getCurTime()])
            self.writeLogLines(['FinishedOrNot', 'True'])





































