# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 09:55:10 2019

@author: zhang-wei
"""


import os
from multiprocessing import cpu_count
from cfDNA_utils import commonError


__metaclass__ = type


class Configure:
    __config = {
        'threads': cpu_count(),
        'genome': None,
        'refdir': None,
        'outdir': None,
        'tmpdir': None,
        'finaldir': None,
        'repdir': None
    }

    def __init__(self,):
        raise commonError('Configure can not be initialized')

    # get configure through name
    @classmethod
    def getConfig(cls, key):
        return cls.__config[key]

    # set configure through name
    @classmethod
    def setConfig(cls, key, val):
        if key == 'threads':
            Configure.setThreads(val)
        elif key == 'genome':
            Configure.setGenome(val)
        elif key == 'outdir':
            Configure.setOutDir(val)
        elif key == 'refdir':
            Configure.setRefDir(val)
        else:
            cls.__config[key] = val

    # set thread
    @classmethod
    def setThreads(cls, val):
        cls.__config['threads'] = val

    # get thread
    @classmethod
    def getThreads(cls):
        return cls.__config['threads']

    # set reference path
    @classmethod
    def setRefDir(cls, folderPath):
        Configure.checkFolderPath(folderPath)
        cls.__config['refdir'] = folderPath

    # get reference path
    @classmethod
    def getRefDir(cls,):
        return cls.__config['refdir']

    # set overall output directory and sub dir
    @classmethod
    def setOutDir(cls, folderPath):
        Configure.checkFolderPath(folderPath)
        cls.__config['outdir'] = folderPath
        cls.__config['tmpdir'] = folderPath + '/intermediate_result'
        cls.__config['finaldir'] = folderPath + '/final_result'
        cls.__config['repdir'] = folderPath + '/report_result'

    # get overall output path
    @classmethod
    def getOutDir(cls, ):
        return cls.__config['outdir']

    # get intermediate result path
    @classmethod
    def getTmpDir(cls, ):
        return cls.__config['tmpdir']

    # get final result path
    @classmethod
    def getFinalDir(cls, ):
        return cls.__config['finaldir']
    
    # get report result path
    @classmethod
    def getRepDir(cls, ):
        return cls.__config['repdir']

    # create intermediate, final and report folder
    @classmethod
    def pipeFolderInit(cls, ):
        Configure.configureCheck()
        if not os.path.exists(cls.__config['tmpdir']):
            os.mkdir(cls.__config['tmpdir'])
        if not os.path.exists(cls.__config['finaldir']):
            os.mkdir(cls.__config['finaldir'])
        if not os.path.exists(cls.__config['repdir']):
            os.mkdir(cls.__config['repdir'])
        Configure.checkFolderPath(cls.__config['tmpdir'])
        Configure.checkFolderPath(cls.__config['finaldir'])
        Configure.checkFolderPath(cls.__config['repdir'])

    # set genome falg
    @classmethod
    def setGenome(cls, val):
        cls.__config['genome'] = val

    # get genome falg
    @classmethod
    def getGenome(cls):
        return cls.__config['genome']

    # check folder legency, existence and sccessibility
    @staticmethod
    def checkFolderPath(folderPath):
        if not os.path.isdir(os.path.abspath(folderPath)):
            raise commonError(folderPath + " is not an folder.")
        if not os.path.exists(folderPath):
            raise commonError(folderPath + " is not exist.")
        if not (os.access(folderPath, os.X_OK) and os.access(folderPath, os.W_OK)):
            raise commonError(folderPath + " is not accessible.")
        return True

    # get intermediate result path
    @classmethod
    def getTmpPath(cls, foldOrFileName):
        if isinstance(foldOrFileName, list):
            result = []
            for name in foldOrFileName:
                result.append(os.path.join(cls.getTmpDir(), name))
            return result
        else:
            return os.path.join(cls.getTmpDir(), foldOrFileName)
        
    # check configure 
    @classmethod
    def configureCheck(cls,):
        if Configure.getGenome() is None:
            raise commonError("Please set genome configure before using.")
        if Configure.getRefDir() is None:
            raise commonError("Please set reference configure before using.")
        if Configure.getConfig('tmpdir') is None:
            raise commonError("Please set Output configure before using.")
        if Configure.getConfig('finaldir') is None:
            raise commonError("Please set Output configure before using.")
        if Configure.getConfig('repdir') is None:
            raise commonError("Please set Output configure before using.")
        if Configure.getConfig('outdir') is None:
            raise commonError("Please set Output configure before using.")





















