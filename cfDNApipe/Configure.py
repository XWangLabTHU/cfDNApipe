# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 09:55:10 2019

@author: zhang-wei
"""


import os, urllib.request
from multiprocessing import cpu_count
from .cfDNA_utils import commonError, un_gz, cmdCall
import pandas as pd


__metaclass__ = type


class Configure:
    __config = {
        'threads': (cpu_count() / 2),
        'genome': None,
        'refdir': None,
        'outdir': None,
        'tmpdir': None,
        'finaldir': None,
        'repdir': None,
        'data': None,
        'type': 'paired'
    }

    def __init__(self,):
        '''
        threads: int, how many thread to use, default: (cpu_count() / 2)
        genome: str, which genome you want to use, 'hg19' or 'hg38'
        refdir: reference folder for aligner (bowtie2 or bismark)
        outdir: overall result folder
        tmpdir: intermediate result folder
        finaldir: most commonly used result folder
        repdir: report result folder
        data: data type, 'WGBS' or 'WGS'
        type: data type, 'paired' or 'single'
        '''
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
        elif key == 'data':
            Configure.setData(val)
        elif key == 'type':
            Configure.setType(val)
        else:
            cls.__config[key] = val

    # set thread
    @classmethod
    def setData(cls, val):
        cls.__config['data'] = val

    # get thread
    @classmethod
    def getData(cls):
        return cls.__config['data']
    
    # set thread
    @classmethod
    def setType(cls, val):
        cls.__config['type'] = val

    # get thread
    @classmethod
    def getType(cls):
        return cls.__config['type']

    # set thread
    @classmethod
    def setThreads(cls, val):
        if Configure.getData() is None:
            raise commonError("Please set data type before using setThreads.")
            
        cls.__config['threads'] = val

    # get thread
    @classmethod
    def getThreads(cls):
        return cls.__config['threads']

    # set reference path
    @classmethod
    def setRefDir(cls, folderPath):
        if Configure.getGenome() is None:
            raise commonError("Please set genome before using setRefDir.")
            
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
        if Configure.getThreads() is None:
            raise commonError("Please set threads before using setGenome.")
        
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
        if Configure.getType() is None:
            raise commonError("Please set type configure before using.")
        if Configure.getData() is None:
            raise commonError("Please set data configure before using.")
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


    # check configure 
    @classmethod
    def refCheck(cls, build = False):
        
        if Configure.getData() == 'WGBS':
            Configure.bismkrefcheck(build)
            print("Bismark reference check finished!")
        elif Configure.getData() == 'WGS':
            Configure.bt2refcheck(build)
            print("Bowtie2 reference check finished!")
        else:
            print("No reference is specified.")


    # ref check
    @classmethod
    def bismkrefcheck(cls, build):
        # check genome reference
        fafile = os.path.join(Configure.getRefDir(), Configure.getGenome() + '.fa')
        if not os.path.exists(fafile):
            print("Reference file " + fafile + " do not exist!")
            if build == True:
                url = "https://hgdownload.soe.ucsc.edu/goldenPath/" + Configure.getGenome() + "/bigZips/" + Configure.getGenome() + ".fa.gz"
                print("Download from URL:" + url + "......")
                urllib.request.urlretrieve(url, os.path.join(Configure.getRefDir(), Configure.getGenome() + '.fa.gz'))
                print("Uncompressing......")
                un_gz(os.path.join(Configure.getRefDir(), Configure.getGenome() + '.fa.gz'))
                print("Finished!")
            
        # check CpG island reference
        CpGislandfile = os.path.join(Configure.getRefDir(), Configure.getGenome() + 'cpgIslandExt.bed')
        if not os.path.exists(CpGislandfile):
            print("Reference file " + CpGislandfile + " do not exist!")
            if build == True:
                url = "http://hgdownload.soe.ucsc.edu/goldenPath/" + Configure.getGenome() + "/database/cpgIslandExt.txt.gz"
                print("Download from URL:" + url + "......")
                urllib.request.urlretrieve(url, os.path.join(Configure.getRefDir(), 'cpgIslandExt.txt.gz'))
                print("Uncompressing......")
                un_gz(os.path.join(Configure.getRefDir(), 'cpgIslandExt.txt.gz'))
                regions = pd.read_csv(os.path.join(Configure.getRefDir(), 'cpgIslandExt.txt'), sep = "\t", header = None)
                output_regions = regions.iloc[:, [1, 2, 3]]
                output_regions.to_csv(CpGislandfile, sep = "\t", header = False, index = False)
                print("Finished!")
        
        # check Bismark reference
        CTfiles = [os.path.join(Configure.getRefDir(), 'Bisulfite_Genome/CT_conversion/' + x) for x in ['BS_CT.1.bt2', 'BS_CT.2.bt2', 'BS_CT.3.bt2', 'BS_CT.4.bt2', 'BS_CT.rev.1.bt2', 'BS_CT.rev.2.bt2', 'genome_mfa.CT_conversion.fa']]
        BAfiles = [os.path.join(Configure.getRefDir(), 'Bisulfite_Genome/GA_conversion/' + x) for x in ['BS_GA.1.bt2', 'BS_GA.2.bt2', 'BS_GA.3.bt2', 'BS_GA.4.bt2', 'BS_GA.rev.1.bt2', 'BS_GA.rev.2.bt2', 'genome_mfa.GA_conversion.fa']]
        bismkRef = CTfiles + BAfiles
        if not all(map(os.path.exists, bismkRef)):
            print('Bismark index file do not exist or missing some files!')
            if build == True:
                cmdline = "bismark_genome_preparation " + Configure.getRefDir()
                print("Start building bismark reference......")
                print("Now, running " + cmdline)
                cmdCall(cmdline)
                print("Finished!")


    # ref check
    @classmethod
    def bt2refcheck(cls, build):
        # check genome reference
        fafile = os.path.join(Configure.getRefDir(), Configure.getGenome() + '.fa')
        if not os.path.exists(fafile):
            print("Reference file " + fafile + " do not exist!")
            if build == True:
                url = "https://hgdownload.soe.ucsc.edu/goldenPath/" + Configure.getGenome() + "/bigZips/" + Configure.getGenome() + ".fa.gz"
                print("Download from URL:" + url + "......")
                urllib.request.urlretrieve(url, os.path.join(Configure.getRefDir(), Configure.getGenome() + '.fa.gz'))
                print("Uncompressing......")
                un_gz(os.path.join(Configure.getRefDir(), Configure.getGenome() + '.fa.gz'))
                print("Finished!")
                
        
        extension = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
        bt2Ref = [Configure.getRefDir() + x for x in extension]
        if not all(map(os.path.exists, bt2Ref)):
            print('Bismark index file do not exist or missing some files!')
            if build == True:
                cmdline = "bowtie2-build -f --threads " + str(Configure.getThreads()) + " " + fafile + " " + os.path.join(Configure.getRefDir(), Configure.getGenome())
                print("Start building Bowtie2 reference......")
                print("Now, running " + cmdline)
                cmdCall(cmdline)
                print("Finished!")




















































