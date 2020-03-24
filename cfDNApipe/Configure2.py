# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 16:45:31 2019

@author: zhang
"""

import os
import urllib.request
from multiprocessing import cpu_count
from .cfDNA_utils import commonError, un_gz, cmdCall
import pandas as pd
from .Configure import Configure

__metaclass__ = type


class Configure2:
    __config = {
        "threads": (cpu_count() / 2),
        "genome": None,
        "refdir": None,
        "outdir": None,
        "data": None,
        "type": "paired",
        "case": "case",
        "ctrl": "ctrl",
    }

    def __init__(self,):
        """
        threads: int, how many thread to use, default: (cpu_count() / 2)
        genome: str, which genome you want to use, 'hg19' or 'hg38'
        refdir: reference folder for aligner (bowtie2 or bismark)
        outdir: overall result folder
        tmpdir: intermediate result folder
        finaldir: most commonly used result folder
        repdir: report result folder
        data: data type, 'WGBS' or 'WGS'
        type: data type, 'paired' or 'single'
        """
        raise commonError("Configure2 can not be initialized")

    # get configure names
    @classmethod
    def getConfigs(cls,):
        return cls.__config.keys()

    # get configure through name
    @classmethod
    def getConfig(cls, key):
        return cls.__config[key]

    # set configure through name
    @classmethod
    def setConfig(cls, key, val):
        if key == "threads":
            Configure2.setThreads(val)
        elif key == "genome":
            Configure2.setGenome(val)
        elif key == "outdir":
            Configure2.setOutDir(val)
        elif key == "refdir":
            Configure2.setRefDir(val)
        elif key == "data":
            Configure2.setData(val)
        elif key == "type":
            Configure2.setType(val)
        elif key == "case":
            Configure2.setCase(val)
        elif key == "ctrl":
            Configure2.setCtrl(val)
        else:
            cls.__config[key] = val

    # set thread
    @classmethod
    def setData(cls, val):
        cls.__config["data"] = val

    # get thread
    @classmethod
    def getData(cls):
        return cls.__config["data"]

    # set thread
    @classmethod
    def setType(cls, val):
        cls.__config["type"] = val

    # get thread
    @classmethod
    def getType(cls):
        return cls.__config["type"]

    # set thread
    @classmethod
    def setThreads(cls, val):
        if Configure2.getData() is None:
            raise commonError("Please set data type before using setThreads.")

        cls.__config["threads"] = val

    # get thread
    @classmethod
    def getThreads(cls):
        return cls.__config["threads"]

    # set reference path
    @classmethod
    def setRefDir(cls, folderPath):
        if Configure2.getGenome() is None:
            raise commonError("Please set genome before using setRefDir.")

        Configure2.checkFolderPath(folderPath)
        cls.__config["refdir"] = folderPath

    # get reference path
    @classmethod
    def getRefDir(cls,):
        return cls.__config["refdir"]

    @classmethod
    def setOutDir(cls, folderPath):
        Configure2.checkFolderPath(folderPath)
        cls.__config["outdir"] = folderPath
        cls.__config["tmpdir"] = os.path.join(
            folderPath, "intermediate_result")
        cls.__config["finaldir"] = os.path.join(folderPath, "final_result")
        cls.__config["repdir"] = os.path.join(folderPath, "report_result")
        cls.__config["casedir"] = os.path.join(
            Configure2.getOutDir(), cls.__config["case"]
        )
        cls.__config["ctrldir"] = os.path.join(
            Configure2.getOutDir(), cls.__config["ctrl"]
        )

    # get overall output path
    @classmethod
    def getOutDir(cls,):
        return cls.__config["outdir"]

    # get intermediate result path
    @classmethod
    def getTmpDir(cls,):
        return cls.__config["tmpdir"]

    # get final result path
    @classmethod
    def getFinalDir(cls,):
        return cls.__config["finaldir"]

    # get report result path
    @classmethod
    def getRepDir(cls,):
        return cls.__config["repdir"]

    # set genome falg
    @classmethod
    def setGenome(cls, val):
        if Configure2.getThreads() is None:
            raise commonError("Please set threads before using setGenome.")

        cls.__config["genome"] = val

    # get genome falg
    @classmethod
    def getGenome(cls):
        return cls.__config["genome"]

    # set case
    @classmethod
    def setCase(cls, val):
        cls.__config["case"] = val

    # get case
    @classmethod
    def getCase(cls):
        return cls.__config["case"]

    # set thread
    @classmethod
    def setCtrl(cls, val):
        cls.__config["ctrl"] = val

    # get thread
    @classmethod
    def getCtrl(cls):
        return cls.__config["ctrl"]

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

    # check folder legency, existence and accessibility
    @staticmethod
    def checkFolderPath(folderPath):
        if not os.path.isdir(os.path.abspath(folderPath)):
            raise commonError(folderPath + " is not an folder.")
        if not os.path.exists(folderPath):
            raise commonError(folderPath + " is not exist.")
        if not (os.access(folderPath, os.X_OK) and os.access(folderPath, os.W_OK)):
            raise commonError(folderPath + " is not accessible.")
        return True

    # create intermediate, final and report folder
    @classmethod
    def pipeFolderInit(cls,):
        Configure2.configureCheck()
        if not os.path.exists(cls.__config["tmpdir"]):
            os.mkdir(cls.__config["tmpdir"])
        if not os.path.exists(cls.__config["finaldir"]):
            os.mkdir(cls.__config["finaldir"])
        if not os.path.exists(cls.__config["repdir"]):
            os.mkdir(cls.__config["repdir"])
        if not os.path.exists(cls.__config["casedir"]):
            os.mkdir(cls.__config["casedir"])
        if not os.path.exists(cls.__config["ctrldir"]):
            os.mkdir(cls.__config["ctrldir"])
        Configure2.checkFolderPath(cls.__config["tmpdir"])
        Configure2.checkFolderPath(cls.__config["finaldir"])
        Configure2.checkFolderPath(cls.__config["repdir"])
        Configure2.checkFolderPath(cls.__config["casedir"])
        Configure2.checkFolderPath(cls.__config["ctrldir"])

    # check configure
    @classmethod
    def configureCheck(cls,):
        if Configure2.getType() is None:
            raise commonError("Please set type configure before using.")
        if Configure2.getData() is None:
            raise commonError("Please set data configure before using.")
        if Configure2.getGenome() is None:
            raise commonError("Please set genome configure before using.")
        if Configure2.getRefDir() is None:
            raise commonError("Please set reference configure before using.")
        if Configure2.getConfig("tmpdir") is None:
            raise commonError("Please set Output configure before using.")
        if Configure2.getConfig("finaldir") is None:
            raise commonError("Please set Output configure before using.")
        if Configure2.getConfig("repdir") is None:
            raise commonError("Please set Output configure before using.")
        if Configure2.getConfig("case") is None:
            raise commonError("Please set case configure before using.")
        if Configure2.getConfig("casedir") is None:
            raise commonError("Please set case configure before using.")
        if Configure2.getConfig("ctrl") is None:
            raise commonError("Please set ctrl configure before using.")
        if Configure2.getConfig("ctrldir") is None:
            raise commonError("Please set ctrl configure before using.")

    # check configure
    @classmethod
    def refCheck(cls, build=False):

        if Configure2.getData() == "WGBS":
            Configure2.bismkrefcheck(build)
            print("Background reference check finished!")
        elif Configure2.getData() == "WGS":
            Configure2.bt2refcheck(build)
            print("Background reference check finished!")
        else:
            print("No reference is specified.")

    # ref check
    @classmethod
    def bismkrefcheck(cls, build):
        # check other reference
        Configure2.genomeRefCheck(build=build)
        Configure2.mappabilityRefCheck(build=build)
        Configure2.cgiRefCheck(build=build)
        Configure2.cytoBandRefCheck(build=build)
        Configure2.ocfRefCheck(build=build)
        # check Bismark reference
        CTfiles = [
            os.path.join(Configure2.getRefDir(),
                         "Bisulfite_Genome/CT_conversion/" + x)
            for x in [
                "BS_CT.1.bt2",
                "BS_CT.2.bt2",
                "BS_CT.3.bt2",
                "BS_CT.4.bt2",
                "BS_CT.rev.1.bt2",
                "BS_CT.rev.2.bt2",
                "genome_mfa.CT_conversion.fa",
            ]
        ]
        BAfiles = [
            os.path.join(Configure2.getRefDir(),
                         "Bisulfite_Genome/GA_conversion/" + x)
            for x in [
                "BS_GA.1.bt2",
                "BS_GA.2.bt2",
                "BS_GA.3.bt2",
                "BS_GA.4.bt2",
                "BS_GA.rev.1.bt2",
                "BS_GA.rev.2.bt2",
                "genome_mfa.GA_conversion.fa",
            ]
        ]
        bismkRef = CTfiles + BAfiles
        if not all(map(os.path.exists, bismkRef)):
            print("Bismark index file do not exist or missing some files!")
            if build:
                cmdline = "bismark_genome_preparation " + Configure2.getRefDir()
                print("Start building bismark reference......")
                print("Now, running " + cmdline)
                cmdCall(cmdline)
                print("Finished!")

    # ref check
    @classmethod
    def bt2refcheck(cls, build):
        # check other reference
        Configure2.genomeRefCheck(build=build)
        Configure2.mappabilityRefCheck(build=build)
        Configure2.cytoBandRefCheck(build=build)
        Configure2.ocfRefCheck(build=build)
        # bowtie2 ref check
        extension = [".1.bt2", ".2.bt2", ".3.bt2",
                     ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
        bt2Ref = [
            os.path.join(Configure2.getRefDir(), Configure2.getGenome() + x)
            for x in extension
        ]
        if not all(map(os.path.exists, bt2Ref)):
            print("Bowtie2 index file do not exist or missing some files!")
            if build:
                cmdline = (
                    "bowtie2-build -f --threads "
                    + str(Configure2.getThreads())
                    + " "
                    + Configure2.getConfig("genome.seq")
                    + " "
                    + os.path.join(Configure2.getRefDir(),
                                   Configure2.getGenome())
                )
                print("Start building Bowtie2 reference......")
                print("Now, running " + cmdline)
                cmdCall(cmdline)
                print("Finished!")

    # check genome reference
    @classmethod
    def genomeRefCheck(cls, build):
        Configure2.setConfig(
            "genome.seq",
            os.path.join(Configure2.getRefDir(),
                         Configure2.getGenome() + ".fa"),
        )
        if not os.path.exists(Configure2.getConfig("genome.seq")):
            print(
                "Reference file "
                + Configure2.getConfig("genome.seq")
                + " do not exist!"
            )
            if build:
                url = (
                    "https://hgdownload.soe.ucsc.edu/goldenPath/"
                    + Configure2.getGenome()
                    + "/bigZips/"
                    + Configure2.getGenome()
                    + ".fa.gz"
                )
                print("Download from URL:" + url + "......")
                urllib.request.urlretrieve(
                    url,
                    os.path.join(
                        Configure2.getRefDir(), Configure2.getGenome() + ".fa.gz"
                    ),
                )
                print("Uncompressing......")
                un_gz(
                    os.path.join(
                        Configure2.getRefDir(), Configure2.getGenome() + ".fa.gz"
                    )
                )
                print("Finished!")

    # check mappability reference
    @classmethod
    def mappabilityRefCheck(cls, build):
        if Configure2.getGenome() == "hg19":
            mapBWfile = os.path.join(
                Configure2.getRefDir(), "wgEncodeCrgMapabilityAlign50mer.bigWig",
            )
            Configure2.setConfig("mappability", mapBWfile)
            if not os.path.exists(mapBWfile):
                print("Reference file " + mapBWfile + " do not exist!")
                if build:
                    url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig"
                    print("Download from URL:" + url + "......")
                    urllib.request.urlretrieve(url, mapBWfile)
        else:
            pass

    # check CpG island reference
    @classmethod
    def cgiRefCheck(cls, build):
        Configure2.setConfig(
            "CGisland",
            os.path.join(
                Configure2.getRefDir(),
                Configure2.getGenome() + "-" + "cpgIslandExt.bed",
            ),
        )
        if not os.path.exists(Configure2.getConfig("CGisland")):
            print(
                "Reference file " +
                Configure2.getConfig("CGisland") + " do not exist!"
            )
            if build:
                url = (
                    "http://hgdownload.soe.ucsc.edu/goldenPath/"
                    + Configure2.getGenome()
                    + "/database/cpgIslandExt.txt.gz"
                )
                print("Download from URL:" + url + "......")
                urllib.request.urlretrieve(
                    url, os.path.join(Configure2.getRefDir(),
                                      "cpgIslandExt.txt.gz")
                )
                print("Uncompressing......")
                un_gz(os.path.join(Configure2.getRefDir(), "cpgIslandExt.txt.gz"))
                regions = pd.read_csv(
                    os.path.join(Configure2.getRefDir(), "cpgIslandExt.txt"),
                    sep="\t",
                    header=None,
                )
                output_regions = regions.iloc[:, [1, 2, 3]]
                output_regions.to_csv(
                    Configure2.getConfig("CGisland"),
                    sep="\t",
                    header=False,
                    index=False,
                )
                print("Finished!")

    # check cytoBand reference
    @classmethod
    def cytoBandRefCheck(cls, build):
        Configure2.setConfig(
            "cytoBand", os.path.join(Configure2.getRefDir(), "cytoBand.txt"),
        )
        if not os.path.exists(Configure2.getConfig("cytoBand")):
            print(
                "Reference file " +
                Configure2.getConfig("cytoBand") + " do not exist!"
            )
            if build:
                url = (
                    "http://hgdownload.cse.ucsc.edu/goldenpath/"
                    + Configure2.getGenome()
                    + "/database/cytoBand.txt.gz"
                )
                print("Download from URL:" + url + "......")
                urllib.request.urlretrieve(
                    url, os.path.join(Configure2.getRefDir(),
                                      "cytoBand.txt.gz"),
                )
                print("Uncompressing......")
                un_gz(os.path.join(Configure2.getRefDir(), "cytoBand.txt.gz"))
                print("Finished!")

    # check OCF reference
    @classmethod
    def ocfRefCheck(cls, build):
        Configure2.setConfig(
            "ocfRef",
            os.path.join(Configure2.getRefDir(),
                         Configure2.getGenome() + ".OCF.bed"),
        )
        if not os.path.exists(Configure2.getConfig("ocfRef")):
            print("Reference file " +
                  Configure2.getConfig("ocfRef") + " do not exist!")
            if build:
                url = (
                    "https://honchkrow.github.io/cfDNApipe/OCF/"
                    + Configure2.getGenome()
                    + ".OCF.bed"
                )
                print("Download from URL:" + url + "......")
                urllib.request.urlretrieve(
                    url, Configure2.getConfig("ocfRef"),
                )
                print("Finished!")


def switchConfigure(confName=None):
    """
    Switch Configure for case and control, these two situation have different output directory
    parameter confName: one of the Configure name defined in Configure2 (case and ctrl)
    """
    Configure.setData(Configure2.getData())
    Configure.setThreads(Configure2.getThreads())
    Configure.setGenome(Configure2.getGenome())
    Configure.setRefDir(Configure2.getRefDir())
    if confName == Configure2.getCase():
        Configure.setOutDir(Configure2.getConfig("casedir"))
    elif confName == Configure2.getCtrl():
        Configure.setOutDir(Configure2.getConfig("ctrldir"))
    else:
        commonError("There is no Configure environment named " + confName + "!")

    Configure.pipeFolderInit()
    Configure.refCheck()


def pipeConfigure(
    threads=(cpu_count() / 2),
    genome=None,
    refdir=None,
    outdir=None,
    data=None,
    type=None,
    case=None,
    ctrl=None,
    build=False,
):
    Configure2.setData(data)
    Configure2.setType(type)
    Configure2.setThreads(threads)
    Configure2.setGenome(genome)
    Configure2.setRefDir(refdir)
    Configure2.setCase(case)
    Configure2.setCtrl(ctrl)
    Configure2.setOutDir(outdir)
    Configure2.pipeFolderInit()
    Configure2.refCheck(build=build)
