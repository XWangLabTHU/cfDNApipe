# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 19:32:30 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""


from cfDNApipe import *

# set global configure
pipeConfigure(
    threads=100,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bismark",
    outdir="/home/zhangwei/pipeline-WGBS-comp",
    data="WGBS",
    type="paired",
    case="case",
    ctrl="ctrl",
    build=True,
)

Configure2.getConfigs()  # get all configures

# global configure
Configure2.getGenome()  # get reference genome
Configure2.getOutDir()  # get global output dir
Configure2.getRefDir()  # get reference dir
Configure2.getRepDir()  # get global report result dir
Configure2.getFinalDir()  # get global final result dir
Configure2.getTmpDir()  # get global intermediate result dir
Configure2.getThreads()  # get threads
Configure2.getData()  # get data type
Configure2.getType()  # get sequence type
Configure2.getCase()  # get case flag
Configure2.getCtrl()  # get ctrl flag

Configure2.getConfig("genome.seq")  # genome sequence fasta file
Configure2.getConfig("mappability")  # genome mappability file, only for hg19
Configure2.getConfig("CGisland")  # CpGisland bed file


# case processing, case configure will not be initialized bedfore running switchConfigure
switchConfigure("case")

Configure.getConfigs()  # get all configures

Configure.getConfig("genome.seq")  # genome sequence fasta file
Configure.getConfig("mappability")  # genome mappability file, only for hg19
Configure.getConfig("CGisland")  # CpGisland bed file

# ctrl processing, ctrl configure will not be initialized bedfore running switchConfigure
switchConfigure("ctrl")









from cfDNApipe import *