# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:26:39 2019

@author: zhang
"""

from .cfDNA_utils import *
from .Configure import *
from .Fun_adapterremoval import *
from .Fun_bam2bed import *
from .Fun_bamsort import *
from .Fun_bismark import *
from .Fun_bowtie2 import *
from .Fun_fastqc import *
from .Fun_fragLen import *
from .Fun_identifyAdapter import *
from .Fun_inputProcess import *
from .Fun_methyl import *
from .Fun_OCF import *
from .Fun_rmDuplicate import *
from .Fun_sequenceTrans import *
from .Fun_addRG import *
from .Pipeline import *
from .StepBase import *
from .report_generator import *


def WGBSConfigure(refDir, outDir, threads=1, genome="hg19", init=True, build=True):
    Configure.setData("WGBS")
    Configure.setThreads(threads)
    Configure.setGenome(genome)
    Configure.setRefDir(refDir)
    Configure.setOutDir(outDir)
    if init:
        Configure.pipeFolderInit()

    if build:
        Configure.refCheck(build=True)
    else:
        Configure.refCheck(build=False)


def cfDNAWGBS(
    inputs=None,
    adapter1=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
    adapter2=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"],
    fastqcOP=None,
    idAdOP=None,
    rmAdOP={"--qualitybase": 33, "--gzip": True},
    bismarkOP={
        "-q": True,
        "--phred33-quals": True,
        "-N": 1,
        "-X": 2000,
        "--bowtie2": True,
        "--no_dovetail": True,
    },
    maxLimit=250,
    MLR=None,
):
    # get inputs
    if isinstance(inputs, list):
        res1 = inputprocess(fqInput1=inputs[0], fqInput2=inputs[1])
    elif isinstance(inputs, str):
        res1 = inputprocess(inputFolder=inputs)
    else:
        commonError("Parameter inputs must be a list (length 2) or string!")

    # fastqc
    res2 = fastqc(upstream=res1, other_params=fastqcOP)

    # identify adapters
    res3 = identifyAdapter(upstream=res1, formerrun=res2, other_params=idAdOP)

    # remove adapters
    res4 = adapterremoval(
        upstream=res3, adapter1=adapter1, adapter2=adapter2, other_params=rmAdOP
    )

    # mapping
    res5 = bismark(upstream=res4, other_params=bismarkOP)

    # sort bam
    res6 = bamsort(upstream=res5)

    # remove duplicates
    res7 = rmduplicate(upstream=res6)

    # bam to bed
    res8 = bam2bed(upstream=res7)

    # plot fragment length distribution
    res9 = fraglenplot(upstream=res8, maxLimit=maxLimit)

    # compute methylation level
    res10 = computemethyl(upstream=res7, formerrun=res9, bedInput=MLR)

    # add read group
    res11 = addRG(upstream=res7, formerrun=res10)

    # generate report
    report_generator(
        fastqcRes=res2,
        identifyAdapterRes=res3,
        bismarkRes=res5,
        rmduplicateRes=res7,
        fraglenplotRes=res9,
    )

    return True
