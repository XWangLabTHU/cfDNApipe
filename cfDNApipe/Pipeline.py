# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:26:39 2019

@author: zhang
"""

from .cfDNA_utils import *
from .Configure import *
from .Configure2 import *
from .Fun_adapterremoval import *
from .Fun_inputProcess import *
from .Fun_bam2bed import *
from .Fun_bamsort import *
from .Fun_bismark import *
from .Fun_bowtie2 import *
from .Fun_fastqc import *
from .Fun_fragLen import *
from .Fun_identifyAdapter import *
from .Fun_inputProcess import *
from .Fun_OCF import *
from .Fun_rmDuplicate import *
from .Fun_sequenceTrans import *
from .Pipeline import *
from .StepBase import *
from .StepBase2 import *
from .report_generator import *


def cfDNAWGBS(
    inputFolder=None,
    fastq1=None,
    fastq2=None,
    adapter1=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
    adapter2=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"],
    fastqcOP=None,
    idAdapter=True,
    idAdOP=None,
    rmAdapter=True,
    rmAdOP={"--qualitybase": 33, "--gzip": True},
    bismarkOP={"-q": True, "--phred33-quals": True, "--bowtie2": True, "--un": True,},
    dudup=True,
    dudupOP=None,
    extractMethyOP={
        "--no_overlap": True,
        "--report": True,
        "--no_header": True,
        "--gzip": True,
        "--bedGraph": True,
        "--zero_based": True,
    },
    methyRegion=None,
    CNV=False,
    fragProfile=False,
    verbose=False,
):
    print("Now, running cell free DNA WGBS data processing pipeline.")

    # set input
    if inputFolder is not None:
        res_inputprocess = inputprocess(inputFolder=inputFolder)
    else:
        res_inputprocess = inputprocess(fqInput1=fastq1, fqInput2=fastq2, verbose=verbose)

    # fastqc
    res_fastqc = fastqc(upstream=res_inputprocess, other_params=fastqcOP, verbose=verbose)

    # identify and remove adapters,  bismark
    if rmAdapter:
        if idAdapter:
            res_identifyAdapter = identifyAdapter(upstream=res_inputprocess, other_params=idAdOP, verbose=verbose)
            res_adapterremoval = adapterremoval(upstream=res_identifyAdapter, other_params=rmAdOP, verbose=verbose)
        else:
            res_adapterremoval = adapterremoval(
                upstream=res_inputprocess, adapter1=adapter1, adapter2=adapter2, other_params=rmAdOP, verbose=verbose
            )

        res_bismark = bismark(upstream=res_adapterremoval, other_params=bismarkOP, verbose=False)

    else:
        print("Skip removing adapters.")
        res_bismark = bismark(upstream=res_inputprocess, other_params=bismarkOP, verbose=False)

    # redup and extract methy
    if dudup:
        res_deduplicate = bismark_deduplicate(upstream=res_bismark, other_params=dudupOP, verbose=False)
        res_methyextract = bismark_methylation_extractor(
            upstream=res_deduplicate, other_params=extractMethyOP, verbose=False
        )
        res_bamsort = bamsort(upstream=res_deduplicate, verbose=False)
    else:
        res_methyextract = bismark_methylation_extractor(
            upstream=res_bismark, other_params=extractMethyOP, verbose=False
        )
        res_bamsort = bamsort(upstream=res_bismark, verbose=False)

    res_compressMethy = compress_methyl(upstream=res_methyextract, verbose=False)
    res_calMethy = calculate_methyl(upstream=res_compressMethy, bedInput=methyRegion, verbose=False)

    res_bam2bed = bam2bed(upstream=res_bamsort, verbose=False)
    res_fraglenplot = fraglenplot(upstream=res_bam2bed, verbose=False)

    if CNV:
        cnv_readCounter = runCounter(upstream=res_bamsort, filetype=1, verbose=False, stepNum="CNV01")
        cnv_gcCounter = runCounter(filetype=0, upstream=True, verbose=False, stepNum="CNV02")
        cnv_GCCorrect = GCCorrect(
            readupstream=cnv_readCounter, gcupstream=cnv_gcCounter, verbose=False, stepNum="CNV03"
        )
    else:
        print("Skip CNV analysis.")

    if fragProfile:
        fp_fragCounter = fpCounter(upstream=res_bam2bed, verbose=False, stepNum="FP02")
        fp_gcCounter = runCounter(filetype=0, binlen=5000000, upstream=True, verbose=False, stepNum="FP01")
        fp_GCCorrect = GCCorrect(
            readupstream=fp_fragCounter, gcupstream=fp_gcCounter, readtype=2, corrkey="-", verbose=False, stepNum="FP03"
        )
