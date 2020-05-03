# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:26:39 2019

@author: zhang Wei
"""

from box import Box
from .Fun_inputProcess import inputprocess
from .Fun_fastqc import fastqc
from .Fun_identifyAdapter import identifyAdapter
from .Fun_adapterremoval import adapterremoval
from .Fun_bismark import bismark
from .Fun_bismarkdedup import bismark_deduplicate
from .Fun_bismarkmethylex import bismark_methylation_extractor
from .Fun_bamsort import bamsort
from .Fun_compressmethyl import compress_methyl
from .Fun_calcmethyl_v2 import calculate_methyl
from .Fun_bam2bed import bam2bed
from .Fun_counter import runCounter
from .Fun_fragLen import fraglenplot
from .Fun_GCcorrect import GCCorrect
from .Fun_fpcount import fpCounter


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
    """
    This function is used for removing duplicates in WGS data.
    Note: this function is calling picard.

    rmduplicate(bamInput=None, outputdir=None, threads=1, stepNum=None, upstream=None, verbose=True)
    {P}arameters:
        bamInput: list, bam file input.
        outputdir: str, output result folder, None means the same folder as input files.
        threads: int, how many thread to use.
        stepNum: int, step number for folder name.
        upstream: upstream output results, used for pipeline.
        verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
    """

    print("Now, running cell free DNA WGBS data processing pipeline.")

    # set final results
    results = {}

    # set input
    if inputFolder is not None:
        res_inputprocess = inputprocess(inputFolder=inputFolder)
    else:
        res_inputprocess = inputprocess(fqInput1=fastq1, fqInput2=fastq2, verbose=verbose)

    results.update({"inputprocess": res_inputprocess})

    # fastqc
    res_fastqc = fastqc(upstream=res_inputprocess, other_params=fastqcOP, verbose=verbose)
    results.update({"fastqc": res_fastqc})

    # identify and remove adapters,  bismark
    if rmAdapter:
        if idAdapter:
            res_identifyAdapter = identifyAdapter(upstream=res_inputprocess, other_params=idAdOP, verbose=verbose)
            res_adapterremoval = adapterremoval(upstream=res_identifyAdapter, other_params=rmAdOP, verbose=verbose)
            results.update({"identifyAdapter": res_identifyAdapter, "adapterremoval": res_adapterremoval})
        else:
            res_adapterremoval = adapterremoval(
                upstream=res_inputprocess, adapter1=adapter1, adapter2=adapter2, other_params=rmAdOP, verbose=verbose
            )
            results.update({"adapterremoval": res_adapterremoval})

        res_bismark = bismark(upstream=res_adapterremoval, other_params=bismarkOP, verbose=False)
        results.update({"bismark": res_bismark})

    else:
        res_bismark = bismark(upstream=res_inputprocess, other_params=bismarkOP, verbose=False)
        results.update({"bismark": res_bismark})

    # redup and extract methy
    if dudup:
        res_deduplicate = bismark_deduplicate(upstream=res_bismark, other_params=dudupOP, verbose=False)
        res_methyextract = bismark_methylation_extractor(
            upstream=res_deduplicate, other_params=extractMethyOP, verbose=False
        )
        res_bamsort = bamsort(upstream=res_deduplicate, verbose=False)
        results.update(
            {
                "bismark_deduplicate": res_deduplicate,
                "bismark_methylation_extractor": res_methyextract,
                "bamsort": res_bamsort,
            }
        )
    else:
        res_methyextract = bismark_methylation_extractor(
            upstream=res_bismark, other_params=extractMethyOP, verbose=False
        )
        res_bamsort = bamsort(upstream=res_bismark, verbose=False)
        results.update(
            {"bismark_methylation_extractor": res_methyextract, "bamsort": res_bamsort,}
        )

    res_compressMethy = compress_methyl(upstream=res_methyextract, verbose=False)
    res_calMethy = calculate_methyl(upstream=res_compressMethy, bedInput=methyRegion, verbose=False)

    res_bam2bed = bam2bed(upstream=res_bamsort, verbose=False)
    res_fraglenplot = fraglenplot(upstream=res_bam2bed, verbose=False)
    results.update(
        {
            "compress_methyl": res_compressMethy,
            "calculate_methyl": res_calMethy,
            "bam2bed": res_bam2bed,
            "fraglenplot": res_fraglenplot,
        }
    )

    if CNV:
        cnv_readCounter = runCounter(upstream=res_bamsort, filetype=1, verbose=False, stepNum="CNV01")
        cnv_gcCounter = runCounter(filetype=0, upstream=True, verbose=False, stepNum="CNV02")
        cnv_GCCorrect = GCCorrect(
            readupstream=cnv_readCounter, gcupstream=cnv_gcCounter, verbose=False, stepNum="CNV03"
        )
        results.update(
            {"cnvReadCounter": cnv_readCounter, "cnvGCCounter": cnv_gcCounter, "cnvGCCorrect": cnv_GCCorrect,}
        )
    else:
        print("Skip CNV analysis.")

    if fragProfile:
        fp_fragCounter = fpCounter(upstream=res_bam2bed, verbose=False, stepNum="FP02")
        fp_gcCounter = runCounter(filetype=0, binlen=5000000, upstream=True, verbose=False, stepNum="FP01")
        fp_GCCorrect = GCCorrect(
            readupstream=fp_fragCounter, gcupstream=fp_gcCounter, readtype=2, corrkey="-", verbose=False, stepNum="FP03"
        )
        results.update(
            {"fpCounter": fp_fragCounter, "fpGCCounter": fp_gcCounter, "fpGCCorrect": fp_GCCorrect,}
        )
    else:
        print("Skip fragmentation analysis.")

    # set all results
    results = Box(results, frozen_box=True)

    return results
