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
from .Configure import Configure
from .cfDNA_utils import logoPrint


def cfDNAWGS(
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
    bowtie2OP={"-q": True, "-N": 1, "--time": True},
    dudup=True,
    CNV=False,
    verbose=False,
):
    print("Now, running cell free DNA WGS data processing pipeline.")

    # set final results
    results = {}


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
    bismarkOP={"-q": True, "--phred33-quals": True, "--bowtie2": True, "--un": True, },
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
    This function is used for processing paired/single end WGBS data.
    Note: User must set Configure or pipeConfigure before using.
          This function provides basic processing steps, like alignment and bis correction.
          If you want to perform comparison process, please use cfDNAWGBS2.

    cfDNAWGBS(inputFolder=None, fastq1=None, fastq2=None,
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
              verbose=False)
    {P}arameters:
        inputFolder: str, input fastq file folder path. Setting this parameter means disable fastq1 and fastq2.
        fastq1: list, fastq1 files.
        fastq2: list, fastq2 files.
        adapter1: list, adapters for fastq1 files, if idAdapter is True, this paramter will be disabled.
        adapter2: list, adapters for fastq2 files, if idAdapter is True, this paramter will be disabled.
        fastqcOP: Other parameters used for fastqc, please see class "fastqc".
        idAdapter: Ture or False, identify adapters or not. This module is not for single end data.
        idAdOP: Other parameters used for AdapterRemoval, please see class "identifyAdapter".
                If idAdapter is False, this paramter will be disabled.
        rmAdapter: Ture or False, remove adapters or not.
        rmAdOP: Other parameters used for AdapterRemoval, please see class "adapterremoval".
                Default: {"--qualitybase": 33, "--gzip": True}.
                If rmAdapter is False, this paramter will be disabled.
        bismarkOP: Other parameters used for Bismark, please see class "bismark".
                   Default: {"-q": True, "--phred33-quals": True, "--bowtie2": True, "--un": True,}.
        dudup: Ture or False, remove duplicates for bismark results or not.
        dudupOP: Other parameters used for bismark_deduplicate, please see class "bismark_deduplicate".
                 If dudup is False, this paramter will be disabled.
        extractMethyOP: Other parameters used for bismark_methylation_extractor, please see class "bismark_methylation_extractor".
                        Default: {"--no_overlap": True, "--report": True, "--no_header": True, "--gzip": True,
                                  "--bedGraph": True, "--zero_based": True,}
        methyRegion: Bed file contains methylation related regions.
        CNV: Compute basic CNV or not.
        fragProfile: Compute basic fragProfile(long short fragement statistics) or not. This module is not for single end data.
        verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
    """

    logoPrint(mess="WGBS Pipeline")

    print("Now, running cell free DNA WGBS data processing pipeline.")

    # set final results
    results = {}

    # set input
    if inputFolder is not None:
        res_inputprocess = inputprocess(inputFolder=inputFolder)
    else:
        res_inputprocess = inputprocess(
            fqInput1=fastq1, fqInput2=fastq2, verbose=verbose
        )

    results.update({"inputprocess": res_inputprocess})

    # fastqc
    res_fastqc = fastqc(
        upstream=res_inputprocess, other_params=fastqcOP, verbose=verbose
    )
    results.update({"fastqc": res_fastqc})

    # identify, remove adapters and bismark
    if rmAdapter:
        if idAdapter and (Configure.getType() == "paired"):
            res_identifyAdapter = identifyAdapter(
                upstream=res_inputprocess, other_params=idAdOP, verbose=verbose
            )
            res_adapterremoval = adapterremoval(
                upstream=res_identifyAdapter, other_params=rmAdOP, verbose=verbose
            )
            results.update(
                {
                    "identifyAdapter": res_identifyAdapter,
                    "adapterremoval": res_adapterremoval,
                }
            )
        else:
            res_adapterremoval = adapterremoval(
                upstream=res_inputprocess,
                adapter1=adapter1,
                adapter2=adapter2,
                other_params=rmAdOP,
                verbose=verbose,
            )
            results.update({"adapterremoval": res_adapterremoval})

        res_bismark = bismark(
            upstream=res_adapterremoval, other_params=bismarkOP, verbose=verbose
        )
        results.update({"bismark": res_bismark})

    else:
        res_bismark = bismark(
            upstream=res_inputprocess, other_params=bismarkOP, verbose=verbose
        )
        results.update({"bismark": res_bismark})

    # redup and extract methy
    if dudup:
        res_deduplicate = bismark_deduplicate(
            upstream=res_bismark, other_params=dudupOP, verbose=verbose
        )
        res_methyextract = bismark_methylation_extractor(
            upstream=res_deduplicate, other_params=extractMethyOP, verbose=verbose
        )
        res_bamsort = bamsort(upstream=res_deduplicate, verbose=verbose)
        results.update(
            {
                "bismark_deduplicate": res_deduplicate,
                "bismark_methylation_extractor": res_methyextract,
                "bamsort": res_bamsort,
            }
        )
    else:
        res_methyextract = bismark_methylation_extractor(
            upstream=res_bismark, other_params=extractMethyOP, verbose=verbose
        )
        res_bamsort = bamsort(upstream=res_bismark, verbose=verbose)
        results.update(
            {"bismark_methylation_extractor": res_methyextract, "bamsort": res_bamsort, }
        )

    res_compressMethy = compress_methyl(upstream=res_methyextract, verbose=verbose)
    res_calMethy = calculate_methyl(
        upstream=res_compressMethy, bedInput=methyRegion, verbose=verbose
    )

    res_bam2bed = bam2bed(upstream=res_bamsort, verbose=verbose)

    results.update(
        {
            "compress_methyl": res_compressMethy,
            "calculate_methyl": res_calMethy,
            "bam2bed": res_bam2bed,
        }
    )

    if Configure.getType() == "paired":
        res_fraglenplot = fraglenplot(upstream=res_bam2bed, verbose=verbose)
        results.update(
            {"fraglenplot": res_fraglenplot, }
        )

    if CNV:
        cnv_readCounter = runCounter(
            upstream=res_bamsort, filetype=1, verbose=verbose, stepNum="CNV01"
        )
        cnv_gcCounter = runCounter(
            filetype=0, upstream=True, verbose=verbose, stepNum="CNV02"
        )
        cnv_GCCorrect = GCCorrect(
            readupstream=cnv_readCounter,
            gcupstream=cnv_gcCounter,
            verbose=verbose,
            stepNum="CNV03",
        )
        results.update(
            {
                "cnvReadCounter": cnv_readCounter,
                "cnvGCCounter": cnv_gcCounter,
                "cnvGCCorrect": cnv_GCCorrect,
            }
        )
    else:
        print("Skip CNV analysis.")

    if fragProfile and (Configure.getType() == "paired"):
        fp_fragCounter = fpCounter(
            upstream=res_bam2bed, verbose=verbose, stepNum="FP02"
        )
        fp_gcCounter = runCounter(
            filetype=0, binlen=5000000, upstream=True, verbose=verbose, stepNum="FP01"
        )
        fp_GCCorrect = GCCorrect(
            readupstream=fp_fragCounter,
            gcupstream=fp_gcCounter,
            readtype=2,
            corrkey="-",
            verbose=verbose,
            stepNum="FP03",
        )
        results.update(
            {
                "fpCounter": fp_fragCounter,
                "fpGCCounter": fp_gcCounter,
                "fpGCCorrect": fp_GCCorrect,
            }
        )
    else:
        print("Skip fragmentation analysis.")

    # set all results
    results = Box(results, frozen_box=True)

    return results
