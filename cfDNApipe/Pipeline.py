# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:26:39 2019

@author: zhang Wei
"""

from box import Box
import time
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
from .Fun_bowtie2 import bowtie2
from .Fun_rmDuplicate import rmduplicate
from .Fun_cnvbatch import cnvbatch
from .Fun_cnvPlot import cnvPlot
from .Fun_cnvTable import cnvTable
from .Fun_cnvHeatmap import cnvHeatmap
from .Fun_bamcount import bamCounter
from .Configure import Configure
from .cfDNA_utils import logoPrint
from .Fun_OCF import computeOCF
from .Fun_OCFplot import OCFplot
from .Fun_CNV import computeCNV
from .Fun_qualimap import qualimap
from .Fun_fragLencomp import fraglenplot_comp
from .Fun_fpplot import fragprofplot
from .Fun_DeconCCN import runDeconCCN
from .Fun_PCA import PCAplot
from .Fun_DMR import computeDMR
from .report_generator import report_generator
from .Configure import *
from .Configure2 import *


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
    armCNV=False,
    fragProfile=False,
    report=False,
    verbose=False,
    box=True,
):
    """
    This function is used for processing paired/single end WGBS data.
    Note: User must set Configure or pipeConfigure before using.
          This function provides basic processing steps, like alignment and bis correction.
          If you want to perform comparison process, please use cfDNAWGBS2.

    cfDNAWGS(inputFolder=None,
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
             armCNV=False,
             fragProfile=False,
             verbose=False,
             box=True)
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
        bowtie2OP: Other parameters used for Bowtie2, please see class "bowtie2".
                   Default: {"-q": True, "-N": 1, "--time": True}.
        dudup: Ture or False, remove duplicates for bowtie2 results or not.
        CNV: Compute basic CNV or not.
        armCNV:  Compute arm level CNV or not.
        fragProfile: Compute basic fragProfile(long short fragement statistics) or not. This module is not for single end data.
        report: Generate user report or not.
        verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        box: output will be a box class. that mean the the results can be specified using '.',
             otherwise, the results will be saved as dict.
    """

    logoPrint(mess="WGS Pipeline")

    time.sleep(3)

    print("")
    print("Now, running cell free DNA WGS data processing pipeline.")
    print("")

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

    # identify, remove adapters and bowtie2
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

        res_bowtie2 = bowtie2(
            upstream=res_adapterremoval, other_params=bowtie2OP, verbose=verbose
        )
        results.update({"bowtie2": res_bowtie2})

    else:
        res_bowtie2 = bowtie2(
            upstream=res_inputprocess, other_params=bowtie2OP, verbose=verbose
        )
        results.update({"bismark": res_bowtie2})

    # sort bam files
    res_bamsort = bamsort(upstream=res_bowtie2, verbose=verbose)
    res_qualimap = qualimap(upstream=res_bamsort, verbose=verbose)
    results.update({"bamsort": res_bamsort, "qualimap": res_qualimap})

    # remove duplicates
    if dudup:
        res_rmduplicate = rmduplicate(upstream=res_bamsort, verbose=verbose)

        # bam2bed
        res_bam2bed = bam2bed(upstream=res_rmduplicate, verbose=verbose)
        results.update({"rmduplicate": res_rmduplicate, "bam2bed": res_bam2bed})
    else:
        # bam2bed
        res_bam2bed = bam2bed(upstream=res_bamsort)
        results.update({"bam2bed": res_bam2bed})

    # fraglenplot
    if Configure.getType() == "paired":
        res_fraglenplot = fraglenplot(upstream=res_bam2bed, verbose=verbose)
        results.update({"fraglenplot": res_fraglenplot})

    # cnv
    if CNV:
        res_cnvbatch = cnvbatch(
            caseupstream=res_rmduplicate,
            access=Configure.getConfig("access-5kb-mappable"),
            annotate=Configure.getConfig("refFlat"),
            verbose=verbose,
            stepNum="CNV01",
        )
        res_cnvPlot = cnvPlot(upstream=res_cnvbatch, verbose=verbose, stepNum="CNV02",)
        res_cnvTable = cnvTable(
            upstream=res_cnvbatch, verbose=verbose, stepNum="CNV03",
        )
        res_cnvHeatmap = cnvHeatmap(
            upstream=res_cnvbatch, verbose=verbose, stepNum="CNV04",
        )
        results.update(
            {
                "cnvPlot": res_cnvPlot,
                "cnvTable": res_cnvTable,
                "cnvHeatmap": res_cnvHeatmap,
            }
        )
    else:
        print("Skip CNV analysis.")

    # arm level CNV
    if armCNV:
        res_bamCounter = bamCounter(
            upstream=res_rmduplicate, verbose=verbose, stepNum="ARMCNV01"
        )
        cnv_gcCounter = runCounter(
            filetype=0, upstream=True, verbose=verbose, stepNum="ARMCNV02"
        )
        cnv_GCCorrect = GCCorrect(
            readupstream=res_bamCounter,
            gcupstream=cnv_gcCounter,
            verbose=verbose,
            stepNum="ARMCNV03",
        )
        results.update(
            {
                "cnvbamCounter": res_bamCounter,
                "cnvGCCounter": cnv_gcCounter,
                "cnvGCCorrect": cnv_GCCorrect,
            }
        )
    else:
        print("Skip arm level CNV analysis.")

    # fragProfile
    if fragProfile and (Configure.getType() == "paired"):
        fp_fragCounter = fpCounter(
            upstream=res_bam2bed, verbose=verbose, stepNum="FP01", processtype=1
        )
        fp_gcCounter = runCounter(
            filetype=0, binlen=5000000, upstream=True, verbose=verbose, stepNum="FP02"
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

    # report
    if report:
        if "fastqc" in results:
            fastqcRes = results["fastqc"]
        else:
            fastqcRes = None
        if "identifyAdapter" in results:
            identifyAdapterRes = results["identifyAdapter"]
        else:
            identifyAdapterRes = None
        if "bismark" in results:
            bismarkRes = results["bismark"]
        else:
            bismarkRes = None
        if "qualimap" in results:
            qualimapRes = results["qualimap"]
        else:
            qualimapRes = None
        if "rmduplicate" in results:
            rmduplicateRes = results["rmduplicate"]
        else:
            rmduplicateRes = None
        if "fraglenplot" in results:
            fraglenplotRes = results["fraglenplot"]
        else:
            fraglenplotRes = None
        if "cnvPlot" in results:
            CNVplotRes = results["cnvPlot"]
        else:
            CNVplotRes = None
        if "cnvHeatmap" in results:
            CNVheatmapRes = results["cnvHeatmap"]
        else:
            CNVheatmapRes = None
        if "cnvGCCorrect" in results:
            CNV_GCcorrectRes = results["cnvGCCorrect"]
        else:
            CNV_GCcorrectRes = None
        if "fpGCCorrect" in results:
            fragprof_GCcorrectRes = results["fpGCCorrect"]
        else:
            fragprof_GCcorrectRes = None
        report_generator(
            fastqcRes=fastqcRes,
            identifyAdapterRes=identifyAdapterRes,
            bismarkRes=bismarkRes,
            qualimapRes=qualimapRes,
            deduplicateRes=None,
            rmduplicateRes=rmduplicateRes,
            fraglenplotRes=fraglenplotRes,
            CNVplotRes=CNVplotRes,
            CNVheatmapRes=CNVheatmapRes,
            CNV_GCcorrectRes=CNV_GCcorrectRes,
            fragprof_GCcorrectRes=fragprof_GCcorrectRes,
            outputdir=None,
        )
    else:
        print("Skip report generation.")

    # set all results
    if box:
        results = Box(results, frozen_box=True)

    return results


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
    bismarkOP={"-q": True, "--phred33-quals": True, "--bowtie2": True, "--un": True},
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
    armCNV=False,
    CNV=False,
    fragProfile=False,
    deconvolution=False,
    verbose=False,
    box=True,
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
              deconvolution=False,
              verbose=False,
              box=True)
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
        armCNV:  Compute arm level CNV or not.
        fragProfile: Compute basic fragProfile(long short fragement statistics) or not. This module is not for single end data.
        deconvolution: Compute tissue proportion for each sample or not.
        verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        box: output will be a box class. that mean the the results can be specified using '.',
             otherwise, the results will be saved as dict.
    """

    logoPrint(mess="WGBS Pipeline")

    time.sleep(3)

    print("")
    print("Now, running cell free DNA WGBS data processing pipeline.")
    print("")

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
        res_qualimap = qualimap(upstream=res_bamsort, verbose=verbose)
        results.update(
            {
                "bismark_deduplicate": res_deduplicate,
                "bismark_methylation_extractor": res_methyextract,
                "bamsort": res_bamsort,
                "qualimap": res_qualimap,
            }
        )
    else:
        res_methyextract = bismark_methylation_extractor(
            upstream=res_bismark, other_params=extractMethyOP, verbose=verbose
        )
        res_bamsort = bamsort(upstream=res_bismark, verbose=verbose)
        res_qualimap = qualimap(upstream=res_bamsort, verbose=verbose)
        results.update(
            {
                "bismark_methylation_extractor": res_methyextract,
                "bamsort": res_bamsort,
                "qualimap": res_qualimap,
            }
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
        results.update({"fraglenplot": res_fraglenplot})

    # cnv
    if CNV:
        res_cnvbatch = cnvbatch(
            caseupstream=res_bamsort,
            access=Configure.getConfig("access-5kb-mappable"),
            annotate=Configure.getConfig("refFlat"),
            verbose=verbose,
            stepNum="CNV01",
        )
        res_cnvPlot = cnvPlot(upstream=res_cnvbatch, verbose=verbose, stepNum="CNV02",)
        res_cnvTable = cnvTable(
            upstream=res_cnvbatch, verbose=verbose, stepNum="CNV03",
        )
        res_cnvHeatmap = cnvHeatmap(
            upstream=res_cnvbatch, verbose=verbose, stepNum="CNV04",
        )
        results.update(
            {
                "cnvPlot": res_cnvPlot,
                "cnvTable": res_cnvTable,
                "cnvHeatmap": res_cnvHeatmap,
            }
        )
    else:
        print("Skip CNV analysis.")

    if armCNV:
        res_bamCounter = bamCounter(
            upstream=res_bamsort, verbose=verbose, stepNum="ARMCNV01"
        )
        cnv_gcCounter = runCounter(
            filetype=0, upstream=True, verbose=verbose, stepNum="ARMCNV02"
        )
        cnv_GCCorrect = GCCorrect(
            readupstream=res_bamCounter,
            gcupstream=cnv_gcCounter,
            verbose=verbose,
            stepNum="ARMCNV03",
        )
        results.update(
            {
                "cnvbamCounter": res_bamCounter,
                "cnvGCCounter": cnv_gcCounter,
                "cnvGCCorrect": cnv_GCCorrect,
            }
        )
    else:
        print("Skip CNV analysis.")

    if fragProfile and (Configure.getType() == "paired"):
        fp_fragCounter = fpCounter(
            upstream=res_bam2bed, verbose=verbose, stepNum="FP01", processtype=1
        )
        fp_gcCounter = runCounter(
            filetype=0, binlen=5000000, upstream=True, verbose=verbose, stepNum="FP02"
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

    if deconvolution:
        res_runDeconCCN = runDeconCCN(upstream=res_calMethy)
        results.update({"runDeconCCN": res_runDeconCCN})

    # set all results
    if box:
        results = Box(results, frozen_box=True)

    return results


def cfDNAWGS2(
    caseFolder=None,
    ctrlFolder=None,
    caseName=None,
    ctrlName=None,
    caseFq1=None,
    caseFq2=None,
    ctrlFq1=None,
    ctrlFq2=None,
    caseAdapter1=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
    caseAdapter2=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"],
    ctrlAdapter1=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
    ctrlAdapter2=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"],
    fastqcOP=None,
    idAdapter=True,
    idAdOP=None,
    rmAdapter=True,
    rmAdOP={"--qualitybase": 33, "--gzip": True},
    bowtie2OP={"-q": True, "-N": 1, "--time": True},
    dudup=True,
    CNV=False,
    armCNV=False,
    fragProfile=False,
    OCF=False,
    verbose=False,
    box=True,
):
    """
    This function is used for case/control analysis of paired/single end WGS data.
    Note: User must set pipeConfigure2 before using.

    cfDNAWGS2(caseFolder=None,
              ctrlFolder=None,
              caseName=None,
              ctrlName=None,
              caseFq1=None,
              caseFq2=None,
              ctrlFq1=None,
              ctrlFq2=None,
              caseAdapter1=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
              caseAdapter2=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"],
              ctrlAdapter1=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
              ctrlAdapter2=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"],
              fastqcOP=None,
              idAdapter=True,
              idAdOP=None,
              rmAdapter=True,
              rmAdOP={"--qualitybase": 33, "--gzip": True},
              bowtie2OP={"-q": True, "-N": 1, "--time": True},
              dudup=True,
              CNV=False,
              armCNV=False,
              fragProfile=False,
              verbose=False,
              box=True)
    {P}arameters:
        caseFolder: str, input case fastq file folder path. Setting this parameter means disable fastq1 and fastq2.
        ctrlFolder: str, input control fastq file folder path. Setting this parameter means disable fastq1 and fastq2.
        caseFq1: list, case fastq1 files.
        caseFq2: list, case fastq2 files.
        ctrlFq1: list, control fastq1 files.
        ctrlFq2: list, control fastq2 files.
        caseAdapter1: list, adapters for case fastq1 files, if idAdapter is True, this paramter will be disabled.
        caseAdapter2: list, adapters control for case fastq2 files, if idAdapter is True, this paramter will be disabled.
        ctrlAdapter1: list, adapters control for fastq1 files, if idAdapter is True, this paramter will be disabled.
        ctrlAdapter2: list, adapters for fastq2 files, if idAdapter is True, this paramter will be disabled.
        fastqcOP: Other parameters used for fastqc, please see class "fastqc".
        idAdapter: Ture or False, identify adapters or not. This module is not for single end data.
        idAdOP: Other parameters used for AdapterRemoval, please see class "identifyAdapter".
                If idAdapter is False, this paramter will be disabled.
        rmAdapter: Ture or False, remove adapters or not.
        rmAdOP: Other parameters used for AdapterRemoval, please see class "adapterremoval".
                Default: {"--qualitybase": 33, "--gzip": True}.
                If rmAdapter is False, this paramter will be disabled.
        bowtie2OP: Other parameters used for Bowtie2, please see class "bowtie2".
                   Default: {"-q": True, "-N": 1, "--time": True}.
        dudup: Ture or False, remove duplicates for bowtie2 results or not.
        CNV: Compute basic CNV or not.
        armCNV:  Compute arm level CNV or not.
        fragProfile: Compute basic fragProfile(long short fragement statistics) or not. This module is not for single end data.
        OCF: Compute OCF or not.
        verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        box: output will be a box class. that mean the the results can be specified using '.',
             otherwise, the results will be saved as dict.
    """

    switchConfigure(caseName)
    mess = "Now, Start processing " + caseName + "......"
    print(mess)
    caseOut_dict = cfDNAWGS(
        inputFolder=caseFolder,
        fastq1=caseFq1,
        fastq2=caseFq2,
        adapter1=caseAdapter1,
        adapter2=caseAdapter2,
        fastqcOP=fastqcOP,
        idAdapter=idAdapter,
        idAdOP=idAdOP,
        rmAdapter=rmAdapter,
        rmAdOP=rmAdOP,
        bowtie2OP=bowtie2OP,
        dudup=dudup,
        CNV=CNV,
        armCNV=armCNV,
        fragProfile=fragProfile,
        verbose=verbose,
        box=False,
    )

    caseOut = Box(caseOut_dict, frozen_box=True)

    switchConfigure(ctrlName)
    mess = "Now, Start processing " + ctrlName + "......"
    print(mess)
    ctrlOut_dict = cfDNAWGS(
        inputFolder=ctrlFolder,
        fastq1=ctrlFq1,
        fastq2=ctrlFq2,
        adapter1=ctrlAdapter1,
        adapter2=ctrlAdapter2,
        fastqcOP=fastqcOP,
        idAdapter=idAdapter,
        idAdOP=idAdOP,
        rmAdapter=rmAdapter,
        rmAdOP=rmAdOP,
        bowtie2OP=bowtie2OP,
        dudup=dudup,
        CNV=CNV,
        armCNV=armCNV,
        fragProfile=fragProfile,
        verbose=verbose,
        box=False,
    )

    ctrlOut = Box(ctrlOut_dict, frozen_box=True)

    switchConfigure(caseName)

    # set comparison results
    results = {}

    # fragment length comparision
    if Configure.getType() == "paired":
        res_fraglenplot_comp = fraglenplot_comp(
            caseupstream=caseOut.bam2bed, ctrlupstream=ctrlOut.bam2bed, verbose=verbose
        )
        results.update({"fraglenplot_comp": res_fraglenplot_comp})

        if fragProfile:
            res_fragprofplot = fragprofplot(
                caseupstream=caseOut.fpGCCorrect,
                ctrlupstream=ctrlOut.fpGCCorrect,
                stepNum="FP04",
            )
            results.update({"fragprofplot": res_fragprofplot})

        if OCF:
            res_computeOCF = computeOCF(
                caseupstream=caseOut.bam2bed,
                ctrlupstream=ctrlOut.bam2bed,
                verbose=verbose,
            )
            res_OCFplot = OCFplot(upstream=res_computeOCF, verbose=verbose)
            results.update({"computeOCF": res_computeOCF, "OCFplot": res_OCFplot})

    # ARM-level CNV plot
    if armCNV:
        res_computeCNV = computeCNV(
            caseupstream=caseOut.cnvGCCorrect,
            ctrlupstream=ctrlOut.cnvGCCorrect,
            stepNum="ARMCNV",
            verbose=verbose,
        )

        results.update({"computeCNV": res_computeCNV})

    # set all results
    caseOut_dict.update(results)
    fi_caseOut = Box(caseOut_dict, frozen_box=True)
    fi_ctrlOut = Box(ctrlOut_dict, frozen_box=True)

    return fi_caseOut, fi_ctrlOut


def cfDNAWGBS2(
    caseFolder=None,
    ctrlFolder=None,
    caseName=None,
    ctrlName=None,
    caseFq1=None,
    caseFq2=None,
    ctrlFq1=None,
    ctrlFq2=None,
    caseAdapter1=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
    caseAdapter2=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"],
    ctrlAdapter1=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
    ctrlAdapter2=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"],
    fastqcOP=None,
    idAdapter=True,
    idAdOP=None,
    rmAdapter=True,
    rmAdOP={"--qualitybase": 33, "--gzip": True},
    bismarkOP={"-q": True, "--phred33-quals": True, "--bowtie2": True, "--un": True},
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
    armCNV=False,
    CNV=False,
    fragProfile=False,
    deconvolution=False,
    OCF=False,
    verbose=False,
    box=True,
):
    """
    This function is used for case/control analysis of paired/single end WGBS data.
    Note: User must set pipeConfigure2 before using.

    cfDNAWGBS2(caseFolder=None,
               ctrlFolder=None,
               caseName=None,
               ctrlName=None,
               caseFq1=None,
               caseFq2=None,
               ctrlFq1=None,
               ctrlFq2=None,
               caseAdapter1=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
               caseAdapter2=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"],
               ctrlAdapter1=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],
               ctrlAdapter2=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"],
               fastqcOP=None,
               idAdapter=True,
               idAdOP=None,
               rmAdapter=True,
               rmAdOP={"--qualitybase": 33, "--gzip": True},
               bismarkOP={"-q": True, "--phred33-quals": True, "--bowtie2": True, "--un": True},
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
               armCNV=False,
               CNV=False,
               fragProfile=False,
               OCF=False,
               verbose=False,
               box=True)
    {P}arameters:
        caseFolder: str, input case fastq file folder path. Setting this parameter means disable fastq1 and fastq2.
        ctrlFolder: str, input control fastq file folder path. Setting this parameter means disable fastq1 and fastq2.
        caseFq1: list, case fastq1 files.
        caseFq2: list, case fastq2 files.
        ctrlFq1: list, control fastq1 files.
        ctrlFq2: list, control fastq2 files.
        caseAdapter1: list, adapters for case fastq1 files, if idAdapter is True, this paramter will be disabled.
        caseAdapter2: list, adapters control for case fastq2 files, if idAdapter is True, this paramter will be disabled.
        ctrlAdapter1: list, adapters control for fastq1 files, if idAdapter is True, this paramter will be disabled.
        ctrlAdapter2: list, adapters for fastq2 files, if idAdapter is True, this paramter will be disabled.
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
        armCNV:  Compute arm level CNV or not.
        CNV: Compute basic CNV or not.
        fragProfile: Compute basic fragProfile(long short fragement statistics) or not. This module is not for single end data.
        deconvolution: Compute tissue proportion for each sample or not.
        OCF: Compute OCF or not.
        verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        box: output will be a box class. that mean the the results can be specified using '.',
             otherwise, the results will be saved as dict.
    """

    switchConfigure(caseName)
    mess = "Now, Start processing " + caseName + "......"
    print(mess)

    caseOut_dict = cfDNAWGBS(
        inputFolder=caseFolder,
        fastq1=caseFq1,
        fastq2=caseFq2,
        adapter1=caseAdapter1,
        adapter2=caseAdapter2,
        fastqcOP=fastqcOP,
        idAdapter=idAdapter,
        idAdOP=idAdOP,
        rmAdapter=rmAdapter,
        rmAdOP=rmAdOP,
        bismarkOP=bismarkOP,
        dudup=dudup,
        dudupOP=dudupOP,
        extractMethyOP=extractMethyOP,
        methyRegion=methyRegion,
        armCNV=armCNV,
        CNV=CNV,
        fragProfile=fragProfile,
        deconvolution=deconvolution,
        verbose=verbose,
        box=False
    )

    caseOut = Box(caseOut_dict, frozen_box=True)

    switchConfigure(ctrlName)
    mess = "Now, Start processing " + ctrlName + "......"
    print(mess)

    ctrlOut_dict = cfDNAWGBS(
        inputFolder=ctrlFolder,
        fastq1=ctrlFq1,
        fastq2=ctrlFq2,
        adapter1=ctrlAdapter1,
        adapter2=ctrlAdapter2,
        fastqcOP=fastqcOP,
        idAdapter=idAdapter,
        idAdOP=idAdOP,
        rmAdapter=rmAdapter,
        rmAdOP=rmAdOP,
        bismarkOP=bismarkOP,
        dudup=dudup,
        dudupOP=dudupOP,
        extractMethyOP=extractMethyOP,
        methyRegion=methyRegion,
        armCNV=armCNV,
        CNV=CNV,
        fragProfile=fragProfile,
        deconvolution=deconvolution,
        verbose=verbose,
        box=False
    )

    ctrlOut = Box(ctrlOut_dict, frozen_box=True)

    switchConfigure(caseName)

    # set comparison results
    results = {}

    # # methylation PCA plot and DMR
    # res_PCA = PCAplot(
    #     caseupstream=caseOut.calculate_methyl, ctrlupstream=ctrlOut.calculate_methyl
    # )

    # res_DMR = computeDMR(
    #     caseupstream=caseOut.calculate_methyl, ctrlupstream=ctrlOut.calculate_methyl
    # )

    # results.update({"PCA": res_PCA, "DMR": res_DMR})

    # fragment length comparision
    if Configure.getType() == "paired":
        res_fraglenplot_comp = fraglenplot_comp(
            caseupstream=caseOut.bam2bed, ctrlupstream=ctrlOut.bam2bed, verbose=verbose
        )
        results.update({"fraglenplot_comp": res_fraglenplot_comp})

        if fragProfile:
            res_fragprofplot = fragprofplot(
                caseupstream=caseOut.fpGCCorrect,
                ctrlupstream=ctrlOut.fpGCCorrect,
                stepNum="FP04",
            )
            results.update({"fragprofplot": res_fragprofplot})

        if OCF:
            res_computeOCF = computeOCF(
                caseupstream=caseOut.bam2bed,
                ctrlupstream=ctrlOut.bam2bed,
                verbose=verbose,
            )
            res_OCFplot = OCFplot(upstream=res_computeOCF, verbose=verbose)
            results.update({"computeOCF": res_computeOCF, "OCFplot": res_OCFplot})

    # ARM-level CNV plot
    if armCNV:
        res_computeCNV = computeCNV(
            caseupstream=caseOut.cnvGCCorrect,
            ctrlupstream=ctrlOut.cnvGCCorrect,
            stepNum="ARMCNV",
            verbose=verbose,
        )

        results.update({"computeCNV": res_computeCNV})

    results.update({"computeCNV": res_computeCNV})

    # set all results
    caseOut_dict.update(results)
    fi_caseOut = Box(caseOut_dict, frozen_box=True)
    fi_ctrlOut = Box(ctrlOut_dict, frozen_box=True)

    return fi_caseOut, fi_ctrlOut
