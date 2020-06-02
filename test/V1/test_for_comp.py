# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:38:07 2019

@author: zhang
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

# case processing
switchConfigure("case")
res_case1 = inputprocess(
    inputFolder=r"/home/zhangwei/pipeline-WGBS-cc/raw/case_small")
res_case2 = fastqc(upstream=res_case1)
res_case3 = identifyAdapter(upstream=res_case1)
res_case4 = adapterremoval(upstream=res_case3)
res_case5 = bismark(upstream=res_case4)
res_case6 = bismark_deduplicate(upstream=res_case5)
res_case7 = bismark_methylation_extractor(upstream=res_case6)
res_case8 = compress_methyl(upstream=res_case7)
res_case9 = calculate_methyl(upstream=res_case8)
res_case10 = bamsort(upstream=res_case6)
res_case11 = bam2bed(upstream=res_case10)
res_case12 = fraglenplot(upstream=res_case11)
res_case13 = readCount(upstream=res_case10)

# ctrl processing
switchConfigure("ctrl")
res_ctrl1 = inputprocess(
    inputFolder=r"/home/zhangwei/pipeline-WGBS-cc/raw/ctrl_small")
res_ctrl2 = fastqc(upstream=res_ctrl1)
res_ctrl3 = identifyAdapter(upstream=res_ctrl1)
res_ctrl4 = adapterremoval(upstream=res_ctrl3)
res_ctrl5 = bismark(upstream=res_ctrl4)
res_ctrl6 = bismark_deduplicate(upstream=res_ctrl5)
res_ctrl7 = bismark_methylation_extractor(upstream=res_ctrl6)
res_ctrl8 = compress_methyl(upstream=res_ctrl7)
res_ctrl9 = calculate_methyl(upstream=res_ctrl8)
res_ctrl10 = bamsort(upstream=res_ctrl6)
res_ctrl11 = bam2bed(upstream=res_ctrl10)
res_ctrl12 = fraglenplot(upstream=res_ctrl11)
res_ctrl13 = readCount(upstream=res_ctrl10)

res1 = computeOCF(caseupstream=res_case11, ctrlupstream=res_ctrl11)
'''
res2 = computeCNV(caseupstream=res_case13, ctrlupstream=res_ctrl13)
'''
resw = computeCNV(ctrlreadInput=['/home/zhangwei/CNV/CUHK/CTR101_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/CTR103_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/CTR104_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/CTR106_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/CTR107_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/CTR108_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/CTR110_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/CTR111_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/CTR113_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/CTR114_sorted.readcount.wig'],
    casereadInput=['/home/zhangwei/CNV/CUHK/HOT151_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/HOT156_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/HOT159_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/HOT162_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/HOT164_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/HOT167_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/HOT170_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/HOT172_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/HOT192_sorted.readcount.wig',
    '/home/zhangwei/CNV/CUHK/HOT197_sorted.readcount.wig'],
    casegcInput='/home/zhangwei/CNV/CUHK/hg19.gc.wig',
    ctrlgcInput='/home/zhangwei/CNV/CUHK/hg19.gc.wig',
    outputdir='/home/zhangwei/CNV/test',
    cytoBandInput='/home/zhangwei/Genome/hg19_bismark/cytoBand.txt',
    )
res3 = fraglenplot_comp(caseupstream=res_case11, ctrlupstream=res_ctrl11)
res4 = PCAplot(caseupstream=res_case9, ctrlupstream=res_ctrl9)
resp = runDeconCCN(mixInput = ["/home/zhangwei/test/DeconCCN/demo.txt",
                               "/home/zhangwei/test/DeconCCN/demo2.txt",
                               "/home/zhangwei/test/DeconCCN/demo3.txt",
                               "/home/zhangwei/test/DeconCCN/demo4.txt",
                               "/home/zhangwei/test/DeconCCN/demo5.txt",
                               "/home/zhangwei/test/DeconCCN/demo6.txt",
                               "/home/zhangwei/test/DeconCCN/demo7.txt",
                              ], 
                   refInput = "/home/zhangwei/test/DeconCCN/ref.txt", 
                   outputdir = "/home/zhangwei/test/DeconCCN/result/",
)
res6 = fragprofplot(
    caseupstream = res_case11, 
    ctrlupstream = res_ctrl11, 
    chromsizeInput = r"/home/zhangwei/test/hg19.chrom.sizes",
)

rep = report_generator_comp(
    case_fastqcRes = res_case2,
    case_identifyAdapterRes = res_case3,
    case_bismarkRes = res_case5,
    case_deduplicateRes = res_case6,
    case_rmduplicateRes = None,
    case_fraglenplotRes = res_case12,
    case_DeconCCNRes = resp,
    ctrl_fastqcRes = res_ctrl2,
    ctrl_identifyAdapterRes = res_ctrl3,
    ctrl_bismarkRes = res_ctrl5,
    ctrl_deduplicateRes = res_ctrl6,
    ctrl_rmduplicateRes = None,
    ctrl_fraglenplotRes = res_ctrl12,
    ctrl_DeconCCNRes = resp,
    OCFRes = res1,
    CNVRes = resw,
    fraglenplotcompRes = res3,
    PCARes = res4,
    fragprofplotRes = res6,
    outputdir = None,
    label = None,
)