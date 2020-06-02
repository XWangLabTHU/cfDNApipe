# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:38:07 2019

@author: zhang
"""

import glob
from cfDNApipe import *

# set global configure
pipeConfigure2(
    threads=60,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bowtie2",
    outdir="/home/zhangwei/pipeline-WGS-comp",
    data="WGS",
    type="paired",
    JavaMem="10G",
    case="normal",
    ctrl="cancer",
    build=True,
)

# ctrl processing
switchConfigure("normal")
res_ctrl6 = bamsort(
    bamInput=glob.glob("/home/zhangwei/pipeline-WGS-comp/bams/case_raw/*.bam"),
    upstream=True,
)
res_ctrl7 = rmduplicate(upstream=res_ctrl6, verbose=False)
res_ctrl8 = addRG(upstream=res_ctrl7, verbose=False)
res_ctrl9 = BaseRecalibrator(
    upstream=res_ctrl8, knownSitesDir=r"/opt/tsinghua/cfDNApipeTest/file/vcf", verbose=False
)
res_ctrl10 = BQSR(upstream=res_ctrl9, verbose=False)
res_ctrl11 = mutect2n(upstream=res_ctrl10, verbose=False)
res_ctrl12 = dbimport(upstream=res_ctrl11, verbose=False)
res_ctrl13 = createPON(upstream=res_ctrl12, verbose=False)
res_ctrl14 = runCounter(upstream=res_ctrl6, verbose=False)
res_ctrl15 = qualimap(upstream=res_ctrl6, verbose=False)
res_ctrl16 = bam2bed(upstream=res_ctrl7, verbose=False)

# ctrl processing
switchConfigure("cancer")
res_case6 = bamsort(
    bamInput=glob.glob("/home/zhangwei/pipeline-WGS-comp/bams/ctrl_raw/*.bam"),
    upstream=True,
)
res_case7 = rmduplicate(upstream=res_case6, verbose=False)
res_case8 = addRG(upstream=res_case7, verbose=False)
res_case9 = BaseRecalibrator(
    upstream=res_case8, knownSitesDir=r"/opt/tsinghua/cfDNApipeTest/file/vcf", verbose=False
)
res_case10 = BQSR(upstream=res_case9, verbose=False)
res_case11 = getPileup(
    upstream=res_case10,
    biallelicvcfInput="/opt/tsinghua/cfDNApipeTest/file/small_exac_common_3_hg19.SNP_biallelic.vcf", verbose=False
)
res_case12 = contamination(upstream=res_case11, verbose=False)
res_case13 = runCounter(upstream=res_case6, verbose=False)
res_case14 = qualimap(upstream=res_case6, verbose=False)
res_case16 = bam2bed(upstream=res_case7, verbose=False)

# snv / InDel
res1 = mutect2t(
    caseupstream=res_case12,
    ctrlupstream=res_ctrl13,
    vcfInput="/opt/tsinghua/cfDNApipeTest/file/af-only-gnomad.raw.sites.hg19.vcf.gz",
)
res2 = filterMutectCalls(upstream=res1)
res3 = gatherVCF(upstream=res2)
res4 = annovar(
    upstream=res3,
    plInput="/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/table_annovar.pl",
    dbdir="/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/humandb/",
)
res5 = annovarStat(upstream=res4,)


# cnv
res6 = cnvbatch(
    caseupstream=res_case7,
    ctrlupstream=res_ctrl7,
    access="/opt/tsinghua/cfDNApipeTest/file/CNVkit/access-5kb-mappable.hg19.bed",
    annotate="/opt/tsinghua/cfDNApipeTest/file/CNVkit/refFlat_hg19.txt",
    stepNum=18,
)
res7 = cnvPlot(upstream=res6, verbose=True)
res8 = cnvTable(upstream=res6)
res9 = cnvHeatmap(upstream=res6)

# OCF
resO = computeOCF(
    caseupstream=res_case16,
    ctrlupstream=res_ctrl16,
    labelInput=["HCC", "CTR"],
    refRegInput=Configure2.getConfig("OCF"),
    verbose=False,
    stepNum=1,
)

res1 = OCFplot(upstream=resO, labelInput=["HCC", "CTR"],)
