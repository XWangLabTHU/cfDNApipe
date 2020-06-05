# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:38:07 2019

@author: zhang
"""

# WGS paired end comparison

from cfDNApipe import *

# set global configure
pipeConfigure2(
    threads=60,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bowtie2",
    outdir="/home/zhangwei/pipeline-WGS-comp",
    data="WGS",
    type="paired",
    case="cancer",
    ctrl="normal",
    JavaMem="8G",
    build=True,
)

# ctrl processing
switchConfigure("normal")
res_ctrl1 = inputprocess(inputFolder=r"/home/zhangwei/pipeline-WGS-comp/fastqs/ctrl")
res_ctrl2 = fastqc(upstream=res_ctrl1, verbose=False)
res_ctrl3 = identifyAdapter(upstream=res_ctrl1, stepNum=3, verbose=False)
res_ctrl4 = adapterremoval(upstream=res_ctrl3, verbose=False)
res_ctrl5 = bowtie2(upstream=res_ctrl4, verbose=False)
res_ctrl6 = bamsort(upstream=res_ctrl5, verbose=False)
res_ctrl7 = rmduplicate(upstream=res_ctrl6, verbose=False)
res_ctrl8 = addRG(upstream=res_ctrl7, verbose=False)
res_ctrl9 = BaseRecalibrator(
    upstream=res_ctrl8,
    knownSitesDir=r"/opt/tsinghua/cfDNApipeTest/file/vcf",
    verbose=False,
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
res_case1 = inputprocess(inputFolder=r"/home/zhangwei/pipeline-WGS-comp/fastqs/case")
res_case2 = fastqc(upstream=res_case1, verbose=False)
res_case3 = identifyAdapter(upstream=res_case1, stepNum=3, verbose=False)
res_case4 = adapterremoval(upstream=res_case3, verbose=False)
res_case5 = bowtie2(upstream=res_case4, verbose=False)
res_case6 = bamsort(upstream=res_case5, verbose=False)
res_case7 = rmduplicate(upstream=res_case6, verbose=False)
res_case8 = addRG(upstream=res_case7, verbose=False)
res_case9 = BaseRecalibrator(
    upstream=res_case8,
    knownSitesDir=r"/opt/tsinghua/cfDNApipeTest/file/vcf",
    verbose=False,
)
res_case10 = BQSR(upstream=res_case9, verbose=False)
res_case11 = getPileup(
    upstream=res_case10,
    biallelicvcfInput="/opt/tsinghua/cfDNApipeTest/file/small_exac_common_3_hg19.SNP_biallelic.vcf",
    verbose=False,
)
res_case12 = contamination(upstream=res_case11, verbose=False)
res_case13 = runCounter(upstream=res_case6, verbose=False)
res_case14 = qualimap(upstream=res_case6, verbose=False)
res_case16 = bam2bed(upstream=res_case7, verbose=False)

# snv / InDel
# res1 = mutect2t(
#     caseupstream=res_case12,
#     ctrlupstream=res_ctrl13,
#     vcfInput="/opt/tsinghua/cfDNApipeTest/file/af-only-gnomad.raw.sites.hg19.vcf.gz",
# )
# res2 = filterMutectCalls(upstream=res1, verbose=True)
# res3 = gatherVCF(upstream=res2)
# res4 = annovar(
#     upstream=res3,
#     plInput="/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/table_annovar.pl",
#     dbdir="/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/humandb/",
# )
# res5 = annovarStat(upstream=res4,)


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
