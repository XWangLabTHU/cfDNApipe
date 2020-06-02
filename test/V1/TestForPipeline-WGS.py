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

# case
switchConfigure("normal")
res6 = bamsort(bamInput=glob.glob("/home/zhangwei/pipeline-WGS-comp/bams/case_raw/*.bam"), upstream=True)
res7 = qualimap(upstream=res6)
res8 = rmduplicate(upstream=res6)
res9 = addRG(upstream=res8)

# WGS-SNV pipeline
res10 = BaseRecalibrator(
    upstream=res9, knownSitesDir=r"/opt/tsinghua/cfDNApipeTest/file/vcf"
)
res11 = BQSR(upstream=res10)

res12 = getPileup(
    upstream=res11,
    biallelicvcfInput="/opt/tsinghua/cfDNApipeTest/file/small_exac_common_3_hg19.SNP_biallelic.vcf",
)
res13 = contamination(upstream=res12)

res14 = mutect2t(
    caseupstream=res13,
    vcfInput="/opt/tsinghua/cfDNApipeTest/file/af-only-gnomad.raw.sites.hg19.vcf.gz",
    ponbedInput="/opt/tsinghua/cfDNApipeTest/file/vcf/pon/somatic-hg19_Mutect2-WGS-panel.vcf.gz",
)
res15 = filterMutectCalls(upstream=res14)
res16 = gatherVCF(upstream=res15)

res17 = annovar(
    upstream=res16,
    plInput="/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/table_annovar.pl",
    dbdir="/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/humandb/",
)
res18 = annovarStat(upstream=res17)


# virus detect
res18 = unmapfasta(
    upstream=res5,
    stepNum=18,
    plInput="/opt/tsinghua/cfDNApipeTest/software/VirusFinder2.0Plus/preprocessPlus_V2.pl",
)
res19 = virusdetect(
    upstream=res18,
    plInput="/opt/tsinghua/cfDNApipeTest/software/VirusFinder2.0Plus/detect_virusPlus.pl",
    virusDB="/opt/tsinghua/cfDNApipeTest/file/virus_genome/viral_REFSEQ.fa",
    blastnIdxH="/opt/tsinghua/cfDNApipeTest/file/hg19/bowtie2/hg19",
    blastnIdxV="/opt/tsinghua/cfDNApipeTest/file/virus_genome/viral_REFSEQ",
    pyscript="/opt/tsinghua/cfDNApipeTest/file/virus_genome/virusID2name.py",
    virusIDfile="/opt/tsinghua/cfDNApipeTest/file/virus_genome/virus_name_list.txt",
)
res20 = seqtk(upstream=res19)
res21 = BSVF(
    upstream1=res4,
    upstream2=res20,
    plInput="/opt/tsinghua/cfDNApipeTest/software/BSVF/BSVF_prepare_configFile.pl",
    bsuit="/opt/tsinghua/cfDNApipeTest/software/BSVF/bsuit",
    hostRef="/opt/tsinghua/cfDNApipeTest/file/hg19/hg19_EBV.fa",
)

# cnv
res22 = cnvbatch(
    caseupstream=res7,
    access="/opt/tsinghua/cfDNApipeTest/file/CNVkit/access-5kb-mappable.hg19.bed",
    annotate="/opt/tsinghua/cfDNApipeTest/file/CNVkit/refFlat_hg19.txt",
    stepNum=22,
)
res23 = cnvPlot(upstream=res22)
res24 = cnvTable(upstream=res22)
res25 = cnvHeatmap(upstream=res22)


print("analysis end...")
