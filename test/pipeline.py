# WGS paired end pipeline

from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bowtie2",
    outdir=r"/home/zhangwei/pipeline-for-paired-WGS",
    data="WGS",
    type="paired",
    JavaMem="8G",
    build=True,
)

res = cfDNAWGS(
    inputFolder=r"/home/zhangwei/pipeline-for-paired-WGBS/raw",
    idAdapter=True,
    rmAdapter=False,
    dudup=True,
    CNV=True,
    fragProfile=True,
    verbose=True,
)
