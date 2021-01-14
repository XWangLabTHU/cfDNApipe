from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bowtie2",
    outdir=r"/opt/tsinghua/zhangwei/Pipeline_test/o_WGS-PE",
    data="WGS",
    type="paired",
    build=True,
    JavaMem="10g",
)

res = cfDNAWGS(
    inputFolder=r"/opt/tsinghua/zhangwei/Pipeline_test/WGS-PE",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    CNV=True,
    armCNV=True,
    fragProfile=True,
    verbose=False,
    report=True,
)