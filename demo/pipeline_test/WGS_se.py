from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bowtie2",
    outdir=r"/opt/tsinghua/zhangwei/Pipeline_test/o_WGS-SE",
    data="WGS",
    type="single",
    JavaMem="8G",
    build=True,
)


res = cfDNAWGS(
    inputFolder=r"/opt/tsinghua/zhangwei/Pipeline_test/WGS-SE",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    CNV=True,
    armCNV=True,
    fragProfile=True,
    verbose=True,
    report=True,
)