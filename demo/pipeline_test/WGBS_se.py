from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bismark",
    outdir=r"/opt/tsinghua/zhangwei/Pipeline_test/o_WGBS-SE",
    data="WGBS",
    type="single",
    build=True,
    JavaMem="10g",
)

res = cfDNAWGBS(
    inputFolder=r"/opt/tsinghua/zhangwei/Pipeline_test/WGBS-SE",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    CNV=True,
    armCNV=True,
    fragProfile=True,
    verbose=True,
    report=True,
)

