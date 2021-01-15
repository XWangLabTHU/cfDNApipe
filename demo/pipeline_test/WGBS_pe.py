from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bismark",
    outdir=r"/opt/tsinghua/zhangwei/Pipeline_test/o_WGBS-PE",
    data="WGBS",
    type="paired",
    build=True,
    JavaMem="10g",
)

res = cfDNAWGBS(
    inputFolder=r"/opt/tsinghua/zhangwei/Pipeline_test/WGBS-PE",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    CNV=True,
    armCNV=True,
    fragProfile=False,
    report=True,
    verbose=False,
)
