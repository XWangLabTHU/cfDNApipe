from cfDNApipe import *

pipeConfigure(
    threads=120,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bismark",
    outdir=r"/opt/tsinghua/zhangwei/302_200723",
    data="WGBS",
    type="paired",
    build=True,
    JavaMem="10g",
)

res = cfDNAWGBS(
    inputFolder=r"/opt/tsinghua/novoseq/novo0720-5G/rawdata/HBV_Healthy",
    idAdapter=True,
    rmAdapter=True,
    bismarkOP={"-q": True, "--phred33-quals": True, "--bowtie2": True, "--un": True, "--score_min": "L,0,-0.6"},
    dudup=True,
    CNV=True,
    armCNV=True,
    fragProfile=True,
    report=True,
    verbose=False,
)


