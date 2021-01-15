from cfDNApipe import *

pipeConfigure2(
    threads=60,
    genome="hg38",
    refdir="/home/zhangwei/Genome/hg38_bowtie2",
    outdir="/opt/tsinghua/zhangwei/Pipeline_test/o_WGS-compare",
    data="WGS",
    type="paired",
    JavaMem="8G",
    case="cancer",
    ctrl="normal",
    build=True,
)

# fragProfile is set to False because the demo data is too small to get valid values
a, b = cfDNAWGS2(
    caseFolder="/opt/tsinghua/zhangwei/Pipeline_test/WGS-compare/case",
    ctrlFolder="/opt/tsinghua/zhangwei/Pipeline_test/WGS-compare/ctrl",
    caseName="cancer",
    ctrlName="normal",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    CNV=True,
    armCNV=True,
    fragProfile=False,
    OCF=True,
    report=True,
    verbose=False,
)