from cfDNApipe import *

pipeConfigure2(
    threads=120,
    genome="hg38",
    refdir="/home/zhangwei/Genome/hg38_bismark",
    outdir="/opt/tsinghua/zhangwei/Pipeline_test/o_WGBS-compare",
    data="WGBS",
    type="paired",
    JavaMem="10G",
    case="cancer",
    ctrl="normal",
    build=True,
)

# fragProfile is set to False because the demo data is too small to get valid values
a, b = cfDNAWGBS2(
    caseFolder="/opt/tsinghua/zhangwei/Pipeline_test/WGBS-compare/case",
    ctrlFolder="/opt/tsinghua/zhangwei/Pipeline_test/WGBS-compare/ctrl",
    caseName="cancer",
    ctrlName="normal",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    armCNV=True,
    CNV=True,
    fragProfile=False,
    deconvolution=True,
    OCF=True,
    report=True,
    verbose=False,
)

DMR = computeDMR(caseupstream=a.calculate_methyl, ctrlupstream=b.calculate_methyl)
