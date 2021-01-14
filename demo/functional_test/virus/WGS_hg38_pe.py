from cfDNApipe import *

pipeConfigure2(
    threads=8,
    genome="hg19",
    refdir="/disk1/wzhang/Genome/hg38",
    outdir="/disk1/wzhang/pipeline_test/hg38_pe",
    data="WGS",
    type="paired",
    JavaMem="8G",
    case="cancer",
    ctrl="normal",
    build=True,
)

fi_caseOut, fi_ctrlOut = cfDNAWGS2(
    caseFolder="/disk1/wzhang/pipeline_test/WGS_pe/case",
    ctrlFolder="/disk1/wzhang/pipeline_test/WGS_pe/ctrl",
    caseName="cancer",
    ctrlName="normal",
    idAdapter=True,
    rmAdapter=True,
    rmAdOP={"--qualitybase": 33, "--gzip": True},
    bowtie2OP={"-q": True, "-N": 1, "--time": True},
    dudup=True,
    CNV=False,
    armCNV=False,
    fragProfile=False,
    OCF=False,
    report=False,
    verbose=False,
)