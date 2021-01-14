from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bowtie2",
    outdir="/opt/tsinghua/zhangwei/Pipeline_SNV/WGS_hg19_se_without_control",
    data="WGS",
    type="single",
    JavaMem="8G",
    build=True,
)

res = cfDNAWGS(
    inputFolder=r"/opt/tsinghua/zhangwei/Pipeline_SNV/WGS_se/case",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    CNV=False,
    armCNV=False,
    fragProfile=False,
    verbose=False,
)

Configure.snvRefCheck(folder="/home/zhangwei/Genome/SNV_hg19", build=True)

# Using bam files directly.
# Of course, the "upstream" of addRG can be from "rmduplicate".
res1 = addRG(upstream=res.rmduplicate)

res2 = BaseRecalibrator(upstream=res1, knownSitesDir=Configure.getConfig("snv.folder"))
res3 = BQSR(upstream=res2)
res4 = getPileup(upstream=res3, biallelicvcfInput=Configure.getConfig("snv.ref")["7"],)
res5 = contamination(upstream=res4)

res6 = mutect2t(
    caseupstream=res5, vcfInput=Configure.getConfig("snv.ref")["6"], ponbedInput=Configure.getConfig("snv.ref")["8"],
)

res7 = filterMutectCalls(upstream=res6)

# ???
res8 = gatherVCF(upstream=res7)

# split somatic mutations
res9 = bcftoolsVCF(upstream=res8, stepNum="somatic")

# split germline mutations
res10 = bcftoolsVCF(upstream=res8, other_params={"-f": "'germline'"}, suffix="germline", stepNum="germline")
