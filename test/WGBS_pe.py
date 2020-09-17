from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bismark",
    outdir=r"/home/zhangwei/pipeline-for-paired-WGBS",
    data="WGBS",
    type="paired",
    JavaMem="8G",
    build=True,
)

res1 = inputprocess(
    fqInput1=[
        "/home/zhangwei/pipeline-for-paired-WGBS/raw/TBR912.read1.fq.gz",
        "/home/zhangwei/pipeline-for-paired-WGBS/raw/TBR901.read1.fq",
    ],
    fqInput2=[
        "/home/zhangwei/pipeline-for-paired-WGBS/raw/TBR912.read2.fq.gz",
        "/home/zhangwei/pipeline-for-paired-WGBS/raw/TBR901.read2.fq",
    ],
)
res2 = fastqc(upstream=res1, verbose=True)
res3 = identifyAdapter(upstream=res1, verbose=True)
res4 = adapterremoval(upstream=res3, verbose=True)
res5 = bismark(upstream=res4, verbose=True)
res6 = bismark_deduplicate(upstream=res5, verbose=True)

res11 = addRG(upstream=res6, verbose=True)
res12 = BaseRecalibrator(
    upstream=res11, knownSitesDir=r"/opt/tsinghua/cfDNApipeTest/file/vcf", verbose=True
)
res13 = BQSR(upstream=res12)
res14 = qualimap(upstream=res13)


res7 = bismark_methylation_extractor(upstream=res6, verbose=False)
res8 = compress_methyl(upstream=res7, verbose=False)
res9 = calculate_methyl(upstream=res8, verbose=False)
res10 = bamsort(upstream=res6, verbose=False)
res11 = bam2bed(upstream=res10, verbose=False)
res12 = fraglenplot(upstream=res11, verbose=False)
res13 = computeDMR(upstream=res9)

# CNV sub step
res13 = runCounter(upstream=res10, filetype=1, verbose=False, stepNum="CNV01")
res14 = runCounter(filetype=0, upstream=True, verbose=False, stepNum="CNV02")
res15 = GCCorrect(readupstream=res13, gcupstream=res14, verbose=False, stepNum="CNV03")

# Fragmentation Profile sub step
res16 = runCounter(
    filetype=0, binlen=5000000, upstream=True, verbose=False, stepNum="FP01"
)
res17 = fpCounter(upstream=res11, verbose=False, stepNum="FP02")
res18 = GCCorrect(
    readupstream=res17,
    gcupstream=res16,
    readtype=2,
    corrkey="-",
    verbose=False,
    stepNum="FP03",
)

report_generator(
    fastqcRes=res2,
    identifyAdapterRes=res3,
    bismarkRes=res5,
    deduplicateRes=res6,
    fraglenplotRes=res12,
    CNV_GCcorrectRes=res15,
    fragprof_GCcorrectRes=res18,
)


from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bismark",
    outdir=r"/home/zhangwei/pipeline-for-paired-WGBS",
    data="WGBS",
    type="paired",
    build=True,
    JavaMem="10g",
)

res = cfDNAWGBS(
    inputFolder=r"/home/zhangwei/pipeline-for-paired-WGBS/raw",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    CNV=True,
    armCNV=True,
    fragProfile=True,
    report=True,
    verbose=False,
)


from cfDNApipe import *

pipeConfigure2(
    threads=20,
    genome="hg19",
    refdir="/home/wzhang/genome/hg19_bismark",
    outdir="/data/wzhang/pipeline_test/pipeline-WGBS-comp",
    data="WGBS",
    type="paired",
    JavaMem="8G",
    case="cancer",
    ctrl="normal",
    build=True,
)

case, ctrl, comp = cfDNAWGBS2(
    caseFolder="/data/wzhang/pipeline_test/pipeline-WGBS-comp/raw/case_large",
    ctrlFolder="/data/wzhang/pipeline_test/pipeline-WGBS-comp/raw/ctrl_large",
    caseName="cancer",
    ctrlName="tumor",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    armCNV=True,
    CNV=True,
    fragProfile=True,
    deconvolution=True,
    OCF=True,
    verbose=True,
)


from cfDNApipe import *

pipeConfigure2(
    threads=60,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bismark",
    outdir="/home/zhangwei/pipeline-WGBS-comp",
    data="WGBS",
    type="paired",
    JavaMem="8G",
    case="cancer",
    ctrl="normal",
    build=True,
)

a, b = cfDNAWGBS2(
    caseFolder="/home/zhangwei/pipeline-WGBS-comp/raw/case_small",
    ctrlFolder="/home/zhangwei/pipeline-WGBS-comp/raw/ctrl_small",
    caseName="cancer",
    ctrlName="normal",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    armCNV=False,
    CNV=False,
    fragProfile=False,
    deconvolution=False,
    OCF=False,
    report=True,
    verbose=True,
)
