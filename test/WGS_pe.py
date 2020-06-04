# paired end WGS

from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bowtie2",
    outdir=r"/home/zhangwei/pipeline-for-paired-WGS",
    data="WGS",
    type="paired",
    JavaMem="8G",
    build=True,
)

res1 = inputprocess(inputFolder=r"/home/zhangwei/pipeline-for-paired-WGS/raw")
res2 = fastqc(upstream=res1, verbose=False)
res3 = identifyAdapter(upstream=res1)
res4 = adapterremoval(upstream=res3)
res5 = bowtie2(upstream=res4)
res6 = bamsort(upstream=res5)
res7 = rmduplicate(upstream=res6, verbose=False)

res8 = bam2bed(upstream=res7)
res9 = fraglenplot(upstream=res8)

# CNV sub step
res13 = runCounter(upstream=res7, filetype=1, verbose=False, stepNum="CNV01")
res14 = runCounter(filetype=0, upstream=True, verbose=False, stepNum="CNV02")
res15 = GCCorrect(readupstream=res13, gcupstream=res14, verbose=False, stepNum="CNV03")

# Fragmentation Profile sub step
res16 = runCounter(filetype=0, binlen=5000000, upstream=True, verbose=False, stepNum="FP01")
res17 = fpCounter(upstream=res8, verbose=False, stepNum="FP02")
res18 = GCCorrect(readupstream=res17, gcupstream=res16, readtype=2, corrkey="-", verbose=False, stepNum="FP03")

report_generator(
    fastqcRes=res2, 
    identifyAdapterRes=res3, 
    rmduplicateRes=res7, 
    fraglenplotRes=res9, 
    CNV_GCcorrectRes=res15,
    fragprof_GCcorrectRes=res18,
)





from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bowtie2",
    outdir=r"/home/zhangwei/pipeline-for-paired-WGS",
    data="WGS",
    type="paired",
    JavaMem="8G",
    build=True,
)

res = cfDNAWGS(
    inputFolder=r"/home/zhangwei/pipeline-for-paired-WGBS/raw",
    idAdapter=True,
    rmAdapter=False,
    dudup=True,
    CNV=True,
    fragProfile=True,
    verbose=True,
)