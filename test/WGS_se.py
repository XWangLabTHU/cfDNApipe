from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bowtie2",
    outdir=r"/home/zhangwei/pipeline-for-single-WGS",
    data="WGS",
    type="single",
    build=True,
)

res1 = inputprocess(inputFolder=r"/home/zhangwei/pipeline-for-single-WGS/raw")
res2 = fastqc(upstream=res1, verbose=False)
res4 = adapterremoval(upstream=res1)
res5 = bowtie2(upstream=res4)
res6 = bamsort(upstream=res5)
res7 = rmduplicate(upstream=res6)
res8 = bam2bed(upstream=res7)

# CNV sub step
res13 = runCounter(upstream=res7, filetype=1, verbose=False, stepNum="CNV01")
res14 = runCounter(filetype=0, upstream=True, verbose=False, stepNum="CNV02")
res15 = GCCorrect(readupstream=res13, gcupstream=res14, verbose=False, stepNum="CNV03")
