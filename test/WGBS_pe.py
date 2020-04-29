from cfDNApipe import *

pipeConfigure(
    threads=20,
    genome="hg19",
    refdir=r"/home/wzhang/genome/hg19_bismark",
    outdir=r"/data/wzhang/pipeline_test/pipeline-for-paired-WGBS",
    data="WGBS",
    type="paired",
    build=True,
)

res1 = inputprocess(inputFolder=r"/data/wzhang/pipeline_test/pipeline-for-paired-WGBS/raw")
res2 = fastqc(upstream=res1, verbose=False)
res3 = identifyAdapter(upstream=res1, verbose=False)
res4 = adapterremoval(upstream=res3, verbose=False)
res5 = bismark(upstream=res4, verbose=False)
res6 = bismark_deduplicate(upstream=res5, verbose=False)
res7 = bismark_methylation_extractor(upstream=res6, verbose=False)
res8 = compress_methyl(upstream=res7, verbose=False)
res9 = calculate_methyl(upstream=res8, verbose=False)
res10 = bamsort(upstream=res6, verbose=False)
res11 = bam2bed(upstream=res10, verbose=False)
res12 = fraglenplot(upstream=res11, verbose=False)

# CNV sub step
res13 = runCounter(upstream=res10, filetype=1, verbose=False, stepNum="CNV01")
res14 = runCounter(filetype=0, upstream=True, verbose=False, stepNum="CNV02")
res15 = GCCorrect(readupstream=res13, gcupstream=res14, verbose=False, stepNum="CNV03")

# Fragmentation Profile sub step
res16 = runCounter(filetype=0, binlen=5000000, upstream=True, verbose=False, stepNum="FP01")
res17 = fpCounter(upstream=res11, verbose=False, stepNum="FP02")
res18 = GCCorrect(readupstream=res17, gcupstream=res16, readtype=2, corrkey="-", verbose=False, stepNum="FP03")




from cfDNApipe import *

pipeConfigure(
    threads=20,
    genome="hg19",
    refdir=r"/home/wzhang/genome/hg19_bismark",
    outdir=r"/data/wzhang/pipeline_test/pipeline-for-paired-WGBS",
    data="WGBS",
    type="paired",
    build=True,
)

cfDNAWGBS(inputFolder=r"/data/wzhang/pipeline_test/pipeline-for-paired-WGBS/raw",
          idAdapter=True, rmAdapter=True, dudup=True, CNV=False,
          fragProfile=True, verbose=True)



