from cfDNApipe import *

pipeConfigure(
    threads=20,
    genome="hg19",
    refdir=r"/home/wzhang/genome/hg19_bismark",
    outdir=r"/data/wzhang/pipeline_test/pipeline-for-single-WGBS",
    data="WGBS",
    type="single",
    build=True,
)

res1 = inputprocess(inputFolder=r"/data/wzhang/pipeline_test/pipeline-for-single-WGBS/raw")
res2 = fastqc(upstream=res1, verbose=False)
res3 = adapterremoval(upstream=res1, other_params={"--qualitybase": 64, "--gzip": True}, verbose=False)
res4 = bismark(upstream=res3, other_params={"-q": True, "--phred64-quals": True, "--bowtie2": True, "--un": True,}, verbose=False)
res5 = bismark_deduplicate(upstream=res4, verbose=False)
res6 = bismark_methylation_extractor(upstream=res5, verbose=False)
res7 = compress_methyl(upstream=res6, verbose=False)
res8 = calculate_methyl(upstream=res7, verbose=False)
res9 = bamsort(upstream=res5, verbose=False)
res10 = bam2bed(upstream=res9, verbose=False)

# CNV sub step
res11 = runCounter(upstream=res9, filetype=1, verbose=False, stepNum="CNV01")
res12 = runCounter(filetype=0, upstream=True, verbose=False, stepNum="CNV02")
res13 = GCCorrect(readupstream=res11, gcupstream=res12, verbose=False, stepNum="CNV03")
