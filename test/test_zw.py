from cfDNApipe import *

Configure.setData("WGBS")
Configure.setThreads(30)
Configure.setGenome("hg19")
Configure.setRefDir(r"/home/wzhang/genome/hg19_bismark")
Configure.setOutDir(r"/data/wzhang/pipeline-for-paired-WGBS")
Configure.pipeFolderInit()
Configure.refCheck(build=True)

res_case1 = inputprocess(inputFolder=r"/data/wzhang/pipeline-for-paired-WGBS/raw")
res_case2 = fastqc(upstream=res_case1)
res_case3 = identifyAdapter(upstream=res_case1)
res_case4 = adapterremoval(upstream=res_case3)
res_case5 = bismark(upstream=res_case4)
res_case6 = bismark_deduplicate(upstream=res_case5)
res_case7 = bismark_methylation_extractor(upstream=res_case6)
res_case8 = compress_methyl(upstream=res_case7)
res_case9 = calculate_methyl(upstream=res_case8)
res_case10 = bamsort(upstream=res_case6)
res_case11 = bam2bed(upstream=res_case10)
res_case12 = fraglenplot(upstream=res_case11)


# single end WGBS

from cfDNApipe import *

Configure.setData("WGBS")
Configure.setType("single")
Configure.setThreads(20)
Configure.setGenome("hg19")
Configure.setRefDir(r"/home/wzhang/genome/hg19_bismark")
Configure.setOutDir(r"/data/wzhang/pipeline-for-single-WGBS")
Configure.pipeFolderInit()
Configure.refCheck(build=True)


res1 = inputprocess(inputFolder=r"/data/wzhang/pipeline-for-single-WGBS/raw")
res2 = fastqc(upstream=res1)
res3 = adapterremoval(
    upstream=res1, formerrun=res2, other_params={"--qualitybase": 64, "--gzip": True}
)
res4 = bismark(
    upstream=res3,
    other_params={
        "-q": True,
        "--phred64-quals": True,
        "-N": 1,
        "--bowtie2": True,
        "--un": True,
    },
)
res5 = bamsort(upstream=res4)
res6 = rmduplicate(upstream=res5)
res7 = bam2bed(upstream=res6)
res8 = computemethyl(upstream=res6, formerrun=res7)
res9 = addRG(upstream=res6, formerrun=res8)


# paired end WGS

from cfDNApipe import *

Configure.setData("WGS")
Configure.setThreads(20)
Configure.setGenome("hg19")
Configure.setRefDir(r"/home/wzhang/genome/hg19")
Configure.setOutDir(r"/data/wzhang/pipeline-for-paired-WGS")
Configure.pipeFolderInit()
Configure.refCheck(build=True)


res1 = inputprocess(inputFolder=r"/data/wzhang/pipeline-for-paired-WGS/raw")
res2 = fastqc(upstream=res1)
res3 = identifyAdapter(upstream=res1, formerrun=res2)
res4 = adapterremoval(upstream=res3)
res5 = bowtie2(upstream=res4)
res6 = bamsort(upstream=res5)
res7 = rmduplicate(upstream=res6)
res8 = bam2bed(upstream=res7)
res9 = fraglenplot(upstream=res8)
res10 = addRG(upstream=res7, formerrun=res9)


# single end WGS

from cfDNApipe import *

Configure.setData("WGS")
Configure.setType("single")
Configure.setThreads(20)
Configure.setGenome("hg19")
Configure.setRefDir(r"/home/wzhang/genome/hg19")
Configure.setOutDir(r"/data/wzhang/pipeline-for-single-WGS")
Configure.pipeFolderInit()
Configure.refCheck(build=True)

res1 = inputprocess(inputFolder=r"/data/wzhang/pipeline-for-single-WGS/raw")
res2 = fastqc(upstream=res1)
res3 = adapterremoval(
    upstream=res1, formerrun=res2, other_params={"--qualitybase": 33, "--gzip": True}
)
res4 = bowtie2(upstream=res3)
res5 = bamsort(upstream=res4)
res6 = rmduplicate(upstream=res5)
res7 = bam2bed(upstream=res6)
res9 = addRG(upstream=res6, formerrun=res7)

