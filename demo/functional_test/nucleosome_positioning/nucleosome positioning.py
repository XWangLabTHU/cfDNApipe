from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bowtie2",
    outdir="/opt/tsinghua/zhangwei/nucleosome_positioning",
    data="WGS",
    type="paired",
    JavaMem="8G",
    build=True,
)

res1 = bam2bed(bamInput=["SRR2130016-rmdup.bam"], upstream=True)

res2 = runWPS(upstream=res1, tsvInput="transcriptAnno-v19.tsv")

