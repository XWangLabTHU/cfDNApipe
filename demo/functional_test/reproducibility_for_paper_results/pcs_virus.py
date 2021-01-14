from cfDNApipe import *
import glob

pipeConfigure(
    threads=20,
    genome="hg19",
    refdir=r"/home/wzhang/genome/hg19",
    outdir=r"/data/wzhang/pcs_final/pcs_virus",
    data="WGS",
    type="paired",
    build=True,
    JavaMem="15g",
)
# virus detect


Configure.virusGenomeCheck(folder="/data/wzhang/pcs_final/pcs_virus/database", build=True)


# paired
fq1 = glob.glob("/data/data/pcs-output/Others/*.fq.1.gz")
fq2 = glob.glob("/data/data/pcs-output/Others/*.fq.2.gz")

virusdetect(seqInput1=fq1, seqInput2=fq2, upstream=True)

