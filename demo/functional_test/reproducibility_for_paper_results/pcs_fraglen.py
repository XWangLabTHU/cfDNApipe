from cfDNApipe import *
import glob
import pysam

pipeConfigure(
    threads=20,
    genome="hg19",
    refdir=r"/home/wzhang/genome/hg19",
    outdir=r"/data/wzhang/pcs_final/pcs_fraglen",
    data="WGS",
    type="paired",
    JavaMem="10G",
    build=True,
)

verbose = False

case_bed = glob.glob("/data/wzhang/pcs_final/HCC/*.bed")
ctrl_bed = glob.glob("/data/wzhang/pcs_final/Healthy/*.bed")

res_fraglenplot_comp = fraglenplot_comp(
    casebedInput=case_bed, ctrlbedInput=ctrl_bed, caseupstream=True, verbose=verbose)
