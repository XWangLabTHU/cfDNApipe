# arm level CNV

from cfDNApipe import *
import glob

pipeConfigure2(
    threads=120,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bowtie2",
    outdir="/home/zhangwei/pipeline_vliad_test/WGS_pe",
    data="WGS",
    type="paired",
    JavaMem="8G",
    case="cancer",
    ctrl="normal",
    build=True,
)

verbose = False


case_bam = glob.glob("/opt/tsinghua/zhangwei/CUHK_PE/HCC/*.bam")
ctrl_bam = glob.glob("/opt/tsinghua/zhangwei/CUHK_PE/Healthy/*.bam")

switchConfigure("cancer")

case1 = bamCounter(
    bamInput=case_bam, verbose=verbose, stepNum="ARMCNV01", upstream=True
)
case2 = runCounter(filetype=0, verbose=verbose, stepNum="ARMCNV02", upstream=True)
case3 = GCCorrect(
    readupstream=case1, gcupstream=case2, verbose=verbose, stepNum="ARMCNV03",
)

switchConfigure("normal")

ctrl1 = bamCounter(
    bamInput=ctrl_bam, verbose=verbose, stepNum="ARMCNV01", upstream=True
)
ctrl2 = runCounter(filetype=0, verbose=verbose, stepNum="ARMCNV02", upstream=True)
ctrl3 = GCCorrect(
    readupstream=ctrl1, gcupstream=ctrl2, verbose=verbose, stepNum="ARMCNV03",
)

switchConfigure("cancer")

armcnv = computeCNV(
    caseupstream=case3, ctrlupstream=ctrl3, stepNum="ARMCNV", verbose=verbose,
)


from cfDNApipe import *
import glob

pipeConfigure2(
    threads=120,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bowtie2",
    outdir="/home/zhangwei/pipeline_vliad_test/WGS_pe",
    data="WGS",
    type="paired",
    JavaMem="8G",
    case="cancer",
    ctrl="normal",
    build=True,
)

verbose = False

# frag profile
case_bed = glob.glob("/opt/tsinghua/zhangwei/CUHK_PE/HCC_bed/*.bed.gz")
ctrl_bed = glob.glob("/opt/tsinghua/zhangwei/CUHK_PE/Healthy_bed/*.bed.gz")

switchConfigure("cancer")

case1 = fpCounter(
    bedgzInput=case_bed, verbose=verbose, stepNum="FP01", processtype=1, upstream=True
)
case2 = runCounter(
    filetype=0, binlen=5000000, upstream=True, verbose=verbose, stepNum="FP02"
)
case3 = GCCorrect(
    readupstream=case1,
    gcupstream=case2,
    readtype=2,
    corrkey="-",
    verbose=verbose,
    stepNum="FP03",
)

switchConfigure("normal")

ctrl1 = fpCounter(
    bedgzInput=ctrl_bed, verbose=verbose, stepNum="FP01", processtype=1, upstream=True
)
ctrl2 = runCounter(
    filetype=0, binlen=5000000, upstream=True, verbose=verbose, stepNum="FP02"
)
ctrl3 = GCCorrect(
    readupstream=ctrl1,
    gcupstream=ctrl2,
    readtype=2,
    corrkey="-",
    verbose=verbose,
    stepNum="FP03",
)

switchConfigure("cancer")

res_fragprofplot = fragprofplot(caseupstream=case3, ctrlupstream=ctrl3, stepNum="FP04",)


# OCF
from cfDNApipe import *
import glob

pipeConfigure2(
    threads=120,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bowtie2",
    outdir="/home/zhangwei/pipeline_vliad_test/WGS_pe",
    data="WGS",
    type="paired",
    JavaMem="8G",
    case="cancer",
    ctrl="normal",
    build=True,
)

verbose = False

# frag profile
case_bed = glob.glob("/opt/tsinghua/zhangwei/CUHK_PE/HCC_bed/*.bed")
ctrl_bed = glob.glob("/opt/tsinghua/zhangwei/CUHK_PE/Healthy_bed/*.bed")

switchConfigure("cancer")

res_computeOCF = computeOCF(
    casebedInput=case_bed,
    ctrlbedInput=ctrl_bed,
    verbose=verbose,
    caseupstream=True,
    stepNum="ocf01"
)

res_OCFplot = OCFplot(upstream=res_computeOCF, verbose=verbose, stepNum="ocf02")


















from cfDNApipe import *
import glob

pipeConfigure2(
    threads=120,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bismark",
    outdir="/home/zhangwei/pipeline_vliad_test/WGBS_se",
    data="WGBS",
    type="single",
    JavaMem="10G",
    case="cancer",
    ctrl="normal",
    build=True,
)

verbose = False

case_bam = glob.glob("/opt/tsinghua/zhangwei/CUHK_cfDNA_rmdup/Others/HOT*.bam")
ctrl_bam = glob.glob("/opt/tsinghua/zhangwei/CUHK_cfDNA_rmdup/Others/CTR*.bam")

switchConfigure("cancer")

case1 = bismark_methylation_extractor(
    bamInput=case_bam,
    other_params={
        "--no_overlap": True,
        "--report": True,
        "--no_header": True,
        "--gzip": True,
        "--bedGraph": True,
        "--zero_based": True,
    },
    upstream=True,
    verbose=verbose,
)
case2 = compress_methyl(upstream=case1, verbose=verbose)
case3 = calculate_methyl(
    upstream=case2, bedInput="plasmaMarkers_hg19.bed", verbose=verbose
)

case4 = calculate_methyl(
    upstream=case2, verbose=verbose, stepNum="cgi"
)

switchConfigure("normal")

ctrl1 = bismark_methylation_extractor(
    bamInput=ctrl_bam,
    other_params={
        "--no_overlap": True,
        "--report": True,
        "--no_header": True,
        "--gzip": True,
        "--bedGraph": True,
        "--zero_based": True,
    },
    upstream=True,
    verbose=verbose,
)
ctrl2 = compress_methyl(upstream=ctrl1, verbose=verbose)
ctrl3 = calculate_methyl(
    upstream=ctrl2, bedInput="plasmaMarkers_hg19.bed", verbose=verbose
)

ctrl4 = calculate_methyl(
    upstream=ctrl2, verbose=verbose, stepNum="cgi"
)

switchConfigure("cancer")

res_PCA = PCAplot(
    caseupstream=case3, ctrlupstream=ctrl3
)

res_DMR = computeDMR(
    caseupstream=case4, ctrlupstream=ctrl4
)


