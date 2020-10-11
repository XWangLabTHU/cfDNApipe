from cfDNApipe import *
import glob

pipeConfigure2(
    threads=120,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bismark",
    outdir="/opt/tsinghua/zhangwei/302_20200902",
    data="WGBS",
    type="paired",
    JavaMem="10G",
    case="Re_HCC",
    ctrl="Re_HBV_Healthy",
    build=True,
)

verbose = False

switchConfigure("Re_HCC")

case_bams = glob.glob(
    "/opt/tsinghua/zhangwei/302_20200902/HCC/intermediate_result/step_06_bamsort/*.bam"
)
case_bed = glob.glob(
    "/opt/tsinghua/zhangwei/302_20200902/HCC/intermediate_result/step_07_bam2bed/*.bed.gz"
)

case_res_bamCounter = bamCounter(
    bamInput=case_bams, upstream=True, verbose=verbose, stepNum="ARMCNV01"
)
case_cnv_gcCounter = runCounter(
    filetype=0, upstream=True, verbose=verbose, stepNum="ARMCNV02"
)
case_cnv_GCCorrect = GCCorrect(
    readupstream=case_res_bamCounter,
    gcupstream=case_cnv_gcCounter,
    verbose=verbose,
    stepNum="ARMCNV03",
)

case_fp_fragCounter = fpCounter(
    bedgzInput=case_bed, upstream=True, verbose=verbose, stepNum="FP01", processtype=1
)
case_fp_gcCounter = runCounter(
    filetype=0, binlen=5000000, upstream=True, verbose=verbose, stepNum="FP02"
)
case_fp_GCCorrect = GCCorrect(
    readupstream=case_fp_fragCounter,
    gcupstream=case_fp_gcCounter,
    readtype=2,
    corrkey="-",
    verbose=verbose,
    stepNum="FP03",
)

switchConfigure("Re_HBV_Healthy")

ctrl_bams = glob.glob(
    "/opt/tsinghua/zhangwei/302_20200902/HBV_Healthy/intermediate_result/step_06_bamsort/*.bam"
)

ctrl_bed = glob.glob(
    "/opt/tsinghua/zhangwei/302_20200902/HBV_Healthy/intermediate_result/step_07_bam2bed/*.bed.gz"
)

ctrl_res_bamCounter = bamCounter(
    bamInput=ctrl_bams, upstream=True, verbose=verbose, stepNum="ARMCNV01",
)
ctrl_cnv_gcCounter = runCounter(
    filetype=0, upstream=True, verbose=verbose, stepNum="ARMCNV02"
)
ctrl_cnv_GCCorrect = GCCorrect(
    readupstream=ctrl_res_bamCounter,
    gcupstream=ctrl_cnv_gcCounter,
    verbose=verbose,
    stepNum="ARMCNV03",
)

ctrl_fp_fragCounter = fpCounter(
    bedgzInput=ctrl_bed, upstream=True, verbose=verbose, stepNum="FP01", processtype=1
)
ctrl_fp_gcCounter = runCounter(
    filetype=0, binlen=5000000, upstream=True, verbose=verbose, stepNum="FP02"
)
ctrl_fp_GCCorrect = GCCorrect(
    readupstream=ctrl_fp_fragCounter,
    gcupstream=ctrl_fp_gcCounter,
    readtype=2,
    corrkey="-",
    verbose=verbose,
    stepNum="FP03",
)

switchConfigure("Re_HCC")

res_computeCNV = computeCNV(
    caseupstream=case_cnv_GCCorrect,
    ctrlupstream=ctrl_cnv_GCCorrect,
    stepNum="ARMCNV",
    verbose=verbose,
)

res_fragprofplot = fragprofplot(
    caseupstream=case_fp_GCCorrect, ctrlupstream=ctrl_fp_GCCorrect, stepNum="FP04",
)
