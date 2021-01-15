from cfDNApipe import *

pipeConfigure2(
    threads=60,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bowtie2",
    outdir="/opt/tsinghua/zhangwei/Pipeline_SNV/WGS_hg19_pe_with_control",
    data="WGS",
    type="paired",
    JavaMem="8G",
    case="cancer",
    ctrl="normal",
    build=True,
)

fi_caseOut, fi_ctrlOut = cfDNAWGS2(
    caseFolder="/opt/tsinghua/zhangwei/Pipeline_SNV/WGS_pe/case",
    ctrlFolder="/opt/tsinghua/zhangwei/Pipeline_SNV/WGS_pe/ctrl",
    caseName="cancer",
    ctrlName="normal",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    CNV=False,
    armCNV=False,
    fragProfile=False,
    OCF=False,
    report=False,
    verbose=False,
)

# with control group, control group is used for to create PON file
Configure2.snvRefCheck(folder="/home/zhangwei/Genome/SNV_hg19", build=True)

# create PON file from normal samples
switchConfigure("normal")

ctrl_addRG = addRG(upstream=fi_ctrlOut.rmduplicate, stepNum="wc_PON01",)
ctrl_BaseRecalibrator = BaseRecalibrator(
    upstream=ctrl_addRG, knownSitesDir=Configure2.getConfig("snv.folder"), stepNum="wc_PON02",
)
ctrl_BQSR = BQSR(upstream=ctrl_BaseRecalibrator, stepNum="wc_PON03",)
ctrl_mutect2n = mutect2n(upstream=ctrl_BQSR, stepNum="wc_PON04",)
ctrl_dbimport = dbimport(upstream=ctrl_mutect2n, stepNum="wc_PON05",)
ctrl_createPon = createPON(upstream=ctrl_dbimport, stepNum="wc_PON06",)

# performing comparison between cancer and normal
switchConfigure("cancer")
case_addRG = addRG(upstream=fi_caseOut.rmduplicate, stepNum="SNV01",)
case_BaseRecalibrator = BaseRecalibrator(
    upstream=case_addRG, knownSitesDir=Configure2.getConfig("snv.folder"), stepNum="wc_SNV02",
)

case_BQSR = BQSR(upstream=case_BaseRecalibrator, stepNum="wc_SNV03")

case_getPileup = getPileup(
    upstream=case_BQSR, biallelicvcfInput=Configure2.getConfig("snv.ref")["7"], stepNum="wc_SNV04",
)
case_contamination = contamination(upstream=case_getPileup, stepNum="wc_SNV05")

# In this step, ponbedInput is ignored by using caseupstream parameter
case_mutect2t = mutect2t(
    caseupstream=case_contamination,
    ctrlupstream=ctrl_createPon,
    vcfInput=Configure2.getConfig("snv.ref")["6"],
    stepNum="wc_SNV06",
)

case_filterMutectCalls = filterMutectCalls(upstream=case_mutect2t, stepNum="wc_SNV07")

case_gatherVCF = gatherVCF(upstream=case_filterMutectCalls, stepNum="wc_SNV08")

# split somatic mutations
case_somatic = bcftoolsVCF(upstream=case_gatherVCF, stepNum="wc_somatic")

# split germline mutations
case_germline = bcftoolsVCF(
    upstream=case_gatherVCF, other_params={"-f": "'germline'"}, suffix="germline", stepNum="wc_germline"
)
