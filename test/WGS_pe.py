# paired end WGS

from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bowtie2",
    outdir=r"/home/zhangwei/pipeline-for-paired-WGS",
    data="WGS",
    type="paired",
    JavaMem="8G",
    build=True,
)

# base processing
verbose = True

res_inputprocess = inputprocess(
    inputFolder=r"/home/zhangwei/pipeline-for-paired-WGS/raw"
)
res_fastqc = fastqc(upstream=res_inputprocess, verbose=verbose)
res_identifyAdapter = identifyAdapter(upstream=res_inputprocess, verbose=verbose)
res_adapterremoval = adapterremoval(upstream=res_identifyAdapter, verbose=verbose)
res_bowtie2 = bowtie2(upstream=res_adapterremoval, verbose=verbose)
res_bamsort = bamsort(upstream=res_bowtie2, verbose=verbose)
res_rmduplicate = rmduplicate(upstream=res_bamsort, verbose=verbose)
res_qualimap = qualimap(upstream=res_rmduplicate, verbose=verbose)
res_addRG = addRG(upstream=res_rmduplicate)

res_bam2bed = bam2bed(upstream=res_rmduplicate)
res_fraglenplot = fraglenplot(upstream=res_bam2bed)

# Arm-level CNV sub step
res_ARMCNV01 = runCounter(
    upstream=res_rmduplicate, filetype=1, verbose=verbose, stepNum="ARMCNV01"
)
res_ARMCNV02 = runCounter(
    filetype=0, upstream=True, verbose=verbose, stepNum="ARMCNV02"
)
res_ARMCNV03 = GCCorrect(
    readupstream=res_ARMCNV01,
    gcupstream=res_ARMCNV02,
    verbose=verbose,
    stepNum="ARMCNV03",
)

# # Fragmentation Profile sub step
# res_FP01 = runCounter(
#     filetype=0, binlen=5000000, upstream=True, verbose=verbose, stepNum="FP01"
# )
# res_FP02 = fpCounter(upstream=res_bam2bed, verbose=verbose, stepNum="FP02")
# res_FP03 = GCCorrect(
#     readupstream=res_FP02,
#     gcupstream=res_FP01,
#     readtype=2,
#     corrkey="-",
#     verbose=False,
#     stepNum="FP03",
# )

# SNV
res_BaseRecalibrator = BaseRecalibrator(
    upstream=res_addRG,
    knownSitesDir=r"/opt/tsinghua/cfDNApipeTest/file/vcf",
    verbose=verbose,
    stepNum="SNV01",
)
res_BQSR = BQSR(upstream=res_BaseRecalibrator, verbose=verbose, stepNum="SNV02")

res_getPileup = getPileup(
    upstream=res_BQSR,
    biallelicvcfInput="/opt/tsinghua/cfDNApipeTest/file/small_exac_common_3_hg19.SNP_biallelic.vcf",
    verbose=verbose,
    stepNum="SNV03",
)
res_contamination = contamination(
    upstream=res_getPileup, verbose=verbose, stepNum="SNV04"
)

res_mutect2t = mutect2t(
    caseupstream=res_contamination,
    vcfInput="/opt/tsinghua/cfDNApipeTest/file/af-only-gnomad.raw.sites.hg19.vcf.gz",
    ponbedInput="/opt/tsinghua/cfDNApipeTest/file/vcf/pon/somatic-hg19_Mutect2-WGS-panel.vcf.gz",
    verbose=verbose,
    stepNum="SNV05",
)
res_filterMutectCalls = filterMutectCalls(
    upstream=res_mutect2t, verbose=verbose, stepNum="SNV06"
)
res_gatherVCF = gatherVCF(
    upstream=res_filterMutectCalls, verbose=verbose, stepNum="SNV07"
)

res_annovar = annovar(
    upstream=res_gatherVCF,
    plInput="/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/table_annovar.pl",
    dbdir="/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/humandb/",
    verbose=verbose,
    stepNum="SNV08",
)
res_annovarStat = annovarStat(upstream=res_annovar)

# virus detect
res_unmapfasta = unmapfasta(
    upstream=res_bowtie2,
    plInput="/opt/tsinghua/cfDNApipeTest/software/VirusFinder2.0Plus/preprocessPlus_V2.pl",
    verbose=verbose,
    stepNum="VD01",
)
res_virusdetect = virusdetect(
    upstream=res_unmapfasta,
    plInput="/opt/tsinghua/cfDNApipeTest/software/VirusFinder2.0Plus/detect_virusPlus.pl",
    virusDB="/opt/tsinghua/cfDNApipeTest/file/virus_genome/viral_REFSEQ.fa",
    blastnIdxH="/opt/tsinghua/cfDNApipeTest/file/hg19/bowtie2/hg19",
    blastnIdxV="/opt/tsinghua/cfDNApipeTest/file/virus_genome/viral_REFSEQ",
    pyscript="/opt/tsinghua/cfDNApipeTest/file/virus_genome/virusID2name.py",
    virusIDfile="/opt/tsinghua/cfDNApipeTest/file/virus_genome/virus_name_list.txt",
    verbose=verbose,
    stepNum="VD02",
)
res_seqtk = seqtk(upstream=res_virusdetect, verbose=verbose, stepNum="VD02")
res_BSVF = BSVF(
    upstream1=res_adapterremoval,
    upstream2=res_seqtk,
    plInput="/opt/tsinghua/cfDNApipeTest/software/BSVF/BSVF_prepare_configFile.pl",
    bsuit="/opt/tsinghua/cfDNApipeTest/software/BSVF/bsuit",
    hostRef="/opt/tsinghua/cfDNApipeTest/file/hg19/hg19_EBV.fa",
    verbose=verbose,
    stepNum="VD02",
)

# cnv
res_cnvbatch = cnvbatch(
    caseupstream=res_rmduplicate,
    access="/opt/tsinghua/cfDNApipeTest/file/CNVkit/access-5kb-mappable.hg19.bed",
    annotate="/opt/tsinghua/cfDNApipeTest/file/CNVkit/refFlat_hg19.txt",
    verbose=verbose,
    stepNum="CNV01",
)
res_cnvPlot = cnvPlot(upstream=res_cnvbatch, verbose=verbose, stepNum="CNV02",)
res_cnvTable = cnvTable(upstream=res_cnvbatch, verbose=verbose, stepNum="CNV03",)
res_cnvHeatmap = cnvHeatmap(upstream=res_cnvbatch, verbose=verbose, stepNum="CNV04",)
