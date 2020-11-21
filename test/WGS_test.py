from cfDNApipe import *

pipeConfigure(
    threads=120,
    genome="hg19",
    refdir=r"/home/zhangwei/Genome/hg19_bowtie2",
    outdir=r"/opt/tsinghua/zhangwei/to_zw",
    data="WGS",
    type="paired",
    build=True,
    JavaMem="15g",
)

res = cfDNAWGS(
    inputFolder=r"/opt/tsinghua/zhangwei/to_zw/raw",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    CNV=True,
    armCNV=True,
    fragProfile=True,
    verbose=False,
)

# virus detect
vir1 = unmapfasta(upstream=res.bowtie2)
vir2 = virusdetect(
    upstream=vir1,
    virusDB="/opt/tsinghua/cfDNApipeTest/file/virus_genome/viral_REFSEQ.fa",
    blastnIdxH="/opt/tsinghua/cfDNApipeTest/file/hg19/bowtie2/hg19",
    blastnIdxV="/opt/tsinghua/cfDNApipeTest/file/virus_genome/viral_REFSEQ",
    virusIDfile="/opt/tsinghua/cfDNApipeTest/file/virus_genome/virus_name_list.txt",
)

# SNV
res8 = addRG(upstream=res.rmduplicate)
res9 = BaseRecalibrator(
    upstream=res8, knownSitesDir=r"/opt/tsinghua/cfDNApipeTest/file/vcf"
)
res10 = BQSR(upstream=res9)
res11 = getPileup(
    upstream=res10,
    biallelicvcfInput="/opt/tsinghua/cfDNApipeTest/file/small_exac_common_3_hg19.SNP_biallelic.vcf",
)
res12 = contamination(upstream=res11)
res13 = mutect2t(
    caseupstream=res12,
    vcfInput="/opt/tsinghua/cfDNApipeTest/file/af-only-gnomad.raw.sites.hg19.vcf.gz",
    ponbedInput="/opt/tsinghua/cfDNApipeTest/file/vcf/pon/somatic-hg19_Mutect2-WGS-panel.vcf.gz",
)
res16 = filterMutectCalls(upstream=res13)

# split somatic mutations
res17 = bcftoolsVCF(upstream=res16)
res18 = annovar(
    upstream=res17,
    plInput=r"/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/table_annovar.pl",
    dbdir=r"/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/humandb/",
)
res19 = annovarStat(upstream=res18)

# split germline mutations
res20 = bcftoolsVCF(
    upstream=res16, other_params={"-f": "'germline'"}, suffix="germline", stepNum=20
)
res21 = annovar(
    upstream=res20,
    plInput=r"/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/table_annovar.pl",
    dbdir=r"/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/humandb/",
)
res22 = annovarStat(upstream=res21)
