from cfDNApipe import *
import argparse
import os

def main(args):
    if not os.path.exists(args.outputDir):
        os.mkdir(args.outputDir)
    
    pipeConfigure(
        threads=10,
        genome="hg19",
        refdir="/home/zhangwei/Genome/hg19_bowtie2",
        outdir=args.outputDir,
        data="WGS",
        type="paired",
        JavaMem="10G",
        build=True,
    )

    Configure.snvRefCheck(folder="/home/zhangwei/Genome/SNV_hg19", build=True)

    bams = [args.input]

    # Using bam files directly.
    # Of course, the "upstream" of addRG can be from "rmduplicate".
    res1 = addRG(bamInput=bams, upstream=True)

    res2 = BaseRecalibrator(upstream=res1, knownSitesDir=Configure.getConfig("snv.folder"))
    res3 = BQSR(upstream=res2)
    res4 = getPileup(upstream=res3, biallelicvcfInput=Configure.getConfig("snv.ref")["7"],)
    res5 = contamination(upstream=res4)

    res6 = mutect2t(
        caseupstream=res5, vcfInput=Configure.getConfig("snv.ref")["6"], ponbedInput=Configure.getConfig("snv.ref")["8"],
    )

    # res7 = filterMutectCalls(upstream=res6, other_params = {"--f-score-beta": 0.5})
    res7 = filterMutectCalls(upstream=res6)

    # ???
    res8 = gatherVCF(upstream=res7)

    # split somatic mutations
    res9 = bcftoolsVCF(upstream=res8, stepNum="somatic")

    # split germline mutations
    res10 = bcftoolsVCF(upstream=res8, other_params={"-f": "'germline'"}, suffix="germline", stepNum="germline")


def parse_arguments(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("outputDir", type=str)
    parser.add_argument("input", type=str)
    return parser.parse_args(argv)


if __name__ == "__main__":
    main(parse_arguments(sys.argv[1:]))