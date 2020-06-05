# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:58:34 2019

@author: Shuying He
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure
import pkg_resources

__metaclass__ = type


class annovar(StepBase):
    def __init__(
        self,
        plInput=None,
        vcfInput=None,
        avinput=None,
        outputdir=None,
        threads=1,
        dbdir=None,
        ref="hg19",
        annodb=["refGene"],
        other_params={
            "--remove": True,
            "--nastring": ".",
            "--vcfinput": True,
            "--intronhgvs": 300,
        },
        stepNum=None,
        upstream=None,
        verbose=True,
        **kwargs
    ):

        super(annovar, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
            if vcfInput:
                self.setInput("annoInput", vcfInput)
            elif avinput:
                self.setInput("annoInput", avinput)
                other_params["--vcfInput"] = False
            else:
                raise commonError(
                    "Please give files for annovar by vcfInput or avinput."
                )
        else:
            # check Configure for running pipeline
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ in ["gatherVCF", "bcftoolsVCF"]:
                # using identified adapters
                self.setInput("annoInput", upstream.getOutput("vcfOutput"))
            else:
                raise commonError(
                    "Parameter upstream must from gatherVCF or bcftoolsVCF."
                )
        self.checkInputFilePath()

        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("annoInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
            self.setParam("threads", threads)
            self.setParam("buildver", ref)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())
            self.setParam("buildver", Configure.getGenome())

        if plInput:
            self.setInput("plInput", plInput)
        else:
            self.setInput(
                "plInput",
                pkg_resources.resource_filename(
                    "cfDNApipe", "data/annovar/annovar/table_annovar.pl"
                ),
            )

        if dbdir:
            self.setInput("dbdir", dbdir)
        else:
            self.setInput(
                "dbdir",
                pkg_resources.resource_filename(
                    "cfDNApipe", "data/annovar/annovar/humandb/"
                ),
            )

        self.setParam("protocol", ",".join(annodb))
        self.setParam("operation", self.GetAnnodbType(annodb))

        self.setParam("other_params", other_params)

        self.setParam(
            "prefix",
            [
                os.path.join(
                    self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)
                )
                for x in self.getInput("annoInput")
            ],
        )

        self.setOutput(
            "txtOutput",
            [
                x + ".%s_multianno.txt" % self.getParam("buildver")
                for x in self.getParam("prefix")
            ],
        )

        multi_run_len = len(self.getInput("annoInput"))
        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(
                [
                    self.getInput("plInput"),
                    self.getInput("annoInput")[i],
                    self.getInput("dbdir"),
                    "--out",
                    self.getParam("prefix")[i],
                    "--buildver",
                    self.getParam("buildver"),
                    "--protocol",
                    self.getParam("protocol"),
                    "--operation",
                    self.getParam("operation"),
                    "--thread",
                    self.getParam("threads"),
                    self.getParam("other_params"),
                ]
            )
            all_cmd.append(tmp_cmd)

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            if verbose:
                self.run(all_cmd)
            else:
                self.multiRun(args=all_cmd, func=None, nCore=1)

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)

    def GetAnnodbType(self, annodb):
        annodb_type = {}
        # Gene annotation
        Gene_anno = ["refGene"]
        for i in Gene_anno:
            annodb_type.setdefault(i, "g")

        # Region annotation
        Region_anno = [
            "targetScanS",
            "rmsk",
            "cytoBand",
            "CpgIslandExt",
            "genomicSuperDups",
            "tfbsConsSites",
            "phastConsElements46way",
            "ucscGenePfam",
        ]
        for i in Region_anno:
            annodb_type.setdefault(i, "r")

        # Filter annotation
        Filter_anno = [
            "snp151",
            "popfreq_all_20150413",
            "gnomad_exome",
            "gnomad_genome",
            "regsnpintron",
            "CGI_biomarkers",
            "CGI_OncogenicMutations",
            "CIViC_ClinicalEvidence",
            "clinvar_20200316",
            "cosmic70",
            "DoCM_V3.2",
            "HGMD_disease_pubmedID",
            "hgmd_id",
            "icgc21",
            "dbnsfp33a",
            "nci60",
            "cg69",
            "hrcr1",
            "intervar_20180118",
            "wgEncodeAwgSegmentationChromhmmGm12878",
            "wgEncodeAwgSegmentationChromhmmH1hesc",
            "wgEncodeAwgSegmentationChromhmmHelas3",
            "wgEncodeAwgSegmentationChromhmmHepg2",
            "wgEncodeAwgSegmentationChromhmmHuvec",
        ]

        for i in Filter_anno:
            annodb_type.setdefault(i, "f")

        # convert the protocol to operation
        for i in annodb:
            if i not in annodb_type:
                raise commonError(i + " not in annovar database, please check!!")

        out = [annodb_type[x] for x in annodb]

        return ",".join(out)
