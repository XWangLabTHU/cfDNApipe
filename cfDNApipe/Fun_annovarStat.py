# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:58:34 2019

@author: Shuying He
"""

import os
import pandas as pd
from collections import Counter
from .StepBase import StepBase
from .cfDNA_utils import commonError
from .Configure import Configure


def gtconvert(s):
    gt = s.split(":")[0]
    if "0" in gt:
        return "Het"
    else:
        return "Hom"


def treat_func_gene(x):
    if ";" in x:
        return x.split(";")[0]
    else:
        return x


class stat_annofile:
    def __init__(self, annofile):
        self.SampName = os.path.basename(annofile).replace(".hg19_multianno.txt", "")
        self.df = pd.read_table(
            annofile,
            sep="\t",
            usecols=["Ref", "Alt", "Func.refGene", "ExonicFunc.refGene", "Otherinfo13"],
        )

        self.snp_stat_dict = {
            "Basic_Statistics": {
                "SNP_Count": 0,
                "Transition": 0,
                "Transversion": 0,
                "Ti / Tv": 0,
                "Hom": 0,
                "Het": 0,
            },
            "Function_Statistics": {
                "intergenic": 0,
                "ncRNA_intronic": 0,
                "ncRNA_exonic": 0,
                "ncRNA_splicing": 0,
                "splicing": 0,
                "intronic": 0,
                "exonic": 0,
                "upstream": 0,
                "downstream": 0,
                "UTR3": 0,
                "UTR5": 0,
            },
            "Exonic_Statistics": {
                "nonsynonymous SNV": 0,
                "synonymous SNV": 0,
                "stopgain": 0,
                "stoploss": 0,
                "startloss": 0,
                "unknown": 0,
            },
        }
        self.indel_stat_dict = {
            "Basic_Statistics": {
                "InDel_Count": 0,
                "Ins": 0,
                "Del": 0,
                "Complex": 0,
                "Hom": 0,
                "Het": 0,
            },
            "Function_Statistics": {
                "intergenic": 0,
                "ncRNA_intronic": 0,
                "ncRNA_exonic": 0,
                "ncRNA_splicing": 0,
                "splicing": 0,
                "intronic": 0,
                "exonic": 0,
                "upstream": 0,
                "downstream": 0,
                "UTR3": 0,
                "UTR5": 0,
            },
            "Exonic_Statistics": {
                "frameshift deletion": 0,
                "frameshift insertion": 0,
                "nonframeshift deletion": 0,
                "nonframeshift insertion": 0,
                "nonframeshift substitution": 0,
                "unknown": 0,
                "stopgain": 0,
                "stoploss": 0,
                "startloss": 0,
            },
        }

        self.statistics = {"SNP": self.snp_stat_dict, "InDel": self.indel_stat_dict}

        self.annodf_convertion()
        self.annodf_statistics()
        print("sample %s stat finish!" % self.SampName)

    def annodf_convertion(self):
        # convert 0/1 -> Het, 1/1 Hom
        self.df["GT"] = self.df["Otherinfo13"].map(lambda x: gtconvert(x))
        # select the first function
        self.df["func"] = self.df["Func.refGene"].map(lambda x: treat_func_gene(x))

        # snp subtype: Ti Tv conversion
        self.df.loc[
            (self.df["Ref"] == "T") & (self.df["Alt"] == "C"), "subtype"
        ] = "Transition"
        self.df.loc[
            (self.df["Ref"] == "C") & (self.df["Alt"] == "T"), "subtype"
        ] = "Transition"
        self.df.loc[
            (self.df["Ref"] == "A") & (self.df["Alt"] == "G"), "subtype"
        ] = "Transition"
        self.df.loc[
            (self.df["Ref"] == "G") & (self.df["Alt"] == "A"), "subtype"
        ] = "Transition"

        self.df.loc[
            (self.df["Ref"] == "T") & (self.df["Alt"] == "A"), "subtype"
        ] = "Transversion"
        self.df.loc[
            (self.df["Ref"] == "T") & (self.df["Alt"] == "G"), "subtype"
        ] = "Transversion"
        self.df.loc[
            (self.df["Ref"] == "C") & (self.df["Alt"] == "A"), "subtype"
        ] = "Transversion"
        self.df.loc[
            (self.df["Ref"] == "C") & (self.df["Alt"] == "G"), "subtype"
        ] = "Transversion"
        self.df.loc[
            (self.df["Ref"] == "A") & (self.df["Alt"] == "T"), "subtype"
        ] = "Transversion"
        self.df.loc[
            (self.df["Ref"] == "A") & (self.df["Alt"] == "C"), "subtype"
        ] = "Transversion"
        self.df.loc[
            (self.df["Ref"] == "G") & (self.df["Alt"] == "T"), "subtype"
        ] = "Transversion"
        self.df.loc[
            (self.df["Ref"] == "G") & (self.df["Alt"] == "C"), "subtype"
        ] = "Transversion"

        # split snp / InDel
        self.df["group"] = "indel"
        self.df.loc[self.df.subtype.notnull(), "group"] = "snp"

        # InDel subtype: ins / del / complex
        self.df.loc[self.df["Ref"] == "-", "subtype"] = "Ins"
        self.df.loc[self.df["Alt"] == "-", "subtype"] = "Del"
        self.df.loc[self.df.subtype.isnull(), "subtype"] = "Complex"

    def annodf_statistics(self):
        # snp statistics
        self.df2 = self.df[self.df["group"] == "snp"]
        snp_count = self.df2.shape[0]
        # Ti / Tv
        snp_subtype = Counter(self.df2["subtype"])
        snp_funcdict = Counter(self.df2["func"])
        snp_gtdict = Counter(self.df2["GT"])

        snp_funcdict["intronic"] = snp_funcdict["intronic"] + snp_funcdict["intron"]
        snp_exondict = Counter(
            self.df2.loc[self.df["func"] == "exonic", "ExonicFunc.refGene"]
        )

        # basic stat
        self.snp_stat_dict["Basic_Statistics"]["SNP_Count"] = snp_count
        self.snp_stat_dict["Basic_Statistics"].update(snp_subtype)
        try:
            self.snp_stat_dict["Basic_Statistics"]["Ti / Tv"] = round(
                self.snp_stat_dict["Basic_Statistics"]["Transition"]
                / self.snp_stat_dict["Basic_Statistics"]["Transversion"],
                2,
            )
        except ZeroDivisionError:
            self.snp_stat_dict["Basic_Statistics"]["Ti / Tv"] = 0
        self.snp_stat_dict["Basic_Statistics"].update(snp_gtdict)
        # function stat
        self.snp_stat_dict["Function_Statistics"].update(snp_funcdict)
        try:
            self.snp_stat_dict["Function_Statistics"].pop("intron")
        except KeyError:
            pass
        self.snp_stat_dict["Exonic_Statistics"].update(snp_exondict)

        # InDel statistics
        self.df2 = self.df[self.df["group"] == "indel"]
        indel_count = self.df2.shape[0]
        indel_subtype = Counter(self.df2["subtype"])
        indel_funcdict = Counter(self.df2["func"])
        indel_funcdict["intronic"] = (
            indel_funcdict["intronic"] + indel_funcdict["intron"]
        )
        indel_exondict = Counter(
            self.df2.loc[self.df["func"] == "exonic", "ExonicFunc.refGene"]
        )
        indel_gtdict = Counter(self.df2["GT"])

        # basic stat
        self.indel_stat_dict["Basic_Statistics"]["InDel_Count"] = indel_count
        self.indel_stat_dict["Basic_Statistics"].update(indel_subtype)
        self.indel_stat_dict["Basic_Statistics"].update(indel_gtdict)
        # function stat
        self.indel_stat_dict["Function_Statistics"].update(indel_funcdict)
        try:
            self.indel_stat_dict["Function_Statistics"].pop("intron")
        except KeyError:
            pass
        # exonic stat
        self.indel_stat_dict["Exonic_Statistics"].update(indel_exondict)


class annovarStat(StepBase):
    def __init__(
        self,
        annoInput=None,  # list
        outputdir=None,
        threads=1,
        stepNum=None,
        upstream=None,
        **kwargs,
    ):

        super(annovarStat, self).__init__(stepNum, upstream)
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("annoInput")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)

            self.setParam("threads", threads)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())

        if (upstream is None) or (upstream is True):
            self.setInput("annoInput", annoInput)
            self.checkInputFilePath()
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            if upstream.__class__.__name__ == "annovar":
                self.setInput("annoInput", upstream.getOutput("txtOutput"))
            else:
                raise commonError("Parameter upstream must from annovar.")

        self.setOutput(
            "SNP-Basic_Statistics-Output",
            os.path.join(self.getOutput("outputdir"), "SNP-Basic-Stat.txt"),
        )

        self.setOutput(
            "SNP-Function_Statistics-Output",
            os.path.join(self.getOutput("outputdir"), "SNP-Func-Stat.txt"),
        )

        self.setOutput(
            "SNP-Exonic_Statistics-Output",
            os.path.join(self.getOutput("outputdir"), "SNP-Exonic-Stat.txt"),
        )

        self.setOutput(
            "InDel-Basic_Statistics-Output",
            os.path.join(self.getOutput("outputdir"), "InDel-Basic-Stat.txt"),
        )

        self.setOutput(
            "InDel-Function_Statistics-Output",
            os.path.join(self.getOutput("outputdir"), "InDel-Func-Stat.txt"),
        )

        self.setOutput(
            "InDel-Exonic_Statistics-Output",
            os.path.join(self.getOutput("outputdir"), "InDel-Exonic-Stat.txt"),
        )

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            overall_stat = {}
            for i in self.getInput("annoInput"):
                tmp_obj = stat_annofile(i)
                if tmp_obj.SampName not in overall_stat:
                    overall_stat.setdefault(tmp_obj.SampName, tmp_obj.statistics)
                else:
                    raise commonError(
                        "Sample %s from annoInput is repeat." % tmp_obj.SampName
                    )

            for i in ["SNP", "InDel"]:
                for j in [
                    "Basic_Statistics",
                    "Function_Statistics",
                    "Exonic_Statistics",
                ]:
                    tmp_dict = {}
                    for (samplename, values) in overall_stat.items():
                        tmp_dict.setdefault(samplename, values[i][j])
                        tmp_df = pd.DataFrame(tmp_dict).T
                        tmp_df.to_csv(
                            self.getOutput(f"{i}-{j}-Output"),
                            sep="\t",
                            index_label="Sample",
                            header=True,
                        )

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
