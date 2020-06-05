# -*- coding: utf-8 -*-
"""
Created on Tue Mar 3 09:37:21 2020

@author: He Shuying

"""


from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
import math
from .Configure import Configure


__metaclass__ = type


class virusbismark(StepBase):
    def __init__(
        self,
        seqInput1=None,  # list
        seqInput2=None,  # list
        ref=None,  # str, virus ref
        OutputDir=None,  # str
        threads=1,
        paired=True,
        other_params={
            "-q": True,
            "--phred33-quals": True,
            "--bowtie2": True,
            "--un": True,
        },
        stepNum=None,
        upstream=None,
        **kwargs
    ):
        """
        do not use prefix paramter in bismark
        """
        super(virusbismark, self).__init__(stepNum, upstream)
        if (upstream is None) or (upstream is True):
            self.setInput("seq1", seqInput1)

            if paired:
                self.setParam("type", "paired")
                self.setInput("seq2", seqInput2)
            else:
                self.setParam("type", "single")
        else:
            Configure.configureCheck()
            upstream.checkFilePath()

            self.setParam("type", Configure.getType())

            if upstream.__class__.__name__ == "lowcomplexityfilter":
                if self.getParam("type") == "paired":
                    self.setInput("seq1", upstream.getOutput("seq1Output"))
                    self.setInput("seq2", upstream.getOutput("seq2Output"))
                elif self.getParam("type") == "single":
                    self.setInput("seq1", upstream.getOutput("seq1Output"))
                else:
                    commonError("Wrong data tpye, must be 'single' or 'paired'!")
            else:
                raise commonError("Parameter upstream must from lowcomplexityfilter.")
            self.checkInputFilePath()

        if upstream is None:
            if OutputDir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("seq1")[0])),
                )
            else:
                self.setOutput("outputdir", OutputDir)
            self.setParam("threads", threads)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())
            self.setParam("threads", Configure.getThreads())

        # check virus reference for bismark
        self.setParam("ref", ref)
        self.bismkrefcheck()

        self.setParam(
            "prefix",
            [
                self.getMaxFileNamePrefixV2(x) + "_bismark_virus"
                for x in self.getInput("seq1")
            ],
        )

        self.setParam(
            "outPrefix",
            [
                os.path.join(self.getOutput("outputdir"), x)
                for x in self.getParam("prefix")
            ],
        )

        if other_params is None:
            self.setParam("other_params", "")
        else:
            self.setParam("other_params", other_params)

        self.setOutput(
            "unmapped-1",
            [x + "_unmapped_reads_1.fq.gz" for x in self.getParam("outPrefix")],
        )

        self.setOutput(
            "bamOutput", [x + "_pe.bam" for x in self.getParam("outPrefix")],
        )

        self.setOutput(
            "bismkRepOutput",
            [x + "_PE_report.txt" for x in self.getParam("outPrefix")],
        )

        if self.getParam("type") == "paired":
            self.setOutput(
                "unmapped-2",
                [x + "_unmapped_reads_2.fq.gz" for x in self.getParam("outPrefix")],
            )

        if len(self.getInput("seq1")) == len(self.getInput("seq2")):
            multi_run_len = len(self.getInput("seq1"))
        else:
            raise commonError("Paired end Input files are not consistent.")

        all_cmd = []

        for i in range(multi_run_len):
            tmp_cmd = [
                "bismark",
                self.getParam("other_params"),
                "-p",
                math.ceil(self.getParam("threads") / 5),
                "-B",
                self.getParam("prefix")[i],
                "--output_dir",
                self.getOutput("outputdir"),
                "--temp_dir",
                self.getOutput("outputdir"),
                "--genome_folder",
                self.getParam("ref"),
            ]
            if self.getParam("type") == "paired":
                tmp_cmd.extend(
                    ["-1", self.getInput("seq1")[i], "-2", self.getInput("seq2")[i]]
                )
            elif self.getParam("type") == "single":
                tmp_cmd.extend([self.getInput("seq1")[i]])
            else:
                raise commonError("Wrong data tpye, must be 'single' or 'paired'!")

            all_cmd.append(self.cmdCreate(tmp_cmd))

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            self.run(all_cmd)

        self.stepInfoRec(cmds=all_cmd, finishFlag=finishFlag)

    # ref check

    def bismkrefcheck(self,):
        fafile = [os.path.join(self.getParam("ref"), "viral_REFSEQ.fa")]
        CTfiles = [
            os.path.join(self.getParam("ref"), "Bisulfite_Genome/CT_conversion/" + x)
            for x in [
                "BS_CT.1.bt2",
                "BS_CT.2.bt2",
                "BS_CT.3.bt2",
                "BS_CT.4.bt2",
                "BS_CT.rev.1.bt2",
                "BS_CT.rev.2.bt2",
                "genome_mfa.CT_conversion.fa",
            ]
        ]
        BAfiles = [
            os.path.join(self.getParam("ref"), "Bisulfite_Genome/GA_conversion/" + x)
            for x in [
                "BS_GA.1.bt2",
                "BS_GA.2.bt2",
                "BS_GA.3.bt2",
                "BS_GA.4.bt2",
                "BS_GA.rev.1.bt2",
                "BS_GA.rev.2.bt2",
                "genome_mfa.GA_conversion.fa",
            ]
        ]
        #    bismkRef = CTfiles + BAfiles
        bismkRef = fafile + CTfiles + BAfiles
        for filePath in bismkRef:
            if not os.path.exists(filePath):
                raise commonError("Bowtie2 index file " + filePath + " don not exist!")
