# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 09:37:21 2019

@author: Wei Zhang

E-mail: w-zhang16@mails.tsinghua.edu.cn
"""


from StepBase import StepBase
from cfDNA_utils import commonError
import os, re
from Configure import Configure


__metaclass__ = type


class bismark(StepBase):
    def __init__(self, 
             seqInput1 = None, # list
             seqInput2 = None, # list
             ref = None, # str
             outputdir = None, # str
             threads = 1,
             other_params = {'-q': True, '--phred33-quals': True, '-N': 1, '-X': 2000, 
                             '--bowtie2': True, '--no_dovetail': True},
             upstream = None,
             formerrun = None,
             **kwargs):
        if upstream is None:
            super(bismark, self).__init__()
            self.setInput('seq1', seqInput1)
            self.setInput('seq2', seqInput2)
            self.checkInputFilePath()
            
            self.setParam('ref', ref)
            
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('fq1')[1])))
            else:
                self.setOutput('outputdir', outputdir)
                
            self.setParam('threads', threads)
            
        else:
            # here to check upstream source!!!!!!!!!!!!!!!!!!!!!!!!!!
            # upstream can come from inputprocess and adapterremoval
            
            if formerrun is None:
                super(bismark, self).__init__(upstream.getStepID())
            else:
                super(bismark, self).__init__(formerrun.getStepID())
                
            # check Configure for running pipeline
            Configure.configureCheck()
            
            upstream.checkFilePath()
            
            if upstream.__class__.__name__ == 'inputprocess':
                self.setInput('seq1', upstream.getOutput('fq1'))
                self.setInput('seq2', upstream.getOutput('fq2'))
            elif upstream.__class__.__name__ == 'adapterremoval':
                self.setInput('seq1', upstream.getOutput('pair1'))
                self.setInput('seq2', upstream.getOutput('pair2'))
            else:
                raise commonError('Parameter upstream must from inputprocess or adapterremoval.')
            
            self.setParam('ref', Configure.getRefDir())
            self.setOutput('outputdir', self.getStepFolderPath())
            self.setParam('threads', Configure.getThreads())
        
        # check reference for bismark
        self.bismkrefcheck()
        
        # generate base name
        prefix = []
        for seq1, seq2 in zip(self.getInput('seq1'), self.getInput('seq2')):
            prefix.append(self.getMaxFileNamePrefix(seq1, seq2))
        
        self.setParam('prefix', prefix)
        
        self.setParam('outPrefix', [os.path.join(self.getOutput('outputdir'), x) for x in self.getParam('prefix')],)
        
        if other_params is None:
            self.setParam('other_params', '')
        else:
            self.setParam('other_params',  other_params)
            
        self.setOutput('unmapped-1', [x + '_unmapped_reads_1.fq.gz' for x in self.getParam('outPrefix')])
        self.setOutput('unmapped-2', [x + '_unmapped_reads_2.fq.gz' for x in self.getParam('outPrefix')])
        self.setOutput('bamOutput', [x + '_pe.bam' for x in self.getParam('outPrefix')])
        self.setOutput('bismkRepOutput', [x + '_PE_report.txt' for x in self.getParam('outPrefix')])
            
        if len(self.getInput('seq1')) == len(self.getInput('seq1')):
            multi_run_len = len(self.getInput('seq1'))
        else:
            raise commonError('Paired end Input files are not consistent.')
        
        all_cmd = []
        
        for i in range(multi_run_len):
            tmp_cmd = self.cmdCreate(["bismark", 
                                       '--genome_folder', self.getParam('ref'),
                                       self.getParam('other_params'),
                                       '--unmapped', self.getParam('prefix')[i],
                                       '--basename', self.getParam('prefix')[i],
                                       '-p', self.getParam('threads'),
                                       '--output_dir', self.getOutput('outputdir'),
                                       '-1', self.getInput('seq1')[i],
                                       '-2', self.getInput('seq2')[i]])
            all_cmd.append(tmp_cmd)
        
        self.setParam('cmd', all_cmd)
        
        finishFlag = self.stepInit(upstream)
        
        self.excute(finishFlag)


    # ref check
    def bismkrefcheck(self, ):
        fafile = [os.path.join(self.getParam('ref'), Configure.getGenome() + '.fa')]
        CTfiles = [os.path.join(self.getParam('ref'), 'Bisulfite_Genome/CT_conversion/' + x) for x in ['BS_CT.1.bt2', 'BS_CT.2.bt2', 'BS_CT.3.bt2', 'BS_CT.4.bt2', 'BS_CT.rev.1.bt2', 'BS_CT.rev.2.bt2', 'genome_mfa.CT_conversion.fa']]
        BAfiles = [os.path.join(self.getParam('ref'), 'Bisulfite_Genome/GA_conversion/' + x) for x in ['BS_GA.1.bt2', 'BS_GA.2.bt2', 'BS_GA.3.bt2', 'BS_GA.4.bt2', 'BS_GA.rev.1.bt2', 'BS_GA.rev.2.bt2', 'genome_mfa.GA_conversion.fa']]
        bismkRef = fafile + CTfiles + BAfiles
        for filePath in bismkRef:
            if not os.path.exists(filePath):
                raise commonError('Bowtie2 index file ' + filePath + ' don not exist!')

