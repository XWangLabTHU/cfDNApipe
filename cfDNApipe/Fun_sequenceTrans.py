# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 18:27:32 2019

@author: Jiqi Huang, Chang Li
"""


from .StepBase import StepBase
from .cfDNA_utils import commonError, isSoftClipped
from .Configure import Configure
from collections import defaultdict
import pysam, os


__metaclass__ = type


class sequencetransfer(StepBase):
    def __init__(self, 
         bamInput = None, # list
         bedInput = None,
         outputdir = None, # str
         threads = 1,
         upstream = None,
         formerrun = None,
         **kwargs):
        bedflag = False
        
        if upstream is None:
            super(sequencetransfer, self).__init__()
            self.setInput('bamInput', bamInput)
            
            if bedInput is not None:
                bedflag = True
                self.setInput('bedInput', bedInput)
                
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput('outputdir', os.path.dirname(os.path.abspath(self.getInput('bamInput')[0])))
            else:
                self.setOutput('outputdir', outputdir)
            
            self.setParam('threads', threads)
        else:
            if formerrun is None:
                super(sequencetransfer, self).__init__(upstream.getStepID())
            else:
                super(sequencetransfer, self).__init__(formerrun.getStepID())
                
            Configure.configureCheck()
            upstream.checkFilePath()
            
            if upstream.__class__.__name__ == 'rmduplicate':
                self.setInput('bamInput', upstream.getOutput('bamOutput'))
            else:
                raise commonError('Parameter upstream must from inputprocess or adapterremoval.')
            
            if bedInput is not None:
                bedflag = True
                self.setInput('bedInput', bedInput)
                
            self.setOutput('outputdir', self.getStepFolderPath())
            self.setParam('threads', Configure.getThreads())
        
        self.setOutput('txtOutput', [os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x)) + '-seq.txt' for x in self.getInput('bamInput')])
        
        if bedflag:
            bed_input = open(self.getInput('bedInput'), 'r')
            bed_lines = bed_input.readlines()
            
        multi_run_len = len(self.getInput('bamInput'))
        for i in range(multi_run_len):
            bam_input = pysam.AlignmentFile(self.getInput('bamInput')[i], 'rb')
            txt_output = open(self.getOutput('txtOutput')[i], 'w')
            if bedflag:
                for j in range(len(bed_lines)):
                    bed_line = bed_lines[j]
                    bed_intv = bed_line.strip().split('\t')
                    pairs = self.read_pair_generator_inputref(bam = bam_input, refs = bed_intv)
                    txt_output.write('Sequences in ' + bed_intv[0] + ', ' + bed_intv[1] + ' - ' + bed_intv[2] + ':\n')
                    self.transfer_output(pairs, txt_output)
                    txt_output.write('--------------------------\n\n')
            else:
                pairs = self.read_pair_generator_inputref(bam = bam_input)
                self.transfer_output(pairs, txt_output)
            
            bam_input.close()
            txt_output.close()
            
        if bedflag:
            bed_input.close()
        
        finishFlag = self.stepInit(upstream)
        
        self.excute(finishFlag, runFlag = False)
            
    def read_pair_generator_inputref(bam, refs = None):
        """
        Generate read pairs in a BAM file or within a reference list split from one line in BED file.
        Reads are added to read_dict until a pair is found.
        reference:
            https://www.biostars.org/p/306041/
        """
        read_dict = defaultdict(lambda: [None, None])
        for read in bam.fetch(reference = refs[0], start = int(refs[1]), end = int(refs[2])):
            # filter reads
            if read.is_unmapped or read.is_qcfail or read.is_duplicate: continue
            if not read.is_paired: continue
            if not read.is_proper_pair: continue
            if read.is_secondary or read.is_supplementary: continue
            if read.mate_is_unmapped: continue
            if read.rnext != read.tid: continue
            if read.template_length == 0: continue
            if isSoftClipped(read.cigar): continue
            
            qname = read.query_name
            if qname not in read_dict:
                if read.is_read1:
                    read_dict[qname][0] = read
                else:
                    read_dict[qname][1] = read
            else:
                if read.is_read1:
                    yield read, read_dict[qname][1]
                else:
                    yield read_dict[qname][0], read
                del read_dict[qname]
    
    def transfer_output(pairs, file_output):
        for r in pairs:
            if r[0].is_read1:
                seq1 = ''
                seq2 = ''
                
                #mark unmethylated C/G
                if r[0].get_tag('XG') == 'GA':
                    for i in range(len(r[0].query_sequence)):
                        if r[0].get_tag('XM')[i] >= 'a' and r[0].get_tag('XM')[i] <= 'z':
                            seq1 += 'g'
                        else:
                            seq1 += r[0].query_sequence[i]
                elif r[0].get_tag('XG') == 'CT':
                    for i in range(len(r[0].query_sequence)):
                        if r[0].get_tag('XM')[i] >= 'a' and r[0].get_tag('XM')[i] <= 'z':
                            seq1 += 'c'
                        else:
                            seq1 += r[0].query_sequence[i]
                if r[1].get_tag('XG') == 'GA':
                    for i in range(len(r[1].query_sequence)):
                        if r[1].get_tag('XM')[i] >= 'a' and r[1].get_tag('XM')[i] <= 'z':
                            seq2 += 'g'
                        else:
                            seq2 += r[1].query_sequence[i]
                elif r[1].get_tag('XG') == 'CT':
                    for i in range(len(r[1].query_sequence)):
                        if r[1].get_tag('XM')[i] >= 'a' and r[1].get_tag('XM')[i] <= 'z':
                            seq2 += 'c'
                        else:
                            seq2 += r[1].query_sequence[i]
                            
                #process if is reverse
                if r[0].is_reverse:
                    for i in range(len(seq1)):
                        if seq1[i] == 'T':
                            seq1[i] == 'A'
                        elif seq1[i] == 'A':
                            seq1[i] == 'T'
                        elif seq1[i] == 'C':
                            seq1[i] == 'G'
                        elif seq1[i] == 'G':
                            seq1[i] == 'C'
                        elif seq1[i] == 'c':
                            seq1[i] == 'g'
                        elif seq1[i] == 'g':
                            seq1[i] == 'c'
                if r[1].is_reverse:
                    for i in range(len(seq2)):
                        if seq2[i] == 'T':
                            seq2[i] == 'A'
                        elif seq2[i] == 'A':
                            seq2[i] == 'T'
                        elif seq2[i] == 'C':
                            seq2[i] == 'G'
                        elif seq2[i] == 'G':
                            seq2[i] == 'C'
                        elif seq2[i] == 'c':
                            seq2[i] == 'g'
                        elif seq2[i] == 'g':
                            seq2[i] == 'c'
                
                offset = 0  #the offset of m due to the insertion or deletion in r
                if r[0].cigartuples != [(0, 75)]:
                    for i in range(len(r[0].cigartuples)):
                        if r[0].cigartuples[i][0] == 1:    #case of insertion; offset is set to positive, which means that m should move right
                            offset += r[0].cigartuples[i][1]
                        elif r[0].cigartuples[i][0] == 2:    #case of deletion; offset is set to negative, which means that m should move left
                            offset -= r[0].cigartuples[i][1]
                
                if abs(r[0].template_length) > 150 - offset:
                    file_output.write(seq1 + '\n')
                    file_output.write(seq2 + '\n\n')
                elif abs(r[0].template_length) <= 150 - offset and r[0].template_length > 75 - offset:
                    seq2_cut = seq2[150 - r[0].template_length - offset : ]
                    file_output.write(seq1 + seq2_cut + '\n\n')
                elif abs(r[0].template_length) < 75 - offset:
                    file_output.write(seq1 + '\n\n')