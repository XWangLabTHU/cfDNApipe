# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 11:37:46 2019

@author: zhang
"""


from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from collections import Iterable
import pysam, pybedtools, os, subprocess, sys
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import gzip, shutil

__metaclass__ = type


class commonError(Exception):
    def __init__(self, message):
        self.message = message
        

# https://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists
def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el
            
def isAlphaOrDigit(x):
    if x.isalpha() or x.isdigit():
        return True
    else:
        return False
    
def rmEndString(x, y):
    for item in y:
        if x.endswith(item):
            x = x.replace(item, '')
    
    return x


def isSoftClipped(cigar):
    """
    see here for more information about this function
    references:
        https://pysam.readthedocs.io/en/latest/api.html
        https://davetang.org/wiki/tiki-index.php?page=SAM
    """
    for (op,count) in cigar:
        if op in [4,5,6]: return True
    return False


def read_pair_generator(bam, region_string = None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    reference:
        https://www.biostars.org/p/306041/
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
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


def bamTobed(bamInput = None, bedOutput = None, compress = True):
    # generate temp file for sorting and indexing
    bedOutput_path = os.path.realpath(bedOutput)
    this_pid = os.getpid()
    tmp_split = os.path.splitext(bedOutput_path)
    tmp_bedOutput = tmp_split[0] + "-temp-" + str(this_pid) + tmp_split[1]
    
    bai = bamInput + ".bai"
    if not os.path.exists(bai):
        message = "Index file " + bai + " do not exist!"
        raise commonError(message)
    
    bedWrite = open(tmp_bedOutput, "w")
    
    input_file = pysam.Samfile(bamInput, "rb")
    chr_reference = input_file.references
    for read1, read2 in read_pair_generator(input_file):
        read1Start = read1.reference_start
        read1End = read1.reference_end
        read2Start = read2.reference_start
        read2End = read2.reference_end
        
        if not read1.is_reverse:  # read1 is forward strand, read2 is reverse strand
            rstart = read1Start  # 0-based left-most site
            rend = read2End
        else:  # read1 is reverse strand, read2 is forward strand
            rstart = read2Start  # 0-based left-most site
            rend = read1End
            
        if (rstart < 0) or (rend < 0) or (rstart >= rend): continue
    
        tmp_str = chr_reference[read1.tid] + "\t" + str(rstart) + "\t" +str(rend) + "\n"
        bedWrite.write(tmp_str)
    
    bedWrite.close()
    print("Fragments generated, waitting for sorting......")
    
    
    bedData = pybedtools.BedTool(tmp_bedOutput)
    bedData.sort(output=bedOutput)
    
    os.remove(tmp_bedOutput)
    
    print("Fragments sorted.")
    
    if compress:
        print("Waitting for compressing and indexing......")
        bedgzfile = bedOutput + ".gz"
        pysam.tabix_compress(bedOutput, bedgzfile, force = False)
        pysam.tabix_index(bedgzfile, preset = "bed", zerobased = True)
        print("Indexing bedgz file finished!")
    
    return True


def bamTobedForSingle(bamInput = None, bedOutput = None, compress = True):
    bam = pybedtools.BedTool(bamInput)
    bed = bam.bam_to_bed()
    bed.sort(output = bedOutput)
    
    print("Bed file generated!")
    
    if compress:
        print("Waitting for compressing and indexing......")
        bedgzfile = bedOutput + ".gz"
        pysam.tabix_compress(bedOutput, bedgzfile, force = False)
        pysam.tabix_index(bedgzfile, preset = "bed", zerobased = True)
        print("Indexing bedgz file finished!")
    
    return True


# bam to bed, tackle the file without "chr", this will raise error
def bam2bedV2(bamInput, bedOutput):
    bamInput = bamInput
    bedOutput = bedOutput
    
    # generate temp file for sorting and indexing
    bedOutput_path = os.path.realpath(bedOutput)
    this_pid = os.getpid()
    tmp_split = os.path.splitext(bedOutput_path)
    tmp_bedOutput = tmp_split[0] + "-temp-" + str(this_pid) + tmp_split[1]
    
    bedWrite = open(tmp_bedOutput, "w")
    
    bai = bamInput + ".bai"
    if not os.path.exists(bai):
        message = "Index file " + bai + " do not exist!"
        raise commonError(message)


    input_file = pysam.Samfile(bamInput, "rb")
    chr_reference = input_file.references
    for read1, read2 in read_pair_generator(input_file):
        read1Start = read1.reference_start
        read1End = read1.reference_end
        read2Start = read2.reference_start
        read2End = read2.reference_end
        
        if not read1.is_reverse:  # read1 is forward strand, read2 is reverse strand
            rstart = read1Start
            rend = read2End
        else:  # read1 is reverse strand, read2 is forward strand
            rstart = read2Start
            rend = read1End
            
        if (rstart < 0) or (rend < 0) or (rstart >= rend): continue
    
        tmp_str = "chr" + chr_reference[read1.tid] + "\t" + str(rstart) + "\t" +str(rend) + "\n"
        bedWrite.write(tmp_str)
    
    bedWrite.close()
    print("Fragments generated, waitting for sorting......")
    
    
    bedData = pybedtools.BedTool(tmp_bedOutput)
    bedData.sort(output=bedOutput)
    
    print("Fragments sorted, waitting for compressing and indexing......")
    
    bedgzfile = bedOutput + ".gz"
    pysam.tabix_compress(bedOutput, bedgzfile, force=False)
    pysam.tabix_index(bedgzfile, preset="bed", zerobased=True)
    
    print("Indexing bedgz file finished!")
    
    os.remove(tmp_bedOutput)
    
    return("Program complete!")


# plot length distribution
def fraglendistribution(bedInput = None, plotOutput = None, binOutput = None, maxLimit = None):
    data = pd.read_table(bedInput, sep="\t", header = None,
                         names = ["chr", "start", "end"])
    len_info = np.asarray(data['end'] - data['start'])
    
    np.save(binOutput, len_info)
    
    len_info = len_info[np.where(len_info <= maxLimit)]
    
    fig = plt.figure(figsize = (10, 8))
    plot_limit = len_info.max() - len_info.min()
    plt.hist(len_info, bins = plot_limit)
    plt.savefig(plotOutput)
    plt.close(fig)
    
    return True
    
    
# calculate methylation level for regions
def calcMethyl(bamInput, bedInput, txtOutput):
    bai = bamInput + ".bai"
    if not os.path.exists(bai):
        message = "Index file " + bai + " do not exist!"
        raise commonError(message)
    
    bam_input = pysam.Samfile(bamInput, "rb")
    regions = pd.read_csv(bedInput, sep = "\t", header = None, names = ["chr", "start", "end"])
    
    
    CXXname = ["unmCpG", "mCpG",
               "unmCHG", "mCHG",
               "unmCHH", "mCHH",
               "unmUNC", "mUNC"]
    d = dict.fromkeys(CXXname, 0)
    regions = regions.assign(**d)
    
    for index, row in regions.iterrows():
        count_data = [0, 0, 0, 0, 0, 0, 0, 0]
        for read in bam_input.fetch(reference = row["chr"], start = row["start"], end = row["end"]):
            CXXinfo = read.get_tag('XM')
            count_data[0] += CXXinfo.count("z")
            count_data[1] += CXXinfo.count("Z")
            count_data[2] += CXXinfo.count("x")
            count_data[3] += CXXinfo.count("X")
            count_data[4] += CXXinfo.count("h")
            count_data[5] += CXXinfo.count("H")
            count_data[6] += CXXinfo.count("u")
            count_data[7] += CXXinfo.count("U")
        
        regions.loc[index, CXXname] = count_data
    
    regions["mlCpG"] = regions["mCpG"] / (regions["mCpG"] + regions["unmCpG"])
    regions["mlCHG"] = regions["mCHG"] / (regions["mCHG"] + regions["unmCHG"])
    regions["mlCHH"] = regions["mCHH"] / (regions["mCHH"] + regions["unmCHH"])
    regions["mlUNC"] = regions["mUNC"] / (regions["mUNC"] + regions["unmUNC"])
    
    regions.to_csv(txtOutput, sep = "\t", header = True, index = False)


# uncompress gz file
def un_gz(gzfile):
    file = gzfile.replace(".gz", "")
    with gzip.open(gzfile, 'r') as f_in, open(file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


# run a single command line
def cmdCall(cmdLine):
    proc = subprocess.Popen(cmdLine, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, universal_newlines = True)
    while True:
        nextline = proc.stdout.readline()
        if (nextline == '') and (proc.poll() is not None):
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()
        
    output, error = proc.communicate()
    exitCode = proc.returncode

    if exitCode != 0:
        raise commonError('**********CMD running error**********')


# compute OCF value for paired end data
def ComputeOCF(bedgz, txtOutput, OCFOutput, regionFile):
    print("Input file:", bedgz)
    print("Output files:", txtOutput)
    print("Output OCF files:", OCFOutput)
    print("Region files:", regionFile)
    
    tbx = pysam.TabixFile(bedgz)
    
    # regions are 1-based
    regions = pd.read_csv(regionFile, sep = "\t", header = None, names = ["chr", "start", "end", "description"])
    
    regions["OCF"] = 0
    cud_output = txtOutput
    
    Udata = np.empty(shape=(0, 241))
    Ddata = np.empty(shape=(0, 241))
    
    for idx, region in regions.iterrows():
        region_Chr, region_Start, region_End = region["chr"], region["start"], region["end"]
        
        if region_Start < 1:
            message = "Start of the region must > 0!"
            raise commonError(message)
            
        covPOS = defaultdict(lambda:[0, 0, 0])
        
        # fetch read in the region
        try:
            fetched_reads = tbx.fetch(region_Chr, region_Start, region_End)
        except ValueError as e:
            continue
        for row in fetched_reads:
            tmp_row = row.split()
            rstart = int(tmp_row[1]) + 1 # convert to 1-based
            rend = int(tmp_row[2]) # end included
            
#            for i in range(rstart, rend + 1):  # for a read, add 1 to every point at which it overlapped, coverage
#                if i >= region_Start and i <= region_End:
#                    covPOS[i][0] += 1
            
            if rstart >= region_Start and rstart <= region_End:  # consider read start point, U end
                covPOS[rstart][1] += 1
                
            if rend >= region_Start and rend <= region_End:  # consider read start point, D end
                covPOS[rend][2] += 1
        
        # finished a region
        midpoint = int((region_Start + region_End) / 2)
        
        Udata = np.vstack((Udata, [covPOS[x][1] for x in range(midpoint - 120, midpoint + 121)]))
        Ddata = np.vstack((Ddata, [covPOS[x][2] for x in range(midpoint - 120, midpoint + 121)]))
        
        left_OCF = sum([covPOS[x][2] for x in range(midpoint - 70, midpoint - 50)]) - sum([covPOS[x][1] for x in range(midpoint - 70, midpoint - 50)])
        right_OCF = sum([covPOS[x][1] for x in range(midpoint + 50, midpoint + 70)]) - sum([covPOS[x][2] for x in range(midpoint + 50, midpoint + 70)])
        regions.loc[idx, "OCF"] = left_OCF + right_OCF
    
    ud_data = np.concatenate([Udata, Ddata], 1)
    
    np.savetxt(cud_output, ud_data, fmt = "%i", delimiter = '\t')
    
    regions.to_csv(OCFOutput, sep = "\t", index = False)
    
    print("Processing finished!")
    
    return("True")


# compute coverage, U-end(upstream end) and D-end(downstream end)
# only count -1000 to 1000 from open region center
def computeCUE(inputFile, refFile, txtOutput, cudOutput):
    inputFile = pybedtools.BedTool(inputFile)
    refFile = pybedtools.BedTool(refFile)
    inputFile.intersect(refFile, wo = True, sorted = True, output = txtOutput)
    
    peak = 60 #distance of U and D peaks from the center
    bin = 10 #half the length of windows
    ocf = []
  
    data = pd.read_csv(txtOutput, sep = "\t", header = None, names = ["read.chr", "read.start", "read.end", "peak.chr", 
                                                                      "peak.start", "peak.end", "description", "overlap"])
    data["peak.start"] = data["peak.start"] + 1
    data["read.start"] = data["read.start"] + 1
    
    save_flag = ["Tcell", "Liver", "Placenta", "Lung", "Breast", "Intestine", "Ovary"]
    flag_num = -1
    for flag in save_flag:
        flag_num += 1
        print("Now, processing " + flag + "......")
        tmp_data= data.loc[data.description == flag, ]
        cov = np.zeros(shape = 2000)
        uend = np.zeros(shape = 2000)
        dend = np.zeros(shape = 2000)
        for idx, row in tmp_data.iterrows():
            if (row["read.start"] < row["peak.start"]) and (row["peak.start"] <= row["read.end"]):
                o_s, o_e = 0, row["read.end"] - row["peak.start"]
                cov[0 : (o_e + 1)] += 1
                dend[o_e] += 1
            elif (row["peak.start"] <= row["read.start"]) and (row["read.end"] <= row["peak.end"]):
                o_s, o_e = row["read.start"] - row["peak.start"], row["read.end"] - row["peak.start"]
                cov[o_s : (o_e + 1)] += 1
                uend[o_s] += 1
                dend[o_e] += 1
            elif (row["read.start"] <= row["peak.end"]) and (row["peak.end"] < row["read.end"]):
                o_s, o_e = row["read.start"] - row["peak.start"], 1999
                cov[o_s : (o_e + 1)] += 1
                uend[o_s] += 1
            else:
                continue
                
        cov_tot, uend_tot, dend_tot = 0, 0, 0
        for i in range(2000):
            cov_tot += cov[i]
            uend_tot += uend[i]
            dend_tot += dend[i]
            
        index = np.array(range(-1000, 1000))
        df = pd.DataFrame.from_dict({'idx': index, 'cov': cov, 'cov%%': cov / cov_tot * 10000, 'uend': uend, 'uend%%': uend / uend_tot * 10000,'dend': dend,'dend%%': dend / dend_tot * 10000})
        df.to_csv(cudOutput[flag_num], sep = "\t", index = False)
        
        trueends = 0
        background = 0
        for i in index:
            if i >= -peak-bin and i <= -peak+bin:
                trueends += dend[i] / dend_tot * 10000
                background += uend[i] / uend_tot * 10000
            elif i >= peak-bin and i <= peak+bin:
                trueends += uend[i] / uend_tot * 10000
                background += dend[i] / dend_tot * 10000
        ocf.append(trueends - background)
        
        print("Processing " + flag + " finished!")
    
    return(ocf)

def OCFplot(ocfcaseinput, ocfctrlinput, output, x_label = ['case', 'control']):
    save_flag = ["Tcell", "Liver", "Placenta", "Lung", "Breast", "Intestine", "Ovary"]
    ocfcaseinput = np.transpose(ocfcaseinput).tolist()
    ocfctrlinput = np.transpose(ocfctrlinput).tolist()
    
    plt.figure()
    bpl = plt.boxplot(ocfcaseinput, positions = np.array(range(len(ocfcaseinput))) * 2.0 - 0.4, sym = '', widths = 0.6, vert = True)
    bpr = plt.boxplot(ocfctrlinput, positions = np.array(range(len(ocfctrlinput))) * 2.0 + 0.4, sym = '', widths = 0.6, vert = True)
    plt.setp(bpl['boxes'], color = 'y')
    plt.setp(bpl['whiskers'], color = 'y')
    plt.setp(bpl['caps'], color = 'y')
    plt.setp(bpl['medians'], color = 'y')
    plt.setp(bpr['boxes'], color = 'b')
    plt.setp(bpr['whiskers'], color = 'b')
    plt.setp(bpr['caps'], color = 'b')
    plt.setp(bpr['medians'], color = 'b')
    plt.plot([], c = 'y', label = x_label[0])
    plt.plot([], c = 'b', label = x_label[1])
    plt.legend()
    plt.xticks(range(0, len(save_flag) * 2, 2), save_flag)
    plt.xlim(-2, len(save_flag) * 2)
    for k in range(7):
        plt.scatter([k * 2.0 - 0.4 for j in range(len(ocfcaseinput[k]))], ocfcaseinput[k], s = 8, c = 'y')
        plt.scatter([k * 2.0 + 0.4 for j in range(len(ocfctrlinput[k]))], ocfctrlinput[k], s = 8, c = 'b')
    plt.ylabel('OCF value')
    plt.tight_layout()
    plt.savefig(output)
    
    return(True)
    
def generate_cudoutput(input, outputdir):
    save_flag = ["Tcell", "Liver", "Placenta", "Lung", "Breast", "Intestine", "Ovary"]
    dict = {}
    prefix = os.path.splitext(os.path.basename(input))[0]
    for flag in save_flag:
        dict[flag] = outputdir + '/' + prefix + '-' + flag + '-cud.txt'
    return dict

# compute WPS(window protection score)
def computeWPS(minInsSize, maxInsSize, protection, outfile, input_region, input_frag, empty):
#    minInsSize, maxInsSize = None, None
    if minInsSize > 0 and maxInsSize > 0 and minInsSize < maxInsSize:
        minInsSize = minInsSize
        maxInsSize = maxInsSize
    
    print(minInsSize, maxInsSize)
    outfile = outfile.strip("""\'""")
    protection = protection//2
    
    validChroms = set(map(str,list(range(1,23))+["X","Y"]))  # human genome
    
    infile = open(input_region)  # input transcript region file
    
    bedgzfile = input_frag
    bedgzfile = bedgzfile.strip("""\'""")
    tbx = pysam.TabixFile(bedgzfile)
    prefix = "chr"
    
    for line in infile.readlines():
        cid, chrom, start, end, strand = line.split() # positions should be 0-based and end non-inclusive
        chrom = chrom.replace("chr","")
        if chrom not in validChroms: continue
        regionStart, regionEnd = int(start), int(end)  # this is 1-based
        
        if regionStart < 1: continue  # invalid region

        posRange = defaultdict(lambda:[0,0])
        filteredReads = Intersecter()
        try:  # if tbx.fetch do not find any row, next row
            for row in tbx.fetch(prefix+chrom, regionStart-protection-1, regionEnd+protection+1):  # all fragments overlaped with this region is collected
                tmp_row = row.split()
                
                rstart = int(tmp_row[1]) + 1 # convert to 1-based
                rend = int(tmp_row[2]) # end included
                lseq = int(tmp_row[2]) - int(tmp_row[1]) # fragment length
                
                if (minInsSize != None) and (maxInsSize != None) and ((lseq < minInsSize) or (lseq > maxInsSize)): continue  # satisfy length requirement
                
                filteredReads.add_interval(Interval(rstart, rend))  # save the fragments overlap with region
                
                for i in range(rstart, rend+1):  # for a single nucleotide site, compute how many reads overlaped span it (include read end point)
                    if i >= regionStart and i <= regionEnd:
                        posRange[i][0] += 1
                
                if rstart >= regionStart and rstart <= regionEnd:  # for a single nucleotide site, compute how many read end point located at this site
                    posRange[rstart][1] += 1
                    
                if rend >= regionStart and rend <= regionEnd:
                    posRange[rend][1] += 1
        except:
            continue
        
        filename = outfile%cid  # name output file by names in transcript file
        outfile = gzip.open(filename,'w')
        cov_sites = 0
        outLines = []
        for pos in range(regionStart,regionEnd+1):
            rstart, rend = pos - protection, pos + protection
            gcount, bcount = 0, 0
            for read in filteredReads.find(rstart,rend):
                if (read.start > rstart) or (read.end < rend): bcount += 1  # fragments located in window
                else: gcount += 1  # fragments spanned window
            covCount,startCount = posRange[pos]
            cov_sites += covCount
            # chrom: chromatin, pos: position in the genome, covCount:how many reads span this site, startCount: how many reads end point located
            # in this site, gcount-bcount: WPS
            outLines.append("%s\t%d\t%d\t%d\t%d\n"%(chrom, pos, covCount, startCount, gcount-bcount))
        
        if strand == "-": outLines = outLines[::-1]  # - strand!!!
        for line in outLines: outfile.write(line.encode())  # write in binary
        outfile.close()
        
        if cov_sites == 0 and not empty:  # remove empty files
            os.remove(filename)
    
    return("True")


# compress .bismark.zero.cov file
def compressMethy(InputFile = None):
    '''
    input must from bismark_methylation_extractor, .bismark.zero.cov file
    '''
    bedgzfile = InputFile + ".gz"
    pysam.tabix_compress(InputFile, bedgzfile, force = False)
    pysam.tabix_index(bedgzfile, preset="bed", zerobased = True)
    
    return("True")



# compress methylation level from .bismark.zero.cov.gz file
def calcMethylV2(tbxInput, bedInput, txtOutput):
    tbi = tbxInput + ".tbi"
    if not os.path.exists(tbi):
        message = "Index file " + tbi + " do not exist!"
        raise commonError(message)
    
    tbx_input = pysam.TabixFile(tbxInput)
    
    regions = pd.read_csv(bedInput, sep = "\t", header = None, names = ["chr", "start", "end"])
    
    CXXname = ["unmCpG", "mCpG"]
    d = dict.fromkeys(CXXname, 0)
    regions = regions.assign(**d)
    
    print("Now, processing fetch and computing CpG level.")
    
    for index, row in regions.iterrows():
        count_data = [0, 0]
        for read in tbx_input.fetch(reference = row["chr"], start = row["start"], end = row["end"]):
            readinfo = read.split()
            count_data[0] += int(readinfo[5])
            count_data[1] += int(readinfo[4])
        
        regions.loc[index, CXXname] = count_data
    
    regions["mlCpG"] = regions["mCpG"] / (regions["mCpG"] + regions["unmCpG"])
    
    regions = regions.fillna(0)
    
    regions.to_csv(txtOutput, sep = "\t", header = True, index = False)
    
    print("finished!")
    
    return(True)




















