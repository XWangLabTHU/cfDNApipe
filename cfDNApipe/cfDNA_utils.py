# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 11:37:46 2019

@author: zhang
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from bx.intervals.intersection import Intersecter, Interval
from collections import Iterable
import pysam
import pybedtools
import os
import subprocess
import sys
import random
from collections import defaultdict
import pandas as pd
import numpy as np

# import matplotlib.pyplot as plt
import gzip
import shutil
import statsmodels.api as sm
import statsmodels.stats.multitest as multi
from scipy import stats, optimize, signal
from scipy.interpolate import interp1d
from sklearn.decomposition import PCA
from sklearn.svm import NuSVR
import seaborn as sns
from tqdm import tqdm
from sklearn.model_selection import KFold
from sklearn.model_selection import LeaveOneOut
from PIL import Image, ImageFont, ImageDraw
import pkg_resources
import collections
import pickle


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
            x = x.replace(item, "")

    return x


def isSoftClipped(cigar):
    """
    see here for more information about this function
    references:
        https://pysam.readthedocs.io/en/latest/api.html
        https://davetang.org/wiki/tiki-index.php?page=SAM
    """
    for (op, count) in cigar:
        if op in [4, 5, 6]:
            return True
    return False


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    reference:
        https://www.biostars.org/p/306041/
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        # filter reads
        if read.is_unmapped or read.is_qcfail or read.is_duplicate:
            continue
        if not read.is_paired:
            continue
        if not read.is_proper_pair:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if read.mate_is_unmapped:
            continue
        if read.rnext != read.tid:
            continue
        if read.template_length == 0:
            continue
        if isSoftClipped(read.cigar):
            continue

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


def bamTobed(bamInput=None, bedOutput=None, fragFilter=False, minLen=None, maxLen=None):
    compress = True
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

        if (rstart < 0) or (rend < 0) or (rstart >= rend):
            continue

        if fragFilter:
            readLen = rend - rstart
            if (readLen <= minLen) or (readLen >= maxLen):
                continue

        tmp_str = chr_reference[read1.tid] + "\t" + str(rstart) + "\t" + str(rend) + "\n"
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
        pysam.tabix_compress(bedOutput, bedgzfile, force=False)
        pysam.tabix_index(bedgzfile, preset="bed", zerobased=True)
        print("Indexing bedgz file finished!")

    return True


def bamTobedForSingle(bamInput=None, bedOutput=None):
    compress = True
    bam = pybedtools.BedTool(bamInput)
    bed = bam.bam_to_bed()
    bed.sort(output=bedOutput)

    print("Bed file generated!")

    if compress:
        print("Waitting for compressing and indexing......")
        bedgzfile = bedOutput + ".gz"
        pysam.tabix_compress(bedOutput, bedgzfile, force=False)
        pysam.tabix_index(bedgzfile, preset="bed", zerobased=True)
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

        if (rstart < 0) or (rend < 0) or (rstart >= rend):
            continue

        tmp_str = "chr" + chr_reference[read1.tid] + "\t" + str(rstart) + "\t" + str(rend) + "\n"
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

    return "Program complete!"


# plot length distribution
def fraglendistribution(bedInput=None, plotOutput=None, pickleOutput=None, maxLimit=None):
    import matplotlib.pyplot as plt

    data = pd.read_table(bedInput, sep="\t", header=None, names=["chr", "start", "end"])
    len_info = np.asarray(data["end"] - data["start"])

    len_info = len_info[np.where(len_info <= maxLimit)]

    a = collections.Counter(len_info)

    with open(pickleOutput, "wb") as f:
        pickle.dump(a, f)

    a = dict(sorted(a.items()))
    keys = np.fromiter(a.keys(), dtype=int)
    vals = np.fromiter(a.values(), dtype=int)

    fig = plt.figure(figsize=(10, 8))
    plt.plot(keys, vals, c="r", linewidth=1)
    plt.tick_params(labelsize=15)
    font = {
        "family": "Times New Roman",
        "weight": "normal",
        "size": 20,
    }
    plt.xlabel("DNA Fragment Size (base pair)", font)
    plt.ylabel("DNA Fragment Counts", font)
    plt.savefig(plotOutput)
    # plt.savefig(os.path.splitext(plotOutput)[0] + ".pdf")
    plt.close(fig)

    return True


def fraglenmultiplot(pickles, plotOutput, txtOutput, ratio):
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(10, 8))
    shortr, longr, peak = [], [], []
    for i in range(len(pickles)):
        with open(pickles[i], "rb") as f:
            dataInput = pickle.load(f)
            dataInput = dict(sorted(dataInput.items()))
            keys = np.fromiter(dataInput.keys(), dtype=int)
            vals = np.fromiter(dataInput.values(), dtype=int)
            vals = vals / np.sum(vals)
            shortr.append(np.sum(vals[np.where(keys < ratio[0])]))
            longr.append(np.sum(vals[np.where(keys > ratio[1])]))
            peak.append(keys[np.where(vals == np.max(vals))[0]][0])
        plt.plot(keys, vals, c="b", linewidth=0.5)
    plt.tick_params(labelsize=15)
    font = {
        "family": "Times New Roman",
        "weight": "normal",
        "size": 20,
    }
    plt.xlabel("DNA Fragment Size (base pair)", font)
    plt.ylabel("Density", font)
    plt.savefig(plotOutput)
    # plt.savefig(os.path.splitext(plotOutput)[0] + ".pdf")
    plt.close(fig)
    stat = [shortr, longr, peak, [2 * j for j in peak]]
    stat_df = pd.DataFrame(
        list(map(list, zip(*stat))),
        index=[os.path.split(k)[1] for k in pickles],
        columns=["Short(<" + str(ratio[0]) + " bp) Rate", "Long(>" + str(ratio[1]) + " bp) Rate", "Peak 1", "Peak 2"],
    )
    stat_df.to_csv(txtOutput, sep="\t", header=True, index=True)
    return True


def fraglencompplot(caseInput, ctrlInput, plotOutput, txtOutput, ratio, labelInput=["case", "control"]):
    import matplotlib.pyplot as plt

    caseprop = []
    ctrlprop = []
    fig = plt.figure(figsize=(10, 8))
    shortr, longr, peak = [], [], []
    for i in range(len(caseInput)):
        with open(caseInput[i], "rb") as f:
            dataInput = pickle.load(f)
            dataInput = dict(sorted(dataInput.items()))
            keys = np.fromiter(dataInput.keys(), dtype=int)
            vals = np.fromiter(dataInput.values(), dtype=int)
            vals = vals / np.sum(vals)
            shortr.append(np.sum(vals[np.where(keys < ratio[0])]))
            longr.append(np.sum(vals[np.where(keys > ratio[1])]))
            peak.append(keys[np.where(vals == np.max(vals))[0]][0])
            caseprop.append(np.sum(vals[:150]))
        (p1,) = plt.plot(keys, vals, c="r", linewidth=0.5)
    for i in range(len(ctrlInput)):
        with open(ctrlInput[i], "rb") as f:
            dataInput = pickle.load(f)
            dataInput = dict(sorted(dataInput.items()))
            keys = np.fromiter(dataInput.keys(), dtype=int)
            vals = np.fromiter(dataInput.values(), dtype=int)
            vals = vals / np.sum(vals)
            shortr.append(np.sum(vals[np.where(keys < ratio[0])]))
            longr.append(np.sum(vals[np.where(keys > ratio[1])]))
            peak.append(keys[np.where(vals == np.max(vals))[0]][0])
            ctrlprop.append(np.sum(vals[:150]))
        (p2,) = plt.plot(keys, vals, c="b", linewidth=0.5)
    plt.tick_params(labelsize=15)
    font = {
        "family": "Times New Roman",
        "weight": "normal",
        "size": 20,
    }
    plt.xlabel("DNA Fragment Size (base pair)", font)
    plt.ylabel("Density", font)

    font_legend = {
        "family": "Times New Roman",
        "weight": "normal",
        "size": 15,
    }
    plt.legend([p1, p2], labelInput, loc="best", prop=font_legend)
    plt.savefig(plotOutput[0])
    # plt.savefig(os.path.splitext(plotOutput[0])[0] + ".pdf")
    plt.close(fig)

    stat = [shortr, longr, peak, [2 * j for j in peak]]
    stat_df = pd.DataFrame(
        list(map(list, zip(*stat))),
        index=[os.path.split(k)[1] for k in caseInput + ctrlInput],
        columns=["Short(<" + str(ratio[0]) + " bp) Rate", "Long(>" + str(ratio[1]) + " bp) Rate", "Peak 1", "Peak 2"],
    )
    stat_df.to_csv(txtOutput, sep="\t", header=True, index=True)

    fig = plt.figure(figsize=(10, 8))
    casegory = [labelInput[0] for i in range(len(caseprop))]
    ctrlgory = [labelInput[1] for i in range(len(ctrlprop))]
    propdf = pd.DataFrame({"category": casegory + ctrlgory, "proportion": caseprop + ctrlprop})
    bp = sns.violinplot(x="category", y="proportion", data=propdf)
    bp.set_xlabel("", fontsize=20)
    bp.set_ylabel("Proportion of fragments below 150bp", fontsize=20)
    bp.tick_params(labelsize=15)
    y, h = propdf["proportion"].max() + 0.1, 0.02
    t, p = stats.ttest_ind(caseprop, ctrlprop, equal_var=False)
    if p >= 0.05:
        text = "p = " + str(p)
    elif p >= 0.01:
        text = "p < 0.05"
    elif p >= 0.001:
        text = "p < 0.01"
    elif p >= 0.0001:
        text = "p < 0.001"
    elif p >= 0:
        text = "p < 0.0001"
    plt.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1, c="k")
    plt.text(0.5, y + h, text, ha="center", va="bottom", color="k", fontdict=font_legend)
    plt.savefig(plotOutput[1])
    # plt.savefig(os.path.splitext(plotOutput[1])[0] + ".pdf")
    plt.close(fig)

    return True


# calculate methylation level for regions
def calcMethyl(bamInput, bedInput, txtOutput):
    bai = bamInput + ".bai"
    if not os.path.exists(bai):
        message = "Index file " + bai + " do not exist!"
        raise commonError(message)

    bam_input = pysam.Samfile(bamInput, "rb")
    regions = pd.read_csv(bedInput, sep="\t", header=None, names=["chr", "start", "end"])

    CXXname = ["unmCpG", "mCpG", "unmCHG", "mCHG", "unmCHH", "mCHH", "unmUNC", "mUNC"]
    d = dict.fromkeys(CXXname, 0)
    regions = regions.assign(**d)

    for index, row in regions.iterrows():
        count_data = [0, 0, 0, 0, 0, 0, 0, 0]
        for read in bam_input.fetch(reference=row["chr"], start=row["start"], end=row["end"]):
            CXXinfo = read.get_tag("XM")
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

    regions.to_csv(txtOutput, sep="\t", header=True, index=False)


# uncompress gz file
def un_gz(gzfile):
    file = gzfile.replace(".gz", "")
    with gzip.open(gzfile, "r") as f_in, open(file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    if os.path.exists(gzfile):
        os.remove(gzfile)
    else:
        raise commonError("**********Decompressing Error**********")

    return True


# run a single command line
def cmdCall(cmdLine):
    proc = subprocess.Popen(
        cmdLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
    )
    while True:
        nextline = proc.stdout.readline()
        if (nextline == "") and (proc.poll() is not None):
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()

    output, error = proc.communicate()
    exitCode = proc.returncode

    if exitCode != 0:
        print(output)
        print(error)
        raise commonError("**********CMD running error**********")


# compute OCF value for paired end data
def ComputeOCF(bedgz, txtOutput, OCFOutput, regionFile):
    print("Input file:", bedgz)
    print("Output files:", txtOutput)
    print("Output OCF files:", OCFOutput)
    print("Region files:", regionFile)

    tbx = pysam.TabixFile(bedgz)

    # regions are 1-based
    regions = pd.read_csv(regionFile, sep="\t", header=None, names=["chr", "start", "end", "description"])

    regions["OCF"] = 0
    cud_output = txtOutput

    Udata = np.empty(shape=(0, 241))
    Ddata = np.empty(shape=(0, 241))

    for idx, region in regions.iterrows():
        region_Chr, region_Start, region_End = (
            region["chr"],
            region["start"],
            region["end"],
        )

        if region_Start < 1:
            message = "Start of the region must > 0!"
            raise commonError(message)

        covPOS = defaultdict(lambda: [0, 0, 0])

        # fetch read in the region
        try:
            fetched_reads = tbx.fetch(region_Chr, region_Start, region_End)
        except ValueError:
            continue
        for row in fetched_reads:
            tmp_row = row.split()
            rstart = int(tmp_row[1]) + 1  # convert to 1-based
            rend = int(tmp_row[2])  # end included

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

        left_OCF = sum([covPOS[x][2] for x in range(midpoint - 70, midpoint - 50)]) - sum(
            [covPOS[x][1] for x in range(midpoint - 70, midpoint - 50)]
        )
        right_OCF = sum([covPOS[x][1] for x in range(midpoint + 50, midpoint + 70)]) - sum(
            [covPOS[x][2] for x in range(midpoint + 50, midpoint + 70)]
        )
        regions.loc[idx, "OCF"] = left_OCF + right_OCF

    ud_data = np.concatenate([Udata, Ddata], 1)

    np.savetxt(cud_output, ud_data, fmt="%i", delimiter="\t")

    regions.to_csv(OCFOutput, sep="\t", index=False)

    print("Processing finished!")

    return True


# compute coverage, U-end(upstream end) and D-end(downstream end)
# only count -1000 to 1000 from open region center
def computeCUE(inputFile, refFile, txtOutput, cudOutput, ocfOutput, flags):
    inputFile = pybedtools.BedTool(inputFile)
    refFile = pybedtools.BedTool(refFile)
    inputFile.intersect(refFile, wo=True, sorted=True, output=txtOutput)

    peak = 60  # distance of U and D peaks from the center
    bin = 10  # half the length of windows
    ocf = []

    data = pd.read_csv(
        txtOutput,
        sep="\t",
        header=None,
        names=["read.chr", "read.start", "read.end", "peak.chr", "peak.start", "peak.end", "description", "overlap",],
    )
    data["peak.start"] = data["peak.start"] + 1
    data["read.start"] = data["read.start"] + 1

    flag_num = -1
    for flag in flags:
        flag_num += 1
        tmp_data = data.loc[
            data.description == flag,
        ]
        cov = np.zeros(shape=2000)
        uend = np.zeros(shape=2000)
        dend = np.zeros(shape=2000)
        for idx, row in tmp_data.iterrows():
            if (row["read.start"] < row["peak.start"]) and (row["peak.start"] <= row["read.end"]):
                o_s, o_e = 0, row["read.end"] - row["peak.start"]
                cov[0 : (o_e + 1)] += 1
                dend[o_e] += 1
            elif (row["peak.start"] <= row["read.start"]) and (row["read.end"] <= row["peak.end"]):
                o_s, o_e = (
                    row["read.start"] - row["peak.start"],
                    row["read.end"] - row["peak.start"],
                )
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
        df = pd.DataFrame.from_dict(
            {
                "idx": index,
                "cov": cov,
                "cov%%": cov / cov_tot * 10000,
                "uend": uend,
                "uend%%": uend / uend_tot * 10000,
                "dend": dend,
                "dend%%": dend / dend_tot * 10000,
            }
        )
        df.to_csv(cudOutput[flag_num], sep="\t", index=False)

        trueends = 0
        background = 0
        for i in index:
            if i >= -peak - bin and i <= -peak + bin:
                trueends += dend[i + 1000] / dend_tot * 10000
                background += uend[i + 1000] / uend_tot * 10000
            elif i >= peak - bin and i <= peak + bin:
                trueends += uend[i + 1000] / uend_tot * 10000
                background += dend[i + 1000] / dend_tot * 10000
        ocf.append(trueends - background)

    ocf_df = pd.DataFrame({"tissue": flags, "OCF": ocf})
    ocf_df.to_csv(ocfOutput, sep="\t", index=None)

    return ocf


def OCF_boxplot(caseocfinput, ctrlocfinput, output, x_label):
    import matplotlib.pyplot as plt

    flag = pd.read_csv(caseocfinput[0], sep="\t", header=0, index_col=None)["tissue"].tolist()
    ocf_case = [[] for x in flag]
    ocf_ctrl = [[] for x in flag]
    for i in range(len(caseocfinput)):
        tmp = pd.read_csv(caseocfinput[i], sep="\t", header=0, index_col=None)["OCF"].tolist()
        for j in range(len(flag)):
            ocf_case[j].append(tmp[j])
    for i in range(len(ctrlocfinput)):
        tmp = pd.read_csv(ctrlocfinput[i], sep="\t", header=0, index_col=None)["OCF"].tolist()
        for j in range(len(flag)):
            ocf_ctrl[j].append(tmp[j])
    plt_fig = plt.figure()
    bpl = plt.boxplot(ocf_case, positions=np.array(range(len(ocf_case))) * 2.0 - 0.4, sym="", widths=0.6, vert=True,)
    bpr = plt.boxplot(ocf_ctrl, positions=np.array(range(len(ocf_ctrl))) * 2.0 + 0.4, sym="", widths=0.6, vert=True,)
    plt.setp(bpl["boxes"], color="y")
    plt.setp(bpl["whiskers"], color="y")
    plt.setp(bpl["caps"], color="y")
    plt.setp(bpl["medians"], color="y")
    plt.setp(bpr["boxes"], color="b")
    plt.setp(bpr["whiskers"], color="b")
    plt.setp(bpr["caps"], color="b")
    plt.setp(bpr["medians"], color="b")
    plt.plot([], c="y", label=x_label[0])
    plt.plot([], c="b", label=x_label[1])
    plt.legend()
    plt.xticks(range(0, len(flag) * 2, 2), flag)
    plt.xlim(-2, len(flag) * 2)
    for k in range(7):
        plt.scatter(
            [k * 2.0 - 0.4 for j in range(len(ocf_case[k]))], ocf_case[k], s=8, c="y",
        )
        plt.scatter(
            [k * 2.0 + 0.4 for j in range(len(ocf_ctrl[k]))], ocf_ctrl[k], s=8, c="b",
        )
    plt.ylabel("OCF value")
    plt.tight_layout()
    plt.savefig(output)
    # plt.savefig(os.path.splitext(output)[0] + ".pdf")
    plt.close(plt_fig)

    return True


def generate_cudoutput(input, outputdir):
    save_flag = ["Tcell", "Liver", "Placenta", "Lung", "Breast", "Intestine", "Ovary"]
    dict = {}
    prefix = os.path.splitext(os.path.basename(input))[0]
    for flag in save_flag:
        dict[flag] = outputdir + "/" + prefix + "-" + flag + "-cud.txt"
    return dict


# compress .bismark.zero.cov file
def compressMethy(InputFile=None, OutputFile=None):
    """
    input must from bismark_methylation_extractor, .bismark.zero.cov file
    """
    pysam.tabix_compress(InputFile, OutputFile, force=False)
    pysam.tabix_index(OutputFile, preset="bed", zerobased=True)

    return True


# compress methylation level from .bismark.zero.cov.gz file
def calcMethylV2(tbxInput, bedInput, txtOutput):
    tbi = tbxInput + ".tbi"
    if not os.path.exists(tbi):
        message = "Index file " + tbi + " do not exist!"
        raise commonError(message)

    tbx_input = pysam.TabixFile(tbxInput)

    regions = pd.read_csv(bedInput, sep="\t", header=None, names=["chr", "start", "end"])

    CXXname = ["unmCpG", "mCpG"]
    d = dict.fromkeys(CXXname, 0)
    regions = regions.assign(**d)

    message = "Now, processing fetch and computing CpG level for file " + tbxInput + "."
    print(message)

    for index, row in regions.iterrows():
        count_data = [0, 0]
        try:
            for read in tbx_input.fetch(reference=row["chr"], start=row["start"], end=row["end"]):
                readinfo = read.split()
                count_data[0] += int(readinfo[5])
                count_data[1] += int(readinfo[4])

        except ValueError:
            regions.loc[index, CXXname] = count_data
            continue
        else:
            regions.loc[index, CXXname] = count_data

    regions["mlCpG"] = regions["mCpG"] / (regions["mCpG"] + regions["unmCpG"])

    regions = regions.fillna(0)

    regions.to_csv(txtOutput, sep="\t", header=True, index=False)

    message = "Processing for file " + tbxInput + " finished."
    print(message)

    return True


def ifvalidchr(
    chr,
    validlist=[
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
    ],
):
    if chr not in validlist:
        return False
    else:
        return True


# process GC correction on read count data from CNV
def correctReadCount(
    readfileInput, gcfileInput, txtOutput, plotOutput, corrkey, readtype, sampleMaxSize=50000,
):
    if readtype == 1:
        readInput = wig2df(readfileInput)
    elif readtype == 2:
        readInput = pd.read_csv(readfileInput, sep="\t", header=0, index_col=None)
    gcInput = wig2df(gcfileInput)
    read_value = readInput["value"].values.tolist()
    read_se = readInput["start-end"].values.tolist()
    gc_value = gcInput["value"].values.tolist()
    gc_se = gcInput["start-end"].values.tolist()
    read_chr = readInput["chrom"].values.tolist()
    gc_chr = gcInput["chrom"].values.tolist()

    # matching different lines of read and gc inputs
    i, gcplus, readplus, nextpos = 0, 0, 0, 0
    reads, gc, chrs, ses = [], [], [], []
    mark = read_chr[0]
    noendappend = False
    while i + gcplus < len(gc_value) and i + readplus < len(read_value):
        if (not ifvalidchr(gc_chr[i + gcplus])) or (not ifvalidchr(read_chr[i + readplus])):
            i += 1
            continue
        # case of gc and read in different chromosomes:
        if gc_chr[i + gcplus] != read_chr[i + readplus]:
            if gc_chr[i + gcplus] != mark:
                while (not ifvalidchr(gc_chr[i + gcplus])) and (i + gcplus < len(gc_chr) - 1):
                    gcplus += 1
                    if i + gcplus == len(gc_chr) - 1:
                        noendappend = True
                if not noendappend:
                    nextpos = read_chr.index(gc_chr[i + gcplus])
                    readplus = nextpos - i
            elif read_chr[i + readplus] != mark:
                while not ifvalidchr(read_chr[i + readplus]) and (i + gcplus < len(read_chr) - 1):
                    readplus += 1
                    if i + readplus == len(read_chr) - 1:
                        noendappend = True
                if not noendappend:
                    nextpos = gc_chr.index(read_chr[i + readplus])
                    gcplus = nextpos - i
            mark = read_chr[i + readplus]
        # case of gc and read in different positions (ignore when distance is under 3 bp):
        if abs(int(gc_se[i + gcplus].split("-")[0]) - int(read_se[i + readplus].split("-")[0])) >= 3:
            while (
                int(gc_se[i + gcplus].split("-")[0]) >= int(read_se[i + readplus].split("-")[0]) + 3
                and i + readplus < len(read_se) - 1
            ):
                readplus += 1
                if i + readplus == len(read_se) - 1:
                    noendappend = True
            while (
                int(gc_se[i + gcplus].split("-")[0]) <= int(read_se[i + readplus].split("-")[0]) - 3
                and i + gcplus < len(gc_se) - 1
            ):
                gcplus += 1
                if i + gcplus == len(gc_se) - 1:
                    noendappend = True
        # second test: case of gc and read in different chromosomes:
        if gc_chr[i + gcplus] != read_chr[i + readplus]:
            if gc_chr[i + gcplus] != mark:
                while (not ifvalidchr(gc_chr[i + gcplus])) and (i + gcplus < len(gc_chr) - 1):
                    gcplus += 1
                    if i + gcplus == len(gc_chr) - 1:
                        noendappend = True
                if not noendappend:
                    nextpos = read_chr.index(gc_chr[i + gcplus])
                    readplus = nextpos - i
            elif read_chr[i + readplus] != mark:
                while not ifvalidchr(read_chr[i + readplus]) and (i + gcplus < len(read_chr) - 1):
                    readplus += 1
                    if i + readplus == len(read_chr) - 1:
                        noendappend = True
                if not noendappend:
                    nextpos = gc_chr.index(read_chr[i + readplus])
                    gcplus = nextpos - i
        if not noendappend:
            reads.append(read_value[i + readplus])
            gc.append(gc_value[i + gcplus])
            chrs.append(read_chr[i + readplus])
            ses.append(read_se[i + readplus])
            mark = read_chr[i + readplus]
        i += 1

    # run the GC_correct function
    correct_reads, correct_reads2, valid = GC_correct(reads, gc, plotOutput, corrkey, sampleMaxSize)

    readOutput = pd.DataFrame(
        {
            "chrom": [chrs[i] for i in range(len(chrs)) if valid[i]],
            "start-end": [ses[i] for i in range(len(ses)) if valid[i]],
            "value": correct_reads,
        }
    )
    readOutput2 = pd.DataFrame({"chrom": chrs, "start-end": ses, "value": correct_reads2})
    readOutput.to_csv(txtOutput, sep="\t", header=True, index=True)

    return readOutput, readOutput2


def GC_correct(readInput, gcInput, plotOutput, corrkey, sampleMaxSize):
    import matplotlib.pyplot as plt

    readl = len(readInput)
    tl = readl
    valid = [True for i in range(readl)]
    for i in range(readl):
        if readInput[i] <= 0 or gcInput[i] < 0:
            valid[i] = False
    ideal = [True for i in range(readl)]
    routlier, doutlier = 0.01, 0.001
    lrange = np.percentile(np.array([readInput[i] for i in range(readl) if valid[i]]), 0)
    rrange = np.percentile(np.array([readInput[i] for i in range(readl) if valid[i]]), (1 - routlier) * 100)
    ldomain = np.percentile(np.array([gcInput[i] for i in range(readl) if valid[i]]), doutlier * 100)
    rdomain = np.percentile(np.array([gcInput[i] for i in range(readl) if valid[i]]), (1 - doutlier) * 100)
    for i in range(readl):
        if (not valid[i]) or (not lrange <= readInput[i] <= rrange) or (not ldomain <= gcInput[i] <= rdomain):
            ideal[i] = False
            tl -= 1
    ideal_prev = np.array([[readInput[i], gcInput[i]] for i in range(readl) if ideal[i]])
    row_rand_array = np.arange(ideal_prev.shape[0])
    np.random.shuffle(row_rand_array)
    ideal_aft = ideal_prev[row_rand_array[: min(tl, sampleMaxSize)]]
    ideal_reads = [ideal_aft[j][0] for j in range(len(ideal_aft))]
    ideal_gc = [ideal_aft[j][1] for j in range(len(ideal_aft))]
    fig, (ax1, ax2) = plt.subplots(figsize=(15, 6), ncols=2)
    ax1.scatter(ideal_gc, ideal_reads, c="deepskyblue", s=0.1)
    ax1.set_xlabel("GC content")
    ax1.set_ylabel("Read Count")
    ax1.set_ylim(bottom=0)
    if corrkey == "0":
        correct_reads = [readInput[i] for i in range(readl) if valid[i]]
        correct_reads2 = readInput
    else:
        med = np.median(ideal_reads)
        rough = sm.nonparametric.lowess(endog=ideal_reads, exog=ideal_gc, frac=0.03, return_sorted=True)
        final = sm.nonparametric.lowess(endog=rough[:, 1], exog=rough[:, 0], frac=0.3, return_sorted=True)
        final_x = list(zip(*final))[0]
        final_y = list(zip(*final))[1]
        f = interp1d(final_x, final_y, bounds_error=False, fill_value="extrapolate")
        if corrkey == "/":
            correct_reads = [(readInput[i] / f(gcInput[i])) for i in range(readl) if valid[i]]
            correct_reads2 = []
            for i in range(readl):
                if valid[i]:
                    correct_reads2.append(readInput[i] / f(gcInput[i]))
                else:
                    correct_reads2.append(None)
            ax2.scatter(
                ideal_gc,
                [(ideal_reads[i] / f(ideal_gc[i])) for i in range(len(ideal_reads))],
                c="mediumaquamarine",
                s=0.1,
            )
        elif corrkey == "-":
            correct_reads, correct_reads2 = [], []
            for i in range(readl):
                if valid[i]:
                    correct_reads.append(readInput[i] - f(gcInput[i]) + med)
                    correct_reads2.append(readInput[i] - f(gcInput[i]) + med)
                else:
                    correct_reads2.append(None)
            ax2.scatter(
                ideal_gc,
                [(ideal_reads[i] - f(ideal_gc[i]) + med) for i in range(len(ideal_reads))],
                c="mediumaquamarine",
                s=0.1,
            )
        ax2.set_xlabel("GC content")
        ax2.set_ylabel("Read Count (corrected)")
        ax2.set_ylim(bottom=0)

    fig.savefig(plotOutput)
    # could not load glyph error
    # fig.savefig(os.path.splitext(plotOutput)[0] + ".pdf")
    plt.close(fig)

    return correct_reads, correct_reads2, valid


# sum the read count data to each chromosome arm
def sumChromarm(txtInput, cytoBandInput):
    dfInput = pd.read_csv(txtInput, sep="\t", header=0, index_col=0,)
    cytoBand = pd.read_csv(cytoBandInput, sep="\t", header=None, names=["chrom", "start", "end", "band", "color"],)
    sumvalue = [0]
    geneflag = [cytoBand["chrom"][0][3:] + cytoBand["band"][0][0]]
    line = 0
    for i in range(len(dfInput)):
        chrom = dfInput["chrom"][i]
        start = int(dfInput["start-end"][i].split("-")[0])
        value = float(dfInput["value"][i])
        ifget = False
        for j in range(line, len(cytoBand)):
            if chrom == cytoBand["chrom"][j] and int(cytoBand["start"][j]) <= start <= int(cytoBand["end"][j]):
                if geneflag[-1] != (chrom[3:] + cytoBand["band"][j][0]):
                    geneflag.append(chrom[3:] + cytoBand["band"][j][0])
                    sumvalue.append(0)
                sumvalue[-1] += value
                line = j
                ifget = True
                break
        if not ifget:
            for j in range(len(cytoBand)):
                if chrom == cytoBand["chrom"][j] and int(cytoBand["start"][j]) <= start <= int(cytoBand["end"][j]):
                    if geneflag[-1] != (chrom[3:] + cytoBand["band"][j][0]):
                        geneflag.append(chrom[3:] + cytoBand["band"][j][0])
                        sumvalue.append(0)
                    sumvalue[-1] += value
                    line = j
                    ifget = True
                    break

    return sumvalue, geneflag


# computer z-score and plot scatter-plot for each sample
def plotCNVscatter(
    caseInput, ctrlInput, casepos, ctrlpos, cytoBandInput, caseplotOutput, ctrlplotOutput,
):
    import matplotlib.pyplot as plt

    cytoBand = pd.read_csv(cytoBandInput, sep="\t", header=None, names=["chrom", "start", "end", "band", "color"],)
    geneflag = [cytoBand["chrom"][i][3:] + cytoBand["band"][i][0] for i in range(len(cytoBand["chrom"]))]
    intv = 100
    mean = ctrlInput.mean(axis=1)
    std = ctrlInput.std(axis=1)
    case_z = caseInput.apply(lambda x: (x - mean) / std)
    ctrl_z = ctrlInput.apply(lambda x: (x - mean) / std)
    for j in range(case_z.shape[1]):
        posnow = -intv
        nextstart = 0
        mark = None
        xpos = [[]]
        xposx = -1
        intvpos = []
        xlabels = []
        for i in range(len(casepos[j][0])):
            if casepos[j][0][i] != mark:
                mark = casepos[j][0][i]
                cytopos = cytoBand["chrom"].tolist().index(casepos[j][0][i])
                xlabels.append(geneflag[cytopos])
                cytonextpos = [
                    cytoBand["band"][k][0]
                    for k in range(len(cytoBand["band"]))
                    if cytoBand["chrom"][k] == casepos[j][0][i]
                ].index("q") + cytopos
                nextstart = int(cytoBand["start"][cytonextpos])
                intvpos.append(posnow - 1)
                posnow += intv
                intvpos.append(posnow)
                xpos.append([])
                xposx += 1
            elif casepos[j][1][i] >= nextstart and nextstart > 0:
                nextstart = -1
                xlabels.append(geneflag[cytonextpos])
                intvpos.append(posnow - 1)
                posnow += intv
                intvpos.append(posnow)
                xpos.append([])
                xposx += 1
            xpos[xposx].append(posnow)
            posnow += 1
        intvpos.append(posnow - 1)
        intvmidpos = [(intvpos[2 * k + 1] + intvpos[2 * k + 2]) / 2 for k in range(int(len(intvpos) / 2))]
        f, (ax) = plt.subplots(figsize=(50, 10))
        for m in range(len(xpos)):
            ax.scatter(
                xpos[m],
                case_z.iloc[:, j].tolist()[
                    sum([len(xpos[k]) for k in range(m)]) : sum([len(xpos[k]) for k in range(m)]) + len(xpos[m])
                ],
                c="black",
                s=0.15,
            )
            # draw the chromosome marking line
            ax.plot(xpos[m], [-5 for k in range(len(xpos[m]))], color="black", linewidth=2)
        ax.set_ylabel("Z-score")
        ax.set_xticks(intvmidpos)
        ax.set_xticklabels(xlabels)
        ax.set_xlim([0, max(max(xpos)) + 1])
        ax.set_ylim([-5, 5])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        plt.savefig(caseplotOutput[j])
        # plt.savefig(os.path.splitext(caseplotOutput[j])[0] + ".pdf")
        plt.close(f)
    for j in range(ctrl_z.shape[1]):
        posnow = -intv
        nextstart = 0
        mark = None
        xpos = [[]]
        xposx = -1
        intvpos = []
        xlabels = []
        for i in range(len(ctrlpos[j][0])):
            if ctrlpos[j][0][i] != mark:
                mark = ctrlpos[j][0][i]
                cytopos = cytoBand["chrom"].tolist().index(ctrlpos[j][0][i])
                xlabels.append(geneflag[cytopos])
                cytonextpos = [
                    cytoBand["band"][k][0]
                    for k in range(len(cytoBand["band"]))
                    if cytoBand["chrom"][k] == ctrlpos[j][0][i]
                ].index("q") + cytopos
                nextstart = int(cytoBand["start"][cytonextpos])
                intvpos.append(posnow - 1)
                posnow += intv
                intvpos.append(posnow)
                xpos.append([])
                xposx += 1
            elif ctrlpos[j][1][i] >= nextstart and nextstart > 0:
                nextstart = -1
                xlabels.append(geneflag[cytonextpos])
                intvpos.append(posnow - 1)
                posnow += intv
                intvpos.append(posnow)
                xpos.append([])
                xposx += 1
            xpos[xposx].append(posnow)
            posnow += 1
        intvpos.append(posnow - 1)
        intvmidpos = [(intvpos[2 * k + 1] + intvpos[2 * k + 2]) / 2 for k in range(int(len(intvpos) / 2))]
        f, (ax) = plt.subplots(figsize=(50, 10))
        for m in range(len(xpos)):
            ax.scatter(
                xpos[m],
                ctrl_z.iloc[:, j].tolist()[
                    sum([len(xpos[k]) for k in range(m)]) : sum([len(xpos[k]) for k in range(m)]) + len(xpos[m])
                ],
                c="black",
                s=0.15,
            )
            # draw the chromosome marking line
            ax.plot(xpos[m], [-5 for k in range(len(xpos[m]))], color="black", linewidth=2)
        ax.set_ylabel("Z-score")
        ax.set_xticks(intvmidpos)
        ax.set_xticklabels(xlabels)
        ax.set_xlim([0, max(max(xpos)) + 1])
        ax.set_ylim([-5, 5])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        plt.savefig(ctrlplotOutput[j])
        # plt.savefig(os.path.splitext(ctrlplotOutput[j])[0] + ".pdf")
        plt.close(f)

    return True


# compute z-score and plot the heatmap for whole samples
def plotCNVheatmap(caseInput, ctrlInput, txtOutput, plotOutput):
    import matplotlib.pyplot as plt

    mean = ctrlInput.mean(axis=1)
    std = ctrlInput.std(axis=1)
    case_z = caseInput.apply(lambda x: (x - mean) / std)
    ctrl_z = ctrlInput.apply(lambda x: (x - mean) / std)
    data = pd.concat([ctrl_z, case_z], axis=1)

    colormap = [
        (0.8622837370242215, 0.42952710495963087, 0.34271434063821604),
        (0.9686274509803922, 0.7176470588235293, 0.5999999999999999),
        (0.982006920415225, 0.9061899269511726, 0.8615916955017301),
        (0.9657054978854287, 0.9672433679354094, 0.9680891964628989),
        (0.9657054978854287, 0.9672433679354094, 0.9680891964628989),
        (0.8838908112264514, 0.9284890426758939, 0.9530180699730872),
        (0.654901960784314, 0.8143790849673205, 0.8941176470588236),
        (0.3234909650134564, 0.6149173394848135, 0.7854671280276817),
    ]

    f, (ax) = plt.subplots(figsize=(20, 20))
    sns.heatmap(data, ax=ax, cmap=colormap, vmin=-4, vmax=4)
    data.to_csv(txtOutput, sep="\t", index=True)
    f.savefig(plotOutput)
    # plt.savefig(os.path.splitext(plotOutput)[0] + ".pdf")
    plt.close(f)

    return True


# read .wig files into dataframe
def wig2df(inputfile):
    f = open(inputfile, "r")
    data = f.readlines()
    chrom_list = []
    startend_list = []
    value_list = []
    for k in data:
        if k.split(" ")[0] == "fixedStep":
            dec = k.split(" ")
            chrom = dec[1].split("=")[1]
            start = int(dec[2].split("=")[1])
            step = int(dec[3].split("=")[1])
            if len(dec) == 5:
                span = int(dec[4].split("=")[1])
            else:
                span = 1
            rownum = -1
        else:
            rownum += 1
            chrom_list.append(chrom)
            startend_list.append(str(rownum * step + start) + "-" + str(rownum * step + start + span - 1))
            value_list.append(float(k.strip("\n")))
    df = pd.DataFrame({"chrom": chrom_list, "start-end": startend_list, "value": value_list})
    return df


def mixedConditionNumber(referenceSelect, markerSelect):
    pinvReferenceSelect = np.linalg.pinv(referenceSelect)
    bNorm = np.linalg.norm(markerSelect)
    maxConditionNumber = 0
    tissueNumber = np.size(referenceSelect, 1)
    for i in range(tissueNumber):
        tq = pinvReferenceSelect[i, :]
        conditionNumber = (bNorm / np.abs(np.dot(tq, markerSelect))) * np.linalg.norm(tq)
        if conditionNumber > maxConditionNumber:
            maxConditionNumber = conditionNumber
    return maxConditionNumber


def markerFromSort(sortIndex, tissueNum, markerNum, selectedNumber):
    selectedMarker = np.zeros(markerNum)
    selectedMarker = selectedMarker.astype(np.bool)
    finalMarker = np.arange(markerNum)
    for i in range(tissueNum):
        if selectedNumber > len(sortIndex[i]):
            selectedMarker[sortIndex[i][:]] = True
        else:
            selectedMarker[sortIndex[i][0:selectedNumber]] = True
    return finalMarker[selectedMarker]


def iterativeWeightedSVR(reference, markerData):
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(markerData, 1)
    # markerNumber = np.size(reference, 0)
    proportionDeconvolution = np.zeros([tissueNum, numOfSamples])
    for i in tqdm(range(numOfSamples)):
        iterNumber = 10
        iterReference = reference
        iterMarkerData = markerData[:, i]
        # Threshold Searching
        mixture = optimize.nnls(iterReference, iterMarkerData)
        test = mixture[0]
        t = test / np.sum(test)
        c1 = np.zeros([np.size(iterReference, 0), 1])
        c1[:, 0] = iterMarkerData[:]
        t1 = np.zeros([tissueNum, 1])
        t1[:, 0] = t
        c2 = np.dot(iterReference, t1)
        error = np.abs(c2 - c1)
        error = error[:, 0]
        error = np.sort(error)
        index = int(np.size(error) * 0.8)
        threshold = error[index]
        for j in range(iterNumber):
            mixture = optimize.nnls(iterReference, iterMarkerData)
            test = mixture[0]
            t = test / np.sum(test)
            c1 = np.zeros([np.size(iterReference, 0), 1])
            c1[:, 0] = iterMarkerData[:]
            t1 = np.zeros([tissueNum, 1])
            t1[:, 0] = t
            c2 = np.dot(iterReference, t1)
            error = np.abs(c2 - c1)
            error = error[:, 0]
            iterSelected = np.arange(np.size(iterReference, 0))
            area = error < threshold
            iterSelected = iterSelected[area]
            iterReference = iterReference[iterSelected, :]
            iterMarkerData = iterMarkerData[iterSelected]
        singleMarker = np.zeros([np.size(iterReference, 0), 1])
        singleMarker[:, 0] = iterMarkerData
        t = nuSVR(iterReference, singleMarker)
        t = t[:, 0]
        proportionDeconvolution[:, i] = t
    return proportionDeconvolution


def predecipher(mixInput, refInput):
    multi_run_len = len(mixInput)
    mix = [[] for i in range(multi_run_len)]
    ref_pd = pd.read_csv(refInput, sep="\t", header=0)
    ref_depos = ref_pd.iloc[:, 1:]
    celltypes = ref_depos.columns.values.tolist()
    ref = [[] for i in range(ref_depos.shape[0])]
    for i in range(ref_depos.shape[0]):
        ref[i] = ref_depos.iloc[i].tolist()
    for i in range(multi_run_len):
        data = pd.read_csv(mixInput[i], sep="\t", header=0, names=["chr", "start", "end", "unmCpG", "mCpG", "mlCpG"],)
        mix[i] = data["mlCpG"].tolist()
    mix = pd.DataFrame(np.array(np.transpose(mix)), columns=[os.path.split(x)[1] for x in mixInput])
    ref = pd.DataFrame(np.array(ref), columns=celltypes)
    return mix, ref


def decipher(
    ref,
    mix,
    save_path="prop_predict.csv",
    marker_path="",
    scale=0.1,
    delcol_factor=10,
    iter_num=10,
    confidence=0.75,
    w_thresh=10,
    unknown=False,
    is_markers=False,
    is_methylation=False,
):
    print("---------------------------------------------")
    print("---------------Deconvolotion-----------------")
    print("---------------------------------------------")
    cell_type = ref.columns.values
    samples = mix.columns.values
    prop = collections.OrderedDict()
    prop["cell types"] = cell_type
    if is_markers:
        reference = []
        mixture = []
        markers = pd.read_csv(marker_path, index_col=0)
        markers = markers.index.values
        for i in range(len(markers)):
            reference.append(ref.loc[markers[i]])
            mixture.append(mix.loc[markers[i]])
    else:
        reference, mixture = ref, mix
    reference = np.asarray(reference)
    mixture = np.asarray(mixture)
    if is_methylation:
        reference, mixture = pre_marker_select(reference, mixture)
    print("Data reading finished!")
    print("RareDecipher Engines Start, Please Wait......")
    prop_predict = DECONVO(
        scale * reference,
        scale * mixture,
        delcol_factor=delcol_factor,
        iter_num=iter_num,
        confidence=confidence,
        w_thresh=w_thresh,
        unknown=unknown,
    )
    print("Deconvo Results Saving!")
    for i in range(len(samples)):
        prop[samples[i]] = []
        for j in range(len(cell_type)):
            prop[samples[i]].append(prop_predict[j, i])
    prop = pd.DataFrame(prop)
    prop.to_csv(save_path, sep="\t", index=False)
    print("Finished!")
    return prop


def DECONVO(ref, mix, delcol_factor=10, iter_num=10, confidence=0.75, w_thresh=10, unknown=False):
    reference, mixtureData = filt_zeros(ref, mix)
    markerNum = np.size(reference, 0)
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(mixtureData, 1)
    bestReference = []
    bestMarker = []
    conditionNumberHistory = []
    bestNumberHistory = []
    proportionDeconvolution = (
        np.zeros([tissueNum, numOfSamples]) if unknown == False else np.zeros([tissueNum + 1, numOfSamples])
    )
    for i in tqdm(range(numOfSamples)):
        selectArea = np.arange(markerNum)
        selectMixture = mixtureData[selectArea, i]
        selectReference = reference[selectArea, :]
        minimumConditionNumber = 10 ** (100)
        endNumber = np.size(selectReference, 0)
        for selectedNumber in range(int(endNumber / delcol_factor)):
            minDistances = 10 ** (50)
            for j in range(tissueNum):
                for k in range(j + 1, tissueNum):
                    distances = selectReference[:, j] - selectReference[:, k]
                    distances = np.sqrt(np.sum(np.multiply(distances, distances)))
                    if distances < minDistances:
                        minDistances = distances
                        closetJ = j
                        closetK = k
            sumData = selectReference[:, closetJ] + selectReference[:, closetK]
            area = sumData == 0
            sumData[area] = 10 ** (-100)
            collinearity = np.abs(selectReference[:, closetJ] - selectReference[:, closetK]) / (sumData)
            collinearityIndex = np.argsort(collinearity)
            area = np.ones(np.size(selectReference, 0))
            area = area.astype(np.bool)
            area[collinearityIndex[0]] = False
            selectArea = np.arange(np.size(selectReference, 0))
            selectArea = selectArea[area]
            selectMixture = selectMixture[selectArea]
            selectReference = selectReference[selectArea, :]
            ConditionNumber = cmpConditionNumber(selectReference, selectMixture)
            conditionNumberHistory.append(ConditionNumber)
            if ConditionNumber < minimumConditionNumber:
                minimumConditionNumber = ConditionNumber
                bestReference = selectReference
                bestMarker = np.zeros([np.size(selectReference, 0), 1])
                bestMarker[:, 0] = selectMixture
                bestNumber = selectedNumber
        t = RobustSVR(
            bestReference, bestMarker, iter_num=iter_num, confidence=confidence, w_thresh=w_thresh, unknown=unknown
        )
        bestNumberHistory.append(bestNumber)
        proportionDeconvolution[:, i] = t[:, 0]
    return proportionDeconvolution


def cmpConditionNumber(referenceSelect, mixtureSelect):
    pinvReferenceSelect = np.linalg.pinv(referenceSelect)
    bNorm = np.linalg.norm(mixtureSelect)
    maxConditionNumber = 0
    tissueNumber = np.size(referenceSelect, 1)
    for i in range(tissueNumber):
        tq = pinvReferenceSelect[i, :]
        conditionNumber = (bNorm / np.abs(np.dot(tq, mixtureSelect))) * np.linalg.norm(tq)
        if conditionNumber > maxConditionNumber:
            maxConditionNumber = conditionNumber
    return maxConditionNumber


def stdRes(res, d, H):
    res_std = np.zeros([np.size(res, 0), np.size(res, 1)])
    s = np.sqrt(np.sum(np.power(res, 2)) / d)
    for i in range(np.size(res, 0)):
        res_std[i, 0] = res[i, 0] / (s * (1 - H[i, i]))
    return res_std


def RobustSVR(reference, mixtureData, iter_num=10, confidence=0.75, w_thresh=10, unknown=False):
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(mixtureData, 1)
    markerNumber = np.size(reference, 0)
    proportionDeconvolution = (
        np.zeros([tissueNum, numOfSamples]) if unknown == False else np.zeros([tissueNum + 1, numOfSamples])
    )
    for i in range(numOfSamples):
        iterReference = reference
        itermixtureData = mixtureData[:, i]
        mixture = sm.RLM(itermixtureData, iterReference).fit()
        test = mixture.params
        t = test / np.sum(test) if unknown == False else test
        c1 = np.zeros([np.size(iterReference, 0), 1])
        c1[:, 0] = itermixtureData[:]
        t1 = np.zeros([tissueNum, 1])
        t1[:, 0] = t
        c2 = np.dot(iterReference, t1)
        res = c1 - c2
        s_2 = np.sum(np.power(res, 2)) / (np.size(iterReference, 0) - np.size(iterReference, 1))
        res_std = np.abs(res / np.sqrt(s_2))
        res_Sort = np.sort(res_std[:, 0])
        T = res_Sort[int(confidence * np.size(res_Sort))]
        memRef = np.zeros([np.size(iterReference, 0), np.size(iterReference, 1)])
        memRef[:, :] = iterReference[:, :]
        memMix = np.zeros(np.size(itermixtureData))
        memMix[:] = itermixtureData[:]
        for j in range(iter_num):
            mixture = sm.RLM(itermixtureData, iterReference).fit()
            test = mixture.params
            t = test / np.sum(test) if unknown == False else test
            c1 = np.zeros([np.size(iterReference, 0), 1])
            c1[:, 0] = itermixtureData[:]
            t1 = np.zeros([tissueNum, 1])
            t1[:, 0] = t
            c2 = np.dot(iterReference, t1)
            res = c1 - c2
            s_2 = np.sum(np.power(res, 2)) / (np.size(iterReference, 0) - np.size(iterReference, 1))
            res_std = res / np.sqrt(s_2)
            iterSelected = np.arange(np.size(iterReference, 0))
            area = np.abs(res_std[:, 0]) <= T
            iterSelected = iterSelected[area]
            iterReference = iterReference[iterSelected, :]
            itermixtureData = itermixtureData[iterSelected]
            if np.size(iterReference, 0) < int(tissueNum):
                iterReference = memRef
                itermixtureData = memMix
                break
            if np.size(iterReference, 0) < int(0.5 * markerNumber):
                break
            memRef = np.zeros([np.size(iterReference, 0), np.size(iterReference, 1)])
            memRef[:, :] = iterReference[:, :]
            memMix = np.zeros(np.size(itermixtureData))
            memMix[:] = itermixtureData[:]
        weights = weightsDesigner(iterReference, itermixtureData, w_thresh=w_thresh)
        t = nuSVR(iterReference, itermixtureData.reshape([-1, 1]), weights, unknown=unknown)
        t = t[:, 0]
        if unknown == False:
            proportionDeconvolution[:, i] = t
        else:
            proportionDeconvolution[0:-1, i] = t
            proportionDeconvolution[-1, i] = max(1 - np.sum(t), 0)
    return proportionDeconvolution


def weightsDesigner(ref, mix, w_thresh=10):
    mixture = sm.RLM(mix, ref).fit()
    test = mixture.params
    x_pre = test / np.sum(test)
    weights = np.abs(np.dot(ref, x_pre))
    for i in range(np.size(weights)):
        weights[i] = 1 / weights[i]
    weights = weights / np.min(weights)
    for i in range(np.size(weights)):
        if weights[i] > w_thresh:
            weights[i] = w_thresh
    return weights / np.mean(weights)


def nuSVR(reference, mixtureData, weights, unknown=False):
    nu = [0.25, 0.50, 0.75]
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(mixtureData, 1)
    proportionDeconvolution = np.zeros([tissueNum, numOfSamples])
    nuHistory = []
    p0 = np.zeros([3, tissueNum, numOfSamples])
    for i in range(0, 3, 1):
        nuI = nu[i]
        clf = NuSVR(nu=nuI, kernel="linear")
        for j in range(numOfSamples):
            clf.fit(reference, mixtureData[:, j], sample_weight=weights)
            t = clf.coef_
            t1 = np.zeros(tissueNum)
            t1[:] = t[0, :]
            area = t1 < 0
            t1[area] = 0
            t1 = t1 / np.sum(t1) if unknown == False else t1
            p0[i, :, j] = t1[:]
    for i in range(numOfSamples):
        minRMSE = 10 ** (50)
        truth = np.zeros([np.size(reference, 0), 1])
        truth[:, 0] = mixtureData[:, i]
        for k in range(0, 3, 1):
            pVector = np.zeros([tissueNum, 1])
            pVector[:, 0] = p0[k, :, i]
            temp = np.dot(reference, pVector)
            error = truth - temp
            error = np.sqrt(np.mean(np.multiply(error, error)))
            if error < minRMSE:
                minRMSE = error
                proportionDeconvolution[:, i] = p0[k, :, i]
                bestNu = k
        nuHistory.append(nu[bestNu])
    return proportionDeconvolution


def pre_marker_select(reference, mixtureData):
    cellTypeNumber = np.size(reference, 1)
    markerNumber = np.size(reference, 0)
    selectedCpG = np.arange(markerNumber)
    area = np.zeros(markerNumber)
    area = area.astype(np.bool)
    for i in range(cellTypeNumber):
        for j in range(i + 1, cellTypeNumber):
            temp = reference[:, [i, j]]
            tempSum = np.sum(temp, axis=1)
            tempSum1 = np.zeros([markerNumber, 2])
            tempSum1[:, 0] = tempSum
            tempSum1[:, 1] = tempSum
            temp = temp / tempSum1
            pairSortIndexIncrease = np.argsort(temp[:, 0])
            pairSortIndexDecrease = np.argsort(temp[:, 1])
            area[pairSortIndexIncrease[0:100]] = True
            area[pairSortIndexDecrease[0:100]] = True
    selectedCpG = selectedCpG[area]
    ref = reference[selectedCpG, :]
    mix = mixtureData[selectedCpG, :]
    return ref, mix


def filt_zeros(ref, mix):
    ref1 = []
    mix1 = []
    for i in range(np.size(ref, 0)):
        if np.max(ref[i, :]) > 0:
            ref1.append(ref[i, :])
            mix1.append(mix[i, :])
    return np.asarray(ref1), np.asarray(mix1)


def deconplot(mixInput, plotOutput, maxSample=5):
    import matplotlib.pyplot as plt

    celltypes = mixInput.iloc[:, 0].tolist()
    mixInput = mixInput.iloc[:, 1:]
    if maxSample > 0 and maxSample < mixInput.shape[1]:
        mixInput = mixInput.iloc[:, :maxSample]
    r = np.arange(mixInput.shape[1])
    bot = [0 for i in range(mixInput.shape[1])]

    for i in range(mixInput.shape[0]):
        color = "#"
        for k in range(6):
            random.seed((i + 1) * (k + 1) / (i * k + 6) + 209)
            tmp = random.choice("0123456789ABCDEF")
            color += tmp
        plt.bar(
            r, mixInput.iloc[i, :].tolist(), bottom=bot, color=color, edgecolor="white", width=0.6, label=celltypes[i],
        )
        for j in range(mixInput.shape[1]):
            bot[j] += mixInput.iloc[i, j]

    plt.xticks(r, mixInput.columns.values.tolist())
    plt.ylabel("proportion")
    plt.legend(
        bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0,
    )
    plt.savefig(plotOutput, bbox_inches="tight")
    # plt.savefig(os.path.splitext(plotOutput)[0] + ".pdf")
    plt.close()

    return True


def processPCA(casetxtInput, ctrltxtInput):
    case_multi_run_len = len(casetxtInput)
    ctrl_multi_run_len = len(ctrltxtInput)
    ml = [[] for i in range(case_multi_run_len + ctrl_multi_run_len)]
    for i in range(case_multi_run_len):
        data = pd.read_csv(
            casetxtInput[i], sep="\t", header=0, names=["chr", "start", "end", "unmCpG", "mCpG", "mlCpG"],
        )
        ml[i] = data["mlCpG"].tolist()
    for i in range(ctrl_multi_run_len):
        data = pd.read_csv(
            ctrltxtInput[i], sep="\t", header=0, names=["chr", "start", "end", "unmCpG", "mCpG", "mlCpG"],
        )
        ml[i + case_multi_run_len] = data["mlCpG"].tolist()
    pca = PCA(n_components=2)
    pca.fit(ml)
    data_2d = pca.fit_transform(ml)
    return data_2d[:case_multi_run_len], data_2d[case_multi_run_len:]


def clusterplot(casedata, ctrldata, plotOutput, labels=["case", "control"]):
    import matplotlib.pyplot as plt

    theta = np.concatenate((np.linspace(-np.pi, np.pi, 50), np.linspace(np.pi, -np.pi, 50)))
    circle = np.array((np.cos(theta), np.sin(theta)))
    casesigma = np.cov(np.array((casedata[:, 0], casedata[:, 1])))
    ctrlsigma = np.cov(np.array((ctrldata[:, 0], ctrldata[:, 1])))
    ed = np.sqrt(stats.chi2.ppf(0.95, 2))
    ell1 = np.transpose(circle).dot(np.linalg.cholesky(casesigma) * ed)
    ell2 = np.transpose(circle).dot(np.linalg.cholesky(ctrlsigma) * ed)
    a1, b1 = np.max(ell1[:, 0]), np.max(ell1[:, 1])
    a2, b2 = np.max(ell2[:, 0]), np.max(ell2[:, 1])
    t = np.linspace(0, 2 * np.pi, 100)
    fig = plt.figure(figsize=(10, 10))
    p1 = plt.scatter(casedata[:, 0], casedata[:, 1], c="y")
    plt.plot(a1 * np.cos(t), b1 * np.sin(t), c="y")
    p2 = plt.scatter(ctrldata[:, 0], ctrldata[:, 1], c="b")
    plt.plot(a2 * np.cos(t), b2 * np.sin(t), c="b")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend([p1, p2], labels, loc="best")
    plt.savefig(plotOutput)
    # plt.savefig(os.path.splitext(plotOutput)[0] + ".pdf")
    plt.close(fig)


def divide_bin_1(chromsize, blacklist, gap, windows, binlen):
    a = pybedtools.BedTool(chromsize)
    a1 = pybedtools.BedTool(blacklist)
    a2 = pybedtools.BedTool(gap)
    bins_init = a.window_maker(w=binlen, g=chromsize)
    bins_bl = bins_init.subtract(a1)
    bins_gap = bins_bl.subtract(a2)
    bins_fin = bins_gap.filter(
        lambda x: x.end - x.start >= binlen / 10
        and x.chrom
        in [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
        ]
    )
    bins_fin.saveas(windows)
    return True


def divide_bin_2(chromsize, windows, binlen):
    a = pybedtools.BedTool(chromsize)
    bins_init = a.window_maker(w=binlen, g=chromsize)
    bins_fin = bins_init.filter(
        lambda x: x.chrom
        in [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
            "chrX",
            "chrY",
        ]
    )
    bins_fin.saveas(windows)
    return True


def count_short_long(windows, bedgz, binlen, domain):
    print("Processing", bedgz, "...")
    f = pysam.Tabixfile(filename=bedgz, mode="r")
    temp_bins = pybedtools.BedTool(windows)
    length = len(temp_bins)
    bins = pybedtools.BedTool(windows)
    prev_start = 1
    prev_end = prev_start + binlen - 1
    first_mark = True
    chrom_mark = None
    shorts, longs = 0, 0
    shorts_data, longs_data, pos_data = [], [], [[], []]
    for iter in range(length):
        bin = bins[iter]
        # case of changing chromosome
        if bin.chrom != chrom_mark and chrom_mark is not None:
            prev_start = 1
            prev_end = prev_start + binlen - 1
            chrom_mark = bin.chrom
        # case of this bin is in the same 5Mb-bin of the previous bin
        if prev_start < bin.start <= prev_end:
            try:
                f.fetch(bin.chrom, bin.start, bin.end)
            except ValueError:
                pass
            else:
                bin_data = []
                for read in f.fetch(bin.chrom, bin.start, bin.end):
                    bin_data.append(int(read.split("\t")[2]) - int(read.split("\t")[1]))
                count = np.bincount(bin_data, minlength=domain[3] + 1)
                shorts += sum(count[domain[0] : domain[1] + 1])
                longs += sum(count[domain[2] : domain[3] + 1])
                first_mark = False
        # case of this bin is in a new 5Mb-bin in the same chromosome
        elif bin.start > prev_end:
            if not first_mark:
                shorts_data.append(shorts)
                longs_data.append(longs)
                pos_data[0].append(bin.chrom)
                pos_data[1].append(prev_start)
            else:
                first_mark = False
            shorts, longs = 0, 0
            try:
                f.fetch(bin.chrom, bin.start, bin.end)
            except ValueError:
                pass
            else:
                bin_data = []
                for read in f.fetch(bin.chrom, bin.start, bin.end):
                    bin_data.append(int(read.split("\t")[2]) - int(read.split("\t")[1]))
                count = np.bincount(bin_data, minlength=domain[3] + 1)
                shorts += sum(count[domain[0] : domain[1] + 1])
                longs += sum(count[domain[2] : domain[3] + 1])
            while bin.start > prev_end:
                prev_start += binlen
                prev_end += binlen
        chrom_mark = bin.chrom
    shorts_data.append(shorts)
    longs_data.append(longs)
    pos_data[0].append(bin.chrom)
    pos_data[1].append(prev_start)
    shorts_df = pd.DataFrame(
        {
            "chrom": pos_data[0],
            "start-end": [
                str(pos_data[1][i]) + "-" + str(pos_data[1][i] + binlen - 1) for i in range(len(pos_data[1]))
            ],
            "value": shorts_data,
        }
    )
    longs_df = pd.DataFrame(
        {
            "chrom": pos_data[0],
            "start-end": [
                str(pos_data[1][i]) + "-" + str(pos_data[1][i] + binlen - 1) for i in range(len(pos_data[1]))
            ],
            "value": longs_data,
        }
    )
    return shorts_df, longs_df


def count_read(windows, bedgz, binlen):
    print("Processing", bedgz, "...")
    f = pysam.Tabixfile(filename=bedgz, mode="r")
    temp_bins = pybedtools.BedTool(windows)
    length = len(temp_bins)
    bins = pybedtools.BedTool(windows)
    prev_start = 1
    prev_end = prev_start + binlen - 1
    first_mark = True
    chrom_mark = None
    reads = 0
    reads_data, pos_data = [], [[], []]
    for iter in range(length):
        bin = bins[iter]
        # case of changing chromosome
        if bin.chrom != chrom_mark and chrom_mark is not None:
            prev_start = 1
            prev_end = prev_start + binlen - 1
            chrom_mark = bin.chrom
        # case of this bin is in the same 100kb-bin of the previous bin
        if prev_start < bin.start <= prev_end:
            try:
                f.fetch(bin.chrom, bin.start, bin.end)
            except ValueError:
                pass
            else:
                bin_data = []
                for read in f.fetch(bin.chrom, bin.start, bin.end):
                    bin_data.append(int(read.split("\t")[2]) - int(read.split("\t")[1]))
                reads += len(bin_data)
                first_mark = False
        # case of this bin is in a new 5Mb-bin in the same chromosome
        elif bin.start > prev_end:
            if not first_mark:
                reads_data.append(reads)
                pos_data[0].append(bin.chrom)
                pos_data[1].append(prev_start)
            else:
                first_mark = False
            reads = 0
            try:
                f.fetch(bin.chrom, bin.start, bin.end)
            except ValueError:
                pass
            else:
                bin_data = []
                for read in f.fetch(bin.chrom, bin.start, bin.end):
                    bin_data.append(int(read.split("\t")[2]) - int(read.split("\t")[1]))
                reads += len(bin_data)
            while bin.start > prev_end:
                prev_start += binlen
                prev_end += binlen
        chrom_mark = bin.chrom
    reads_data.append(reads)
    pos_data[0].append(bin.chrom)
    pos_data[1].append(prev_start)
    reads_df = pd.DataFrame(
        {
            "chrom": pos_data[0],
            "start-end": [
                str(pos_data[1][i]) + "-" + str(pos_data[1][i] + binlen - 1) for i in range(len(pos_data[1]))
            ],
            "value": reads_data,
        }
    )
    return reads_df


def count_fragprof(
    bedgzInput=None, bedOutput=None, txtOutput=None, domain=[100, 150, 151, 220], binlen=5000000, type=None,
):
    if type == 1:
        shorts_df, longs_df = count_short_long(bedOutput, bedgzInput, binlen, domain)
        shorts_df.to_csv(txtOutput[0], sep="\t", header=True, index=None)
        longs_df.to_csv(txtOutput[1], sep="\t", header=True, index=None)
    elif type == 2:
        reads_df = count_read(bedOutput, bedgzInput, binlen)
        reads_df.to_csv(txtOutput, sep="\t", header=True, index=None)
    return True


def fragProfileplot(
    casetxtInput, ctrltxtInput, cytoBandInput, plotOutput, txtOutput, labels=["case", "control"],
):
    import matplotlib.pyplot as plt

    caseInput = []
    ctrlInput = []
    for i in range(len(casetxtInput)):
        caseInput.append(pd.read_csv(casetxtInput[i], sep="\t", header=0))
    for i in range(len(ctrltxtInput)):
        ctrlInput.append(pd.read_csv(ctrltxtInput[i], sep="\t", header=0))
    case_fp = []
    ctrl_fp = []
    casepos = [
        caseInput[0]["chrom"].tolist(),
        [
            int(caseInput[0]["start-end"].tolist()[i].split("-")[0])
            for i in range(len(caseInput[0]["start-end"].tolist()))
        ],
    ]
    for i in range(int(len(caseInput) / 2)):
        shorts_temp = caseInput[2 * i]["value"].tolist()
        longs_temp = caseInput[2 * i + 1]["value"].tolist()
        case_fp.append([shorts_temp[j] / longs_temp[j] for j in range(len(shorts_temp))])
    for i in range(int(len(ctrlInput) / 2)):
        shorts_temp = ctrlInput[2 * i]["value"].tolist()
        longs_temp = ctrlInput[2 * i + 1]["value"].tolist()
        ctrl_fp.append([shorts_temp[j] / longs_temp[j] for j in range(len(shorts_temp))])
    data_df = pd.DataFrame({"chrom": casepos[0], "start-end": caseInput[0]["start-end"].tolist()})
    for i in range(int(len(caseInput) / 2)):
        data_df = pd.concat(
            [data_df, pd.DataFrame({casetxtInput[2 * i].split("/")[-1].split("_short")[0]: case_fp[i]}),], axis=1,
        )
    for i in range(int(len(ctrlInput) / 2)):
        data_df = pd.concat(
            [data_df, pd.DataFrame({ctrltxtInput[2 * i].split("/")[-1].split("_short")[0]: ctrl_fp[i]}),], axis=1,
        )
    data_df.to_csv(txtOutput, sep="\t", header=True, index=None)

    cytoBand = pd.read_csv(cytoBandInput, sep="\t", header=None, names=["chrom", "start", "end", "band", "color"],)
    try:
        geneflag = [cytoBand["chrom"][i][3:] + cytoBand["band"][i][0] for i in range(len(cytoBand["chrom"]))]
    except TypeError:
        raise commonError("Low sequencing depth detected, cannot get valid value for specific genome region.")
    intv = 1
    mark = None
    posnow = -intv
    nextstart = 0
    xpos = [[]]
    xposx = -1
    intvpos = []
    xlabels = []
    for i in range(len(casepos[0])):
        if casepos[0][i] != mark:
            mark = casepos[0][i]
            cytopos = cytoBand["chrom"].tolist().index(casepos[0][i])
            xlabels.append(geneflag[cytopos])
            cytonextpos = [
                cytoBand["band"][k][0] for k in range(len(cytoBand["band"])) if cytoBand["chrom"][k] == casepos[0][i]
            ].index("q") + cytopos
            nextstart = int(cytoBand["start"][cytonextpos])
            intvpos.append(posnow - 1)
            posnow += intv
            intvpos.append(posnow)
            xpos.append([])
            xposx += 1
        elif casepos[1][i] >= nextstart and nextstart > 0:
            nextstart = -1
            xlabels.append(geneflag[cytonextpos])
            intvpos.append(posnow - 1)
            posnow += intv
            intvpos.append(posnow)
            xpos.append([])
            xposx += 1
        xpos[xposx].append(posnow)
        posnow += 1
    intvpos.append(posnow - 1)
    intvmidpos = [(intvpos[2 * k + 1] + intvpos[2 * k + 2]) / 2 for k in range(int(len(intvpos) / 2))]

    f, (ax1, ax2) = plt.subplots(2, 1, figsize=[50, 10])
    for i in range(len(case_fp)):
        case_fp[i] -= np.mean(case_fp[i])
        for j in range(len(xpos)):
            ax1.plot(
                xpos[j],
                case_fp[i][
                    sum([len(xpos[k]) for k in range(j)]) : sum([len(xpos[k]) for k in range(j)]) + len(xpos[j])
                ],
                color="black",
                marker=".",
                markersize=0.3,
                linewidth=0.3,
            )
            # draw the chromosome marking line
            if i == 0:
                ax1.plot(
                    xpos[j], [-0.2 for k in range(len(xpos[j]))], color="black", linewidth=2,
                )
    ax1.set_ylabel(labels[0])
    ax1.set_xticks(intvmidpos)
    ax1.set_xticklabels(xlabels)
    ax1.set_xlim([0, max(max(xpos)) + 1])
    ax1.set_ylim([-0.2, 0.2])
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["bottom"].set_visible(False)

    for i in range(len(ctrl_fp)):
        ctrl_fp[i] -= np.mean(ctrl_fp[i])
        for j in range(len(xpos)):
            ax2.plot(
                xpos[j],
                ctrl_fp[i][
                    sum([len(xpos[k]) for k in range(j)]) : sum([len(xpos[k]) for k in range(j)]) + len(xpos[j])
                ],
                color="black",
                marker=".",
                markersize=0.3,
                linewidth=0.3,
            )
            # draw the chromosome marking line
            if i == 0:
                ax2.plot(
                    xpos[j], [-0.2 for k in range(len(xpos[j]))], color="black", linewidth=2,
                )
    ax2.set_ylabel(labels[1])
    ax2.set_xticks(intvmidpos)
    ax2.set_xticklabels(xlabels)
    ax2.set_xlim([0, max(max(xpos)) + 1])
    ax2.set_ylim([-0.2, 0.2])
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.spines["bottom"].set_visible(False)
    plt.savefig(plotOutput)
    # plt.savefig(os.path.splitext(plotOutput)[0] + ".pdf")
    plt.close(f)


def count_bam(
    bamInput, chromsize, bedOutput, txtOutput, binlen,
):
    if not os.path.exists(bedOutput):
        divide_bin_2(chromsize, bedOutput, binlen)
    bedtool = pybedtools.BedTool(bedOutput)
    result = bedtool.multi_bam_coverage(bams=[bamInput])
    cov_result = []
    for intv in result:
        cov_result.append(intv[-1])
    cov_result = np.transpose(cov_result)
    pos = [[], []]
    bins = pybedtools.BedTool(bedOutput)
    for bin in bins:
        pos[0].append(bin.chrom)
        pos[1].append(str(bin.start + 1) + "-" + str(bin.end))
    cov_df = pd.DataFrame({"chrom": pos[0], "start-end": pos[1], "value": cov_result})
    cov_df.to_csv(txtOutput, sep="\t", header=True, index=None)
    return True


def processWPS(bedgzInput, tsvInput, protectInput, outputfile, empty, minInsSize, maxInsSize):
    if minInsSize > 0 and maxInsSize > 0 and minInsSize < maxInsSize:
        pass
    else:
        minInsSize, maxInsSize = None, None
    outputfile = outputfile.strip("""\'""")
    protection = protectInput // 2
    validChroms = set(map(str, list(range(1, 23)) + ["X", "Y"]))  # human genome
    infile = open(tsvInput)  # input transcript region file
    bedgzfile = bedgzInput.strip("""\'""")
    tbx = pysam.TabixFile(bedgzfile)
    prefix = "chr"
    # print input files
    print("input file:", bedgzInput, tsvInput)
    for line in infile.readlines():
        (cid, chrom, start, end, strand,) = line.split()  # positions should be 0-based and end non-inclusive
        chrom = chrom.replace("chr", "")
        if chrom not in validChroms:
            continue
        regionStart, regionEnd = int(start), int(end)  # this is 1-based
        if regionStart < 1:
            continue  # invalid region
        posRange = defaultdict(lambda: [0, 0])
        filteredReads = Intersecter()
        try:  # if tbx.fetch do not find any row, next row
            for row in tbx.fetch(
                prefix + chrom, regionStart - protection - 1, regionEnd + protection + 1
            ):  # all fragments overlaped with this region is collected
                tmp_row = row.split()
                rstart = int(tmp_row[1]) + 1  # convert to 1-based
                rend = int(tmp_row[2])  # end included
                lseq = int(tmp_row[2]) - int(tmp_row[1])  # fragment length
                if (
                    (minInsSize is not None)
                    and (maxInsSize is not None)
                    and ((lseq < minInsSize) or (lseq > maxInsSize))
                ):
                    continue  # satisfy length requirement
                filteredReads.add_interval(Interval(rstart, rend))  # save the fragments overlap with region
                for i in range(
                    rstart, rend + 1
                ):  # for a single nucleotide site, compute how many reads overlaped span it (include read end point)
                    if i >= regionStart and i <= regionEnd:
                        posRange[i][0] += 1
                if (
                    rstart >= regionStart and rstart <= regionEnd
                ):  # for a single nucleotide site, compute how many read end point located at this site
                    posRange[rstart][1] += 1
                if rend >= regionStart and rend <= regionEnd:
                    posRange[rend][1] += 1
        except Exception:
            continue
        filename = outputfile % cid  # name output file by names in transcript file
        outfile = gzip.open(filename, "w")
        cov_sites = 0
        outLines = []
        for pos in range(regionStart, regionEnd + 1):
            rstart, rend = pos - protection, pos + protection
            gcount, bcount = 0, 0
            for read in filteredReads.find(rstart, rend):
                if (read.start > rstart) or (read.end < rend):
                    bcount += 1  # fragments located in window
                else:
                    gcount += 1  # fragments spanned window
            covCount, startCount = posRange[pos]
            cov_sites += covCount
            # chrom: chromatin, pos: position in the genome, covCount:how many reads span this site, startCount: how many reads end point located
            # in this site, gcount-bcount: WPS
            outLines.append("%s\t%d\t%d\t%d\t%d\n" % (chrom, pos, covCount, startCount, gcount - bcount))
        if strand == "-":
            outLines = outLines[::-1]  # - strand!!!
        for line in outLines:
            outfile.write(line.encode())  # write in binary
        outfile.close()
        if cov_sites == 0 and not empty:  # remove empty files
            os.remove(filename)

    return True


def processDMR(mlInput, caselen, method, max_dif, txtOutput):
    allsum = mlInput.iloc[:, 3:].sum(axis=1).tolist()
    val = [True if allsum[i] > 0 else False for i in range(len(allsum))]
    ml_val = mlInput.iloc[np.where(np.array(val))[0], :]
    ml_val.index = pd.Series(np.arange(ml_val.shape[0]))
    casemean = ml_val.iloc[:, 3 : caselen + 3].mean(axis=1).tolist()
    ctrlmean = ml_val.iloc[:, caselen + 3 :].mean(axis=1).tolist()
    dif = [abs(casemean[i] - ctrlmean[i]) for i in range(len(casemean))]
    is_dmr = [True if dif[i] >= max_dif else False for i in range(len(dif))]
    updown = ["up" if casemean[i] >= ctrlmean[i] else "down" for i in range(len(casemean))]
    p_value = []
    for i in range(ml_val.shape[0]):
        caseml = ml_val.iloc[i, 3 : caselen + 3].values.tolist()
        ctrlml = ml_val.iloc[i, caselen + 3 :].values.tolist()
        t, p = stats.ttest_ind(caseml, ctrlml, equal_var=False)
        p_value.append(p)
    adj_res = multi.multipletests(p_value, method=method)
    marker = [
        ml_val["chr"].tolist()[i] + ":" + str(ml_val["start"].tolist()[i]) + "-" + str(ml_val["end"].tolist()[i])
        for i in range(ml_val.shape[0])
    ]
    ml = pd.concat(
        [
            pd.DataFrame({"Chr:start-end": marker}),
            ml_val.iloc[:, 3:],
            pd.DataFrame(
                {"is_DMR": is_dmr, "case_higher_or_lower": updown, "p_value": p_value, "p_value_adj": adj_res[1]}
            ),
        ],
        axis=1,
    )
    ml.to_csv(txtOutput, sep="\t", header=True, index=False)

    return True


def get_cross_validation(algorithm, X, y):
    trueLabel = np.array([])
    predLabel = np.array([])
    if len(y) > 10000:
        print("Data contains more than 10000 rows. Using 10-fold cross validation...")
        splitter = KFold(n_splits=10, shuffle=True, random_state=4294967295)  # 2**32 - 1
    else:
        print("Using leave-one-out cross validation...")
        splitter = LeaveOneOut()

    for train_indices, test_indices in splitter.split(X):
        train_X, train_y = X[train_indices], y[train_indices]
        test_X, test_y = X[test_indices], y[test_indices]

        model = algorithm.fit(train_X, train_y)
        trueLabel = np.append(trueLabel, test_y)
        predLabel = np.append(predLabel, model.predict(test_X))

    return trueLabel, predLabel


def mapBitToChar(im, col, row):
    if im.getpixel((col, row)):
        return " "
    else:
        return "#"


def logoPrint(
    mess="cfDNApipe", fontType=pkg_resources.resource_filename("cfDNApipe", "data/Verdana.ttf"), size=13,
):
    """
    from https://stackoverflow.com/questions/9632995/how-to-easily-print-ascii-art-text
    """
    ShowText = mess
    font = ImageFont.truetype(fontType, size)
    size = font.getsize(ShowText)
    image = Image.new("1", size, 1)
    draw = ImageDraw.Draw(image)
    draw.text((0, 0), ShowText, font=font)
    for r in range(size[1]):
        print("".join([mapBitToChar(image, c, r) for c in range(size[0])]))


def indexCheck(filePath=None, fileSuffix=None):
    if not os.path.exists(filePath):
        mess = "File " + filePath + " is not exist!"
        raise commonError(mess)

    index_file = filePath + fileSuffix

    if os.path.exists(index_file):
        return True
    else:
        return False


def getFileNamePrefix(fileName):
    return os.path.split(fileName)[-1]


def getMaxFileNamePrefix(file1, file2):
    file1 = getFileNamePrefix(file1)
    file2 = getFileNamePrefix(file2)
    len1 = len(file1)
    len2 = len(file2)
    for i in range(min(len1, len2)):
        if file1[i] != file2[i]:
            break
    if i == 0:
        return ""
    elif i == (min(len1, len2) - 1):
        tmp_str = file1[: i + 1]
    else:
        tmp_str = file1[:i]

    for k in reversed(range(len(tmp_str))):
        if isAlphaOrDigit(tmp_str[k]):
            final_name = tmp_str[: k + 1]
            final_name = rmEndString(final_name, [".pair", "R", "_", ".", "*"])
            return final_name
        else:
            k = k - 1

    raise commonError("File names must contain at least one alphabet or number.")


def maxCore(nCore=None):
    if nCore > 16:
        return 16
        print("The thread number is forced to 16!")
    else:
        return nCore
