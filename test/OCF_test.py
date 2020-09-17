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
import matplotlib.pyplot as plt
import gzip
import shutil
import statsmodels.api as sm
import statsmodels.stats.multitest as multi
from scipy import stats, optimize
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


inputFile = r"/opt/tsinghua/zhangwei/302_20200902/HCC/intermediate_result/step_07_bam2bed/"
refFile
txtOutput
cudOutput
ocfOutput
flags

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
    names=[
        "read.chr",
        "read.start",
        "read.end",
        "peak.chr",
        "peak.start",
        "peak.end",
        "description",
        "overlap",
    ],
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
        if (row["read.start"] < row["peak.start"]) and (
            row["peak.start"] <= row["read.end"]
        ):
            o_s, o_e = 0, row["read.end"] - row["peak.start"]
            cov[0 : (o_e + 1)] += 1
            dend[o_e] += 1
        elif (row["peak.start"] <= row["read.start"]) and (
            row["read.end"] <= row["peak.end"]
        ):
            o_s, o_e = (
                row["read.start"] - row["peak.start"],
                row["read.end"] - row["peak.start"],
            )
            cov[o_s : (o_e + 1)] += 1
            uend[o_s] += 1
            dend[o_e] += 1
        elif (row["read.start"] <= row["peak.end"]) and (
            row["peak.end"] < row["read.end"]
        ):
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



