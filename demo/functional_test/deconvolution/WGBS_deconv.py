# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:38:07 2019

@author: zhang
"""

import glob
from cfDNApipe import *

pipeConfigure2(
    threads=100,
    genome="hg19",
    refdir="/home/zhangwei/Genome/hg19_bismark",
    outdir="/opt/tsinghua/zhangwei/WGBS_single/test",
    data="WGBS",
    type="single",
    case="HCC",
    ctrl="Healthy",
    JavaMem="10G",
    build=True,
)

hcc = glob.glob("/opt/tsinghua/zhangwei/WGBS_single/HCC/intermediate_result/step_06_compress_methyl/*.gz")
ctr = glob.glob("/opt/tsinghua/zhangwei/WGBS_single/Healthy/intermediate_result/step_06_compress_methyl/*.gz")

verbose = False

switchConfigure("HCC")
hcc1 = calculate_methyl(tbxInput=hcc, bedInput="plasmaMarkers_hg19.bed", upstream=True, verbose=verbose)
hcc2 = deconvolution(upstream=hcc1)



switchConfigure("Healthy")
ctr1 = calculate_methyl(tbxInput=ctr, bedInput="plasmaMarkers_hg19.bed", upstream=True, verbose=verbose)
ctr2 = deconvolution(upstream=ctr1)


