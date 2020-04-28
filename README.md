# cfDNApipe

[![Build Status](https://travis-ci.org/pages-themes/minimal.svg?branch=master)](https://travis-ci.org/pages-themes/minimal) [![Gem Version](https://badge.fury.io/rb/jekyll-theme-minimal.svg)](https://badge.fury.io/rb/jekyll-theme-minimal)

cfDNApipe(<u>C</u>ell <u>F</u>ree <u>DNA</u> <u>Pipe</u>line) is an integrated pipeline for analyzing [cell-free DNA](https://en.wikipedia.org/wiki/Circulating_free_DNA) WGBS/WGS data. It contains many cfDNA quality control and feature extration algorithms. Also we collected some useful cell free DNA references and provide them [here](https://honchkrow.github.io/cfDNAReferences/).

The whole pipeline was established based on processing graph principle. Users can use the inside integrated pipeline for WGBS/WGS data as well as build their own analysis pipeline from any intermediate data like bam files. The main functions are as the following picture.

![cfDNApipe Functions](./cfDNApipe_picture1.jpg)

## Section 1: Installation Tutorial

### Section 1.1: System requirement

The popular WGBS/WGS analysis softwares are released on Unix/Linux system, based on different program language, like Bowtie2/Bismark and picard. Therefore, it's very difficult to rewrite all the software in one language. Fortunately, [conda](https://docs.conda.io/en/latest/)/[bioconda](http://bioconda.github.io/) program collected many prevalent python mudules and bioinformatics software, so we can install all the dependencies through [conda](https://docs.conda.io/en/latest/)/[bioconda](http://bioconda.github.io/) and arrange pipelines using python.

We recommend using conda environment and python >= 3.6. If you did not install conda before, please follow [this tutorial](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) to install conda first.

After installation, you can create a new virtual environment for cfDNA analysis. Virtual environment management means that you can install all the dependencies in this virtual environment and delete them easily by removing this virtual environment.

you can create a new virtual environment using the following command:

```shell
# create a conda environment named cfDNApipe with python3.6
conda create -n cfDNApipe python=3.6

# enter the environment
conda activate cfDNApipe
```

### Section 1.2: Install Dependencies

Please download this repository and put it in your working directory.

```shell
wget https://raw.githubusercontent.com/Honchkrow/cfDNApipe/master/sysCheck.sh
```

Then, run the following command and following the instruction to install the dependencies.

```shell
chmod +x sysCheck.sh
./sysCheck.sh
```

If your computer fulfills the requirement, you will see the following message.

```shell
The environment configuration is done!
```

This step will install the following software and path packages:

+ FASTQC
+ Bowtie2
+ Bismark
+ Picard
+ Samtools
+ Bedtools
+ AdapterRemoval
+ hmmcopy_utils
+ pysam
+ numpy
+ matplotlib
+ pybedtools
+ yattag
+ statsmodels
+ seaborn
+ sklearn

### Section 1.3: Install cfDNApipe

cfDNApipe can be downloaded from pypi, users can install easily by running the following command in "cfDNApipe" environment.

```shell
pip install cfDNApipe
```

Once the package is installed, user can enter python or write scripts to processing cell free DNA WGBS/WGS data.

## Section 2: A Quick Tutorial for Analysis WGBS data

In this section, we will demonstrate how to perform a quick analysis for paired end WGBS data using the build-in pipeline.

### Section 2.1: Set Global Reference Configures

First, user must set some important configure, for example, which genome to be used, how many threads should be used and where to put the analysis results. cfDNApipe provide a configure function for user to set these parameters. Below is an instance.

```Python
from cfDNApipe import *

pipeConfigure(
    threads=20,
    genome="hg19",
    refdir=r"./genome/hg19_bismark",
    outdir=r"./pipeline-for-paired-WGBS",
    data="WGBS",
    type="paired",
    build=True,
)
```

pipeConfigure function 7 necessary parameters as input. 
Parameter 'threads' defines the max threads user want to use. 
Parameter 'genome' shows which genome to be used, must be 'hg19' or 'hg38'. 
'refdir' means where to find genome reference files like sequence fasta file and CpG island ananotation files. 

.
├── pipeline-for-paired-WGBS/
    ├── final_result/
    ├── report_result/
    |   ├── Cell_Free_DNA_WGBS_Analysis_Report.html
    |   └── Other files and folders
    └── intermediate_result/
        ├── step_01_inputprocess
        ├── step_02_fastqc
        ├── step_02_identifyAdapter
        └── Other processing folders

There will be 3 major ouput folder, named final_result, report_result and intermediate_result. Folder 'final_result' is an empty folder for users to save any result for this analysis. Folder 'report_result' save a html report and related data which shows some visualization results like quality control and figures. Folder 'intermediate_result' contains many folder named by every single step, all the intermediate results and processing record will be save in each folder.

Parameter 'data' and 'type' show data type and sequencing type.
Parameter 'build' means whether to download and build references. For example, if human genome 'hg19' is specified and there is no this reference genome file in refdir, then hg19.fa will be downloaded from UCSC and other annotation files will be downloaded from [cfDNAReferences](https://honchkrow.github.io/cfDNAReferences/). This step is necessary but put things right once and for all. If user already build references for Bismark (a folder contains Bisulfite_Genome and hg19.fa), then just set this folder as refdir, the program will skip download hg19.fa and rebuild Bismark reference. cfDNApipe will only download other references and this will save lots of times.

### Section 2.2: Perform build-in WGBS Analysis Pipeline


