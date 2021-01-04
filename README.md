# cfDNApipe

* [Introduction](#introduction)
* [Section 1: Installation Tutorial](#section-1-installation-tutorial)
    * [Section 1.1: System requirement](#section-11-system-requirement)
    * [Section 1.2: Create environment and Install Dependencies](#section-12-create-environment-and-install-dependencies)
    * [Section 1.3: Enter Environment and Use cfDNApipe](#section-13-enter-environment-and-use-cfdnapipe)
* [Section 2: cfDNApipe Highlights](#section-2-cfdnapipe-highlights)
    * [Section 2.1: Dataflow Graph for WGS and WGBS Data Processing](#section-21-dataflow-graph-for-wgs-and-wgbs-data-processing)
    * [Section 2.2: Reference Auto Download and Building](#section-22-reference-auto-download-and-building)
    * [Section 2.3: Output Folder Arrangement](#section-23-output-folder-arrangement)
    * [Section 2.4: Analysis Report](#section-24-analysis-report)
    * [Section 2.5: Breakpoint Detection](#section-25-breakpoint-detection)
    * [Section 2.6: Other Mechanisms](#section-26-other-mechanisms)
* [Section 3: A Quick Tutorial for Analysis WGBS data](#section-3-a-quick-tutorial-for-analysis-wgbs-data)
    * [Section 3.1: Set Global Reference Configures](#section-31-set-global-reference-configures)
    * [Section 3.2: Execute build-in WGBS Analysis Pipeline](#section-32-execute-build-in-wgbs-analysis-pipeline)
* [Section 4: Perform Case-Control Analysis for WGBS data](#section-4-perform-case-control-analysis-for-wgbs-data)
* [Section 5: How to Build Customized Pipepline using cfDNApipe](#section-5-how-to-build-customized-pipepline-using-cfdnapipe)
* [Section 6: Additional function: WGS SNV/InDel Analysis](#section-6-additional-function-wgs-snvindel-analysis)
    * [Section 6.1: Reference Files Preparation](#section-61-reference-files-preparation)
    * [Section 6.2: Performing Single Group Samples SNV Analysis](#section-62-performing-single-group-samples-snv-analysis)
    * [Section 6.3: Performing Case-Control SNV Analysis](#section-63-performing-case-control-snv-analysis)
* [Section 7: Additional Function: Virus Detection](#section-7-additional-function-virus-detection)

## Introduction

**cfDNApipe(<u>c</u>ell <u>f</u>ree <u>DNA</u> <u>Pipe</u>line)** is an integrated pipeline for analyzing [cell-free DNA](https://en.wikipedia.org/wiki/Circulating_free_DNA) WGBS/WGS data. It contains many cfDNA quality control and feature extration algorithms. Also we collected some useful cell free DNA references and provide them [here](https://honchkrow.github.io/cfDNAReferences/).

The whole pipeline was established based on processing graph principle. Users can use the inside integrated pipeline for WGBS/WGS data as well as build their own analysis pipeline from any intermediate data like bam files. The main functions are as the following picture.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./pics/pipeline.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">cfDNApipe Functions</div>
</center>


## Section 1: Installation Tutorial

### Section 1.1: System requirement

The popular WGBS/WGS analysis softwares are released on Unix/Linux system, based on different program language, like FASTQC and Bowtie2. Therefore, it's very difficult to rewrite all the software in one language. Fortunately, [conda](https://docs.conda.io/en/latest/)/[bioconda](http://bioconda.github.io/) program collected many prevalent python mudules and bioinformatics software, so we can install all the dependencies through [conda](https://docs.conda.io/en/latest/)/[bioconda](http://bioconda.github.io/) and arrange pipelines using python.

We recommend using conda and create an virtual environment to manage all the dependencies. If you did not install conda before, please follow [this tutorial](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) to install conda first.

After installation, you can create a new virtual environment for cfDNA analysis. Virtual environment management means that you can install all the dependencies in this virtual environment and delete them easily by removing this virtual environment.

### Section 1.2: Create environment and Install Dependencies

We tested our pipeline using different version of software and provide an  environment yml file for users. Users can download this file and create the environment in one command line. 

First, please download the yml file.

``` shell
wget https://raw.githubusercontent.com/Honchkrow/cfDNApipe/master/environment.yml
```

Then, run the following command. The environment will be created and all the dependencies as well as the latest cfDNApipe will be installed.

``` shell
conda env create -n cfDNApipe -f environment.yml
```

*<font color=red>Note:</font> The environment name can be changed by replacing "-n cfDNApipe" to "-n environment_name".*

### Section 1.3: Enter Environment and Use cfDNApipe

Once the environment is created, user can enter environment using the foloowing command.

``` shell
conda activate cfDNApipe
```

 Now, just open python and process **cell free DNA WGBS/WGS paired/single end** data. For more detailed explanation for each function and parameters, please see [here](***************).


## Section 2: cfDNApipe Highlights

cfDNApipe is a highly integrated cfDNA related WGS/WGBS data processing pipeline. We designed many useful build-in mechanism. Here, we introduce some of them to the users.

### Section 2.1: Dataflow Graph for WGS and WGBS Data Processing

cfDNApipe is organized by a built-in dataflow with strictly defined up- and down-stream data interface. The following two figures shows the up- and down-stream relationship between functions in WGS and WGBS data processing.

<br/>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./pics/WGBS_pipeline.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">WGBS Dataflow Overview</div>
</center>

<br/>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./pics/WGS_pipeline.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">WGS Dataflow Overview</div>
</center>

### Section 2.2: Reference Auto Download and Building

For any HTS data analysis, the initial step is to set reference files like genome sequence and annotation files. Here, we introduced global reference configure function in cfDNApipe to download and build reference files automatically.

cfDNApipe contains 2 types of global reference configure function, **pipeConfigure** and **pipeConfigure2**. Function **pipeConfigure** is for single group data analysis (without control group). Function **pipeConfigure2** is for case and control analysis. Either function will check the reference files, such as bowtie2 and bismark references. If not detected, references will be downloaded and built. This step is **<font color=red>necessary</font>** and puts things right once and for all.

*<font color=red>Note:</font> Users should use the correct configure function **pipeConfigure** and **pipeConfigure2**. The output folder arrangement stategy is totally different for these two function. In addition, some default files can only be accessed through pipeConfigure or pipeConfigure2. Therefore, if a single group data analysis is needed, using **pipeConfigure**. If a case-control comparison analysis is needed, using **pipeConfigure2**. If users want to switch analysis from single group to case-control group and vice versa, the customized pipeline can achieve the seamless linking between output and input of different functions.*

The folowing is a simple **pipeConfigure** example for building WGBS reference files.

```Python
from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"path_to_reference/hg19_bismark",
    outdir=r"path_to_output/WGBS",
    data="WGBS",
    type="paired",
    JavaMem="10g",
    build=True,
)

```

pipeConfigure function takes 8 necessary parameters as input. 

* **'threads'**: the max threads user want to be used. 
* **'genome'**: which genome to be used, must be 'hg19' or 'hg38'. 
* **'refdir'**: where to find genome reference files like sequence fasta file and CpG island ananotation files. 
* **'outdir'**: where to put all the output files.
* **'data'**: "WGBS" or "WGS".
* **'type'**: "paired" or "single".
* **'JavaMem'**: maximum memory java program can used.
* **'build'**: download and build reference or not after reference checking.

Like the above example, if refdir is empty, cfDNApipe will download hg19.fa and other annotation files automatically. Once done, the program will print "Background reference check finished!", then users can do the analyzing steps.

We also provide a simple **pipeConfigure2** example for building WGBS reference files.

```Python
from cfDNApipe import *

pipeConfigure2(
    threads=60,
    genome="hg19",
    refdir=r"path_to_reference/hg19_bismark",
    outdir=r"path_to_output/WGBS",
    data="WGBS",
    type="paired",
    case="cancer",
    ctrl="normal",
    JavaMem="10g",
    build=True,
)


```

*<font color=red>Note:</font> The download procudure is always time-consuming. Therefore cfDNApipe can detect the reference files which are already existed in refdir. For instance, users can just put hg19.fa file into refdir and cfDNApipe will not download it again. Other reference files can be got from [here](https://github.com/Honchkrow/cfDNAReferences). Downlaoding, uncompressing and putting them into refdir will be much faster.*


### Section 2.3: Output Folder Arrangement

Generally, the cell free DNA analysis contains many steps, which will generate lots of output files. cfDNApipe arrange the outputs into every functinal specific folders. Based on analysis stategy (with or without control), the output folders are arranged as follows.

- Analysis Results Without Control Samples

``` 
output_folders/  
├── final_result/  
├── report_result/  
│   ├── Cell_Free_DNA_WGBS_Analysis_Report.html  
│   └── Other files and folders  
└── intermediate_result/  
    ├── step_01_inputprocess  
    ├── step_02_fastqc  
    ├── step_02_identifyAdapter  
    └── Other processing folders  
```

- Analysis Results With Control Samples (assume case and control)

``` 
output_folders/
├──case/
│   ├── final_result/  
│   ├── report_result/  
│   │   ├── Analysis_Report.html  
│   │   └── Other files and folders  
│   └── intermediate_result/  
│       ├── step_01_inputprocess  
│       ├── step_02_fastqc  
│       ├── step_02_identifyAdapter  
│       └── Other processing folders  
├──control/
    ├── final_result/  
    ├── report_result/  
    │   ├── Analysis_Report.html  
    │   └── Other files and folders  
    └── intermediate_result/  
        ├── step_01_inputprocess  
        ├── step_02_fastqc  
        ├── step_02_identifyAdapter  
        └── Other processing folders  
```

There will be 3 major ouput folder for every sample group, named **"final_result"**, **"report_result"**, and **"intermediate_result"**. 

Folder **"intermediate_result"** contains folders named by every single step, all the intermediate results and processing record will be save in each folder. User can accsee any files they want. This folder is evry large since all the intermediate files are saved in this folder. Users can move some results to the folder **"final_result"** and deleted **"intermediate_result"** after all the analysis is finished.

Folder **"report_result"** save a pretty html report and related data which shows some visualization results like quality control and analysis figures. 

Folder **"final_result"** is an empty folder for users to save specific results from intermediate_result folder. 


### Section 2.4: Analysis Report

Folder **"report_result"** can visualize analysis results, like DNA fragment length distribution and mapping statistics. The report folder can be copied to any where. Here is an [example]() showing the final report.

We try our best to plot every figure ready to publish. If users want to make some changes like changing colors, they can access figure data saved at every step foler in  **"intermediate_result"**.


### Section 2.5: Breakpoint Detection

Sometimes, the program may be interrupted by irresistible reasons like computer crash. cfDNApipe provide **breakpoint detection mechanism**, which compute md5 code for inputs, outputs, as well as all parameters. Therefore, user do not warry about any interrupt situation. Re-running the same program, the finished step will show message like below and be skipped automatically.

``` shell
************************************************************
                bowtie2 has been completed!
************************************************************
```

### Section 2.6: Other Mechanisms

* Parallel Computing
* Case and Control Analysis
* Inputs Legality Checking
* ......

## Section 3: A Quick Tutorial for Analysis WGBS data

In this section, we will demonstrate how to perform a quick analysis for paired end WGBS data using the build-in pipeline.

### Section 3.1: Set Global Reference Configures

First, user must set some important configure, for example, which genome to be used, how many threads should be used and where to put the analysis results. cfDNApipe provide a configure function for user to set these parameters. Below is an instance.

``` Python
from cfDNApipe import *

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"path_to_reference/hg19_bismark",
    outdir=r"path_to_output/WGBS",
    data="WGBS",
    type="paired",
    build=True,
    JavaMem="10g",
)
```



### Section 3.2: Execute build-in WGBS Analysis Pipeline

cfDNApipe provides an integrated pipeline for paired/single end WGBS/WGS data, user can use it easily by assigning fastq sequencing files as the input of the pipeline. All the parameters used in pipeline are carefully selected after numerous tests.

``` Python
res = cfDNAWGBS(inputFolder=r"path_to_fastqs",
                idAdapter=True,
                rmAdapter=True,
                dudup=True,
                CNV=False,
                armCNV=False,
                fragProfile=False,
                report=True,
                verbose=False)
```

In the above example, user just pass the input folder which contains all the raw fastq files to the function, then the processing will start and all results will be saved in output folder mentioned in the former section. What's more, "report=True" will generate a html report for users.

In addition, cfDNApipe also provides **case-control** comparison analysis for WGBS/WGS data. For using this function, please see the section 4 and function **cfDNAWGS2** and **cfDNAWGBS2**.


## Section 4: Perform Case-Control Analysis for WGBS data

The analysis steps for case-control analysis are the same as section 2.1 and 2.2. First, set global configure. Second, run analysis command.

Setting global configure is a little bit different from section 2.1. Below is an example.

``` Python
from cfDNApipe import *

pipeConfigure2(
    threads=60,
    genome="hg19",
    refdir=r"path_to_reference/hg19_bismark",
    outdir=r"path_to_output/WGBS",
    data="WGBS",
    type="paired",
    JavaMem="10G",
    case="cancer",
    ctrl="normal",
    build=True,
)
```

Here, 2 more parameters are used. Parameter **"case"** and **"ctrl"** is the name flag for case and control data. These two parameters control the output for case and control samples.

Next, using function **cfDNAWGBS2** to processing case and control analysis.

``` Python
case, ctrl, comp = cfDNAWGBS2(
    caseFolder=r"case_fastqs",
    ctrlFolder=r"ctrl_fastqs",
    caseName="cancer",
    ctrlName="normal",
    idAdapter=True,
    rmAdapter=True,
    dudup=True,
    armCNV=True,
    CNV=True,
    fragProfile=True,
    verbose=False,
)
```

After analysis, user can get all the output as well as reports for case and control. Of course, the comparison results will be saved in **case folder**.

## Section 5: How to Build Customized Pipepline using cfDNApipe

Some users are familiar with cfDNA processing and want to customize their own pipelines. cfDNApipe provide a flexible pipeline framework for building customized pipeline. The following is an example of how to build pipeline from intermediate steps.

Assume that we have some WGS samples and all the samples have already been aligned. Now, we want to perform CNA analysis to these data compared with default sequence and get gene level annotation.

First, set global configure.

``` Python
from cfDNApipe import *
import glob

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"path_to_reference/hg19_bowtie2",
    outdir=r"path_to_output/WGS",
    data="WGS",
    type="paired",
    JavaMem="10G",
    build=True,
)
```

Second, take all the aligned bam as inputs.

``` Python
# get all bam files
bams = glob.glob("samples/*.bam")

# sort bam and remove duplicates
res_bamsort = bamsort(bamInput=bams, upstream=True)
res_rmduplicate = rmduplicate(upstream=res_bamsort)

# perform CNV analysis, stepNum is set as flag of every step
res_cnvbatch = cnvbatch(
    caseupstream=res_rmduplicate,
    access=Configure.getConfig("access-mappable"),
    annotate=Configure.getConfig("refFlat"),
    stepNum="CNV01",
)
res_cnvPlot = cnvPlot(upstream=res_cnvbatch, stepNum="CNV02")
res_cnvTable = cnvTable(upstream=res_cnvbatch, stepNum="CNV03")
res_cnvHeatmap = cnvHeatmap(upstream=res_cnvbatch, stepNum="CNV04")
```

In the above codes, **"upstream=True"** means puts all the results to the output folder mentioned in section 3.1. CNV analysis needs two reference files is control samples are not provided. These two reference files are already included in cfDNApipe reference data, user can access them easily.

Once finished, user can get CNV related files like below.

## Section 6: Additional function: WGS SNV/InDel Analysis

*<font color=red>Note:</font> This function is <font color=red>only</font> supported for processing WGS data.*

We wrapped classical software [**GATK4**](https://gatk.broadinstitute.org/hc/en-us) to call WGS mutations. Detecting mutations needs addtional references related to human genome. These references are provided by [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) and not suit for auto-downloading. Therfore, users should download the reference files manually. [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) provides different ways to download reference files like lftp and ftp. We recommend using **lftp** to download the VCF references for convenient.

If **lftp** is not installed, users can install it from conda easily.

```shell
conda install -c conda-forge lftp
```

Before the analysis, we recommend users to create a new folder for saving the snv related referencef files. For example, create a folder in your genome reference folder and name it based on the genome version like hg19_snv or hg38_snv. Then, enter the folder to download the snv reference files.

### Section 6.1: Reference Files Preparation

Users can login in GATK resource bundle and download the dependent VCF files.

```shell
lftp gsapubftp-anonymous@ftp.broadinstitute.org:/
#just click "Enter" for password
lftp gsapubftp-anonymous@ftp.broadinstitute.org:/> ls
```

Data hosted in GATK resource bundle is shown as follows.

```shell
drwxrwxr-x  12 ebanks   wga           235 Jul 30  2018 bundle
drwxrwxr-x   2 ebanks   wga           219 Nov 28  2011 DePristoNatGenet2011
drwxr-xr-x   2 droazen  broad         253 Oct 24  2014 foghorn140
-rw-r--r--   1 vdauwera broad     2141164 Jul 15  2014 gatkdocs-3_1_v_3_2.zip
drwxrwxr-x   2 ebanks   wga            34 Aug 12  2011 HLA
drwxrwxr-x   2 ebanks   wga           293 May 21  2019 Liftover_Chain_Files
-rw-r--r--   1 vdauwera broad    11645558 Aug  3  2017 mutect-1.1.7.jar
-rw-r--r--   1 vdauwera broad    10547264 Aug  3  2017 mutect-1.1.7.jar.zip
drwxrwxr-x   2 thibault broad         169 May 11  2015 travis
drwxrwxr-x   4 ebanks   wga            71 Aug  7  2016 triosForGA4GH
-rw-rw-r--   1 vdauwera broad    739681240 Oct 22  2013 tutorial_files.zip
drwxrwxr-x   7 vdauwera broad         126 Aug  3  2017 tutorials
drwxrwxr-x   5 vdauwera broad          72 Jun 24  2016 workshops
```

Then download the dependent files based on related genome version.

For **hg19**:

*<font color=red>Note:</font> Some files for hg19 is not provided by GATK, therefore we should convert them from b37 version. cfDNApipe will do the conversion <font color=red>automatically</font>.*

```shell
glob -- pget -c -n 12 bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz
glob -- pget -c -n 12 bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
glob -- pget -c -n 12 bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
glob -- pget -c -n 12 bundle/hg19/dbsnp_138.hg19.vcf.gz
glob -- pget -c -n 12 bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
glob -- pget -c -n 12 bundle/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz
glob -- pget -c -n 12 bundle/Mutect2/GetPileupSummaries/small_exac_common_3_b37.vcf.gz
glob -- pget -c -n 12 Liftover_Chain_Files/b37tohg19.chain
```

<br/>

For **hg38**:

```shell
glob -- pget -c -n 12 bundle/hg38/1000G_omni2.5.hg38.vcf.gz
glob -- pget -c -n 12 bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
glob -- pget -c -n 12 bundle/hg38/dbsnp_146.hg38.vcf.gz
glob -- pget -c -n 12 bundle/hg38/hapmap_3.3.hg38.vcf.gz
glob -- pget -c -n 12 bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
glob -- pget -c -n 12 bundle/Mutect2/af-only-gnomad.hg38.vcf.gz.gz
glob -- pget -c -n 12 bundle/Mutect2/GetPileupSummaries/small_exac_common_3.hg38.vcf.gz
```

For analyzing **single group samples**, we need [Public GATK Panel of Normals](https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-) (PON) file. Here we provide the google storage download address directly.

For **hg19**, [somatic-b37_Mutect2-WGS-panel-b37.vcf](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-b37;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).


<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./pics/hg19_snv_download.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">click the dwonload button in red mark</div>
</center>

<br/>

For **hg38**, [somatic-hg38_1000g_pon.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./pics/hg38_snv_download.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">click the dwonload button in red mark</div>
</center>

<br/>

Just put the PON file in the same folder as the other snv reference files download before.

Finally, uncompressing all the files.

``` shell 
gunzip *.gz
```

### Section 6.2: Performing Single Group Samples SNV Analysis

When finish preparing all the files, we can use "snvRefCheck" function in cfDNApipe to achieve genome version conversion from b37 to hg19 and indexing. Here we use hg19 to show the demo for single group samples SNV analysis.

The following are the whole scripts for single group samples SNV analysis with simple annotations. For detailed explanation, please "help" function in python.

```Python
from cfDNApipe import *
import glob

pipeConfigure(
    threads=60,
    genome="hg19",
    refdir=r"path_to_reference/hg19",
    outdir=r"path_to_output/snv_output",
    data="WGS",
    type="paired",
    build=True,
    JavaMem="10G",
)

# just set build=True to finish all the works
Configure.snvRefCheck(folder="path_to_reference/hg19/hg19_snv", build=True)

# see all the SNV related files
Configure.getConfig("snv.folder")
Configure.getConfig("snv.ref")

# indexed bam files after remove duplicates
bams = glob.glob("path_to_samples/*.bam")

res1 = addRG(bamInput=bams, upstream=True)

res2 = BaseRecalibrator(
    upstream=res1, knownSitesDir=Configure.getConfig("snv.folder")
)
res3 = BQSR(upstream=res2)

res4 = getPileup(
    upstream=res3,
    biallelicvcfInput=Configure.getConfig('snv.ref')["7"],
)

res5 = contamination(upstream=res4)

# In this step, files are split to chromatin
res6 = mutect2t(
    caseupstream = res5,
    vcfInput=Configure.getConfig('snv.ref')["6"],
    ponbedInput=Configure.getConfig('snv.ref')["8"],
)

# 
res7 = filterMutectCalls(upstream=res6)

# put all chromatin files together
res8 = gatherVCF(upstream=res7)

# split somatic mutations
res9 = bcftoolsVCF(upstream=res8, stepNum="somatic")

# split germline mutations
res10 = bcftoolsVCF(
    upstream=res8, other_params={"-f": "'germline'"}, suffix="germline", stepNum="germline"
)

```


*<font color=red>Note:</font> User can adjust the parameter <font color=red>"--f-score-beta"</font> in function <font color=red>filterMutectCalls</font> for a very strict filtering. For detailed information, please see [filterMutectCalls manual](https://gatk.broadinstitute.org/hc/en-us/articles/360037225412-FilterMutectCalls).*

The output vcf file from function <font color=blue>bcftoolsVCF</font> can be annotated by other software such as [annovar](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/).

### Section 6.3: Performing Case-Control SNV Analysis

Here we also use hg19 to show the demo for case-control SNV analysis.

```Python
from cfDNApipe import *
import glob

# set global configure
pipeConfigure2(
    threads=100,
    genome="hg19",
    refdir=r"path_to_reference/hg19",
    outdir=r"path_to_output/snv_output",
    data="WGS",
    type="paired",
    case="cancer",
    ctrl="normal",
    JavaMem="10g",
    build=True,
)

Configure2.snvRefCheck(folder="path_to_reference/hg19/SNV_hg19", build=True)

case_bams = glob.glob("path_to_samples/cancer/*.bam")
ctrl_bams = glob.glob("path_to_samples/normal/*.bam")

# create PON file from normal samples
switchConfigure("normal")

ctrl_addRG = addRG(bamInput=ctrl_bams, upstream=True, stepNum="PON00",)
ctrl_BaseRecalibrator = BaseRecalibrator(
        upstream=ctrl_addRG,
        knownSitesDir=Configure2.getConfig("snv.folder"),
        stepNum="PON01",
        )
ctrl_BQSR = BQSR(upstream = ctrl_BaseRecalibrator, stepNum="PON02",)
ctrl_mutect2n = mutect2n(upstream = ctrl_BQSR, stepNum="PON03",)
ctrl_dbimport = dbimport(upstream = ctrl_mutect2n, stepNum="PON04",)
ctrl_createPon = createPON(upstream = ctrl_dbimport, stepNum="PON05",)

# performing comparison between cancer and normal 
switchConfigure("cancer")
case_addRG = addRG(bamInput=case_bams, upstream=True, stepNum="SNV00",)
case_BaseRecalibrator = BaseRecalibrator(
    upstream=case_addRG,
    knownSitesDir=Configure2.getConfig("snv.folder"),
    stepNum="SNV01",
)

case_BQSR = BQSR(
    upstream=case_BaseRecalibrator, 
    stepNum="SNV02")

case_getPileup = getPileup(
    upstream=case_BQSR,
    biallelicvcfInput=Configure2.getConfig('snv.ref')["7"],
    stepNum="SNV03",
)
case_contamination = contamination(
    upstream=case_getPileup,  
    stepNum="SNV04"
)

# In this step, ponbedInput is ignored by using caseupstream parameter
case_mutect2t = mutect2t(
    caseupstream=case_contamination,
    ctrlupstream=ctrl_createPon,
    vcfInput=Configure2.getConfig('snv.ref')["6"],
    stepNum="SNV05",
)

case_filterMutectCalls = filterMutectCalls(
    upstream=case_mutect2t,
    stepNum="SNV06"
)

case_gatherVCF = gatherVCF(
    upstream=case_filterMutectCalls, 
    stepNum="SNV07"
)

# split somatic mutations
case_somatic = bcftoolsVCF(upstream=case_gatherVCF, stepNum="somatic")

# split germline mutations
case_germline = bcftoolsVCF(
    upstream=case_gatherVCF, other_params={"-f": "'germline'"}, suffix="germline", stepNum="germline"
)


```

*<font color=red>Note:</font> User can adjust the parameter <font color=red>"--f-score-beta"</font> in function <font color=red>filterMutectCalls</font> for a very strict filtering. For detailed information, please see [filterMutectCalls manual](https://gatk.broadinstitute.org/hc/en-us/articles/360037225412-FilterMutectCalls).*


The output vcf file from function <font color=blue>bcftoolsVCF</font> can be annotated by other software such as [annovar](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/).

## Section 7: Additional Function: Virus Detection

*<font color=red>Note:</font> This function is <font color=red>only</font> supported for processing WGS data.*

cfDNApipe wraps [centrifuge](https://ccb.jhu.edu/software/centrifuge/) to detect virus. Centrifuge is a very rapid metagenomic classification toolkit. The unmapped reads are used for virus detection. 

Virus detection needs additional virus genome reference for DNA read classification. Downloading and building reference is time-consuming. Therefore, we provide extra function to tackle this problem. Users can do this in a single command as follows.

```Python
# set global parameters
from cfDNApipe import *
import glob

pipeConfigure(
    threads=20,
    genome="hg19",
    refdir=r"path_to_reference/hg19",
    outdir=r"path_to_output/virus_output",
    data="WGS",
    type="paired",
    build=True,
    JavaMem="10g",
)

# Download and Build Virus Genome
Configure.virusGenomeCheck(folder="path_to_reference/virus_database", build=True)

```

Method "virusGenomeCheck" will download and build virus reference automatically. Then, we can do virus detection.

```Python
# paired data
fq1 = glob.glob("path_to_unmapped/*.fq.1.gz")
fq2 = glob.glob("path_to_unmapped/*.fq.2.gz")

virusdetect(seqInput1=fq1, seqInput2=fq2, upstream=True)
```

The output for every sample will be 2 files. One file with suffix "output" saves the classification results for every unmapped read. Another file with suffix "report" reports statistics for virus detection like below.

| name | taxID | taxRank | genomeSize | numReads | numUniqueReads | abundance |
| :----: | :----: | :----: | :----: | :----: | :----: | :----: |
| Human endogenous retrovirus K113 | 166122 | leaf | 9472 | 257 | 257 | 0.00283657 |
| Escherichia virus phiX174 | 10847 | species | 5386 | 3813 | 3809 | 0.997073 |

