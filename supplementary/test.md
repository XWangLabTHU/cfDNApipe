
# 4 cfDNApipe Highlights

&emsp;We designed many useful build-in mechanism in cfDNApipe. Here, we introduce some of them to the users. 

<br />


## 4.1 Output Folder Arrangement

&emsp;Generally, the cell free DNA analysis contains many steps, which will generate lots of output files. cfDNApipe arrange the outputs into every functinal specific folders. Based on analysis stategy (with or without control), the output folders are arranged as follows.

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

&emsp;There will be 3 major ouput folder, named **"final_result"**, **"report_result"**, and **"intermediate_result"**. 

&emsp;Folder **"final_result"** is an empty folder for users to save any result for this analysis. 

&emsp;Folder **"report_result"** save a pretty html report and related data which shows some visualization results like quality control and analysis figures. 

&emsp;Folder **"intermediate_result"** contains folders named by every single step, all the intermediate results and processing record will be save in each folder. User can accsee any files they want.

## 4.2 Analysis Report

&emsp;Folder **"report_result"** many visable analysis results, like DNA fragment length distribution and mapping statistics. The report folder can be copied to any where. Here is an [example]() showing the final report.

&emsp;We try our best to plot every figure ready to publish. If users want to make some changes like changing colors, they can access figure data saved at every step foler in  **"intermediate_result"**.

### Section 3.3: Reference Auto Download and Build Function

&emsp;In the section2.2, We introduced global reference configure function, in which parameter **'build'** means whether to download and build references.

&emsp;cfDNApipe contains 2 type of global reference configure function, **pipeConfigure** and **pipeConfigure2**. Function **pipeConfigure** is for single type data analysis. Function **pipeConfigure2** is for case and control analysis. Either function will check the reference files, such as bowtie2 and bismark references. If not detected, references will be downloaded and built. For example, if human genome 'hg19' is specified and there is no this reference genome file in refdir, then hg19.fa will be downloaded from UCSC and other annotation files will be downloaded from [cfDNAReferences](https://honchkrow.github.io/cfDNAReferences/). 

&emsp;This step is **necessary** but put things right once and for all. If user already build references for Bismark (a folder contains Bisulfite_Genome and hg19.fa), then just set this folder as refdir, the program will **skip** download hg19.fa and rebuild Bismark reference. cfDNApipe will only download other references and this will save lots of times.

### Section 3.4: Breakpoint Detection

Sometimes, the program may be interrupted by irresistible reasons like computer crash. cfDNApipe provide **breakpoint detection mechanism**, which compute md5 code for inputs, outputs, as well as all parameters. Therefore, user do not warry about any interrupt situation. Re-running the same program, the finished step will show finished message like below.

``` shell
************************************************************
                bowtie2 has been completed!
************************************************************
```

### Section 3.5: Other Mechanisms

* Parallel Computing
* Memory Control
* Dataflow Graph
* Case and Control Analysis
* Numerous QC functions
* Inputs Legality Checking
* ......




