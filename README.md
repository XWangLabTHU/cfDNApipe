# cfDNApipe
Cell Free DNA Sequencing Analysis Pipeline

## Section 1: Installation Tutorial

### Section 1.1: System requirement
Unix/Linux system, conda environment and python >= 3.6

you can create a new virtual environment using the following command, it will not influence your default invironment and be removed easily:

```shell
# create a conda environment named cfDNApipe and install python3.6
conda create -n cfDNApipe python=3.6
# enetr the environment
conda activate cfDNApipe
```

### Section 1.2: Install Dependencies
Please download this repository and put it in your working directory.

```shell
git clone https://github.com/Honchkrow/cfDNApipe.git
```

Then, run the following command to install the dependencies.

```shell
cd cfDNApipe
chmod +x sysCheck
./sysCheck
```

### Section 1.3: Install cfDNAipe
Install cfDNApipe module.

```shell
pip install ./dist/cfDNApipe-0.0.4.tar.gz
```

## Section 2: WGBS Data Pipeline Demo and Parameters
```Python
from cfDNApipe import *

Configure.setGenome("hg19")
Configure.setRefDir(r'/home/wzhang/genome/hg19_bismark')
Configure.setThreads(20)
Configure.setOutDir(r'/data/wzhang/pipeline-test')
Configure.pipeFolderInit()

res1 = inputprocess(inputFolder = r"/data/wzhang/pipeline-test/raw-data")
res2 = fastqc(upstream = res1)
res3 = identifyAdapter(upstream = res1, formerrun = res2)
res4 = adapterremoval(upstream = res3)
res5 = bismark(upstream = res4)
res6 = bamsort(upstream = res5)
res7 = rmduplicate(upstream = res6)
res8 = bam2bed(upstream = res7)
res9 = fraglenplot(upstream = res8)
```
