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

Then, run the following command and following the instrucion to install the dependencies.

```shell
cd cfDNApipe
chmod +x sysCheck
./sysCheck
```

If your computer fulfills the requirement, you will see the following message.

```shell
The environment configuration is done!
```

### Section 1.3: Install cfDNAipe
Install cfDNApipe module.

```shell
pip install ./dist/cfDNApipe-0.0.4.tar.gz
```

## Section 2: WGBS Data Pipeline Demo
```Python
from cfDNApipe import *

Configure.setData('WGBS')
Configure.setThreads(20)
Configure.setGenome("hg19")
Configure.setRefDir(r'/home/wzhang/genome/hg19_bismark')
Configure.setOutDir(r'/data/wzhang/pipeline-test')
Configure.pipeFolderInit()

# Check references for pipeline, 'build = True' means download and build references which don't exist.
# If you don't want execute this step, just ignore this line or change build to Flase.
# This command is recommend for the first run, because it puts things right once and for all.
Configure.refCheck(build = True)

# Just put all your sequence data files in a folder
res1 = inputprocess(inputFolder = r"/data/wzhang/pipeline-test/raw-data")

# Quality Control
res2 = fastqc(upstream = res1)

# Identify adapters
res3 = identifyAdapter(upstream = res1, formerrun = res2)

# Remove adapters
res4 = adapterremoval(upstream = res3)

# Alignment using bismark
res5 = bismark(upstream = res4)

# Sort bam files
res6 = bamsort(upstream = res5)

# Remove duplicates, also you can use deduplicate_bismark, they are doing the same things.
res7 = rmduplicate(upstream = res6)

# Convert bam to 3 column bed files, this step will merge 2 reads to a single DNA fragment.
res8 = bam2bed(upstream = res7)

# Plot fragment length distribution
res9 = fraglenplot(upstream = res8)

# Compute methylation level for regions, default is CpG island from UCSC
res10 = computemethyl(upstream = res7, formerrun = res9)

# Group reads using picard
res11 = addRG(upstream = res7, formerrun = res10)
```
