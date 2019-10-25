# cfDNA-Pipeline
Cell Free DNA Sequencing Analysis Pipeline


# installation requirement
# The required environment is conda and python3.6
# you can create a new virtual environment, it will not influence your default invironment:
conda create -n cfDNApipe python=3.6
conda activate cfDNApipe

# enter cfDNA-Pipeline folder and run the following command to install the dependencies
chmod +x sysCheck
./sysCheck

#### WGBS data demo and parameters ####

# enter src folder and start python

# for WGS data

# set default parameters
Configure.Configure.setGenome("hg19")  # which genome to be used
Configure.Configure.setRefDir(r'/home/wzhang/genome/hg19_bismark') # genome reference path
Configure.Configure.setThreads(20) # how many thread to be used
Configure.Configure.setOutDir(r'/data/wzhang/pipeline-test') # output file path
Configure.Configure.pipeFolderInit() # setup pipeline

# put all your input files in a folder, then pass to inputprocess function
res1 = Fun_inputProcess.inputprocess(inputFolder = r"/data/wzhang/pipeline-test/raw-data")

# fastqc
# fastqInput: a list, contains all the input fastq files; if has a upstream input, please set None
# upstream: upstream output
# fastqcOutputDir: Output folder path
# other_params: a dict, contains other parameter. {"param1' : "True", "--param2' : True, "param3' : 20}
res2 = Fun_fastqc.fastqc(upstream = res1)

# identifyAdapter
# fqInput1, fqInput2: list of inputs
# upstream: upstream output
# formerrun: the former step, used to adjust folder number
# outputdir: output path
res3 = Fun_identifyAdapter.identifyAdapter(upstream = res1, formerrun = res2)

# adapterremoval
# fqInput1, fqInput2: list of inputs
# adapter1, adapter2: list of adapters
res4 = Fun_adapterremoval.adapterremoval(upstream = res3)

# bowtie2
# seqInput1, seqInput2: list of inputs
# ref: reference path
# other_params: {'-q': True, '-N': 1, '-X': 2000, '--no-mixed': True, 
#               '--no-discordant': True, '--dovetail': True, '--time': True}
res5 = Fun_bowtie2.bowtie2(upstream = res4)

# bamsort
# bamInput: list of inputs
res6 = Fun_bamsort.bamsort(upstream = res5)

# rmduplicate
# bamInput: list of inputs
res7 = Fun_rmDuplicate.rmduplicate(upstream = res6)

# bam2bed
# bamInput: list of inputs
res8 = Fun_bam2bed.bam2bed(upstream = res7)

# fraglenplot
# bedInput: list of inpu
res9 = Fun_fragLen.fraglenplot(upstream = res8)



