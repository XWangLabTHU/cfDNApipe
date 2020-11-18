#!/usr/bin/perl -w

#Modified from VirusFinder2

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use IO::Handle();
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use Mosaic;

my @usage;
push @usage, "Usage: perl detect_integrationPlus.pl [options]\n\n";
push @usage, "Options:\n";
push @usage, "  --help          Usage information\n";
push @usage, "  --fq1           Sequencing reads are FASTQ file. This file is mandatory unless there is a file\n";
push @usage, "                    unmapped.1.fq in the directory specified by the user. The file unmapped.1.fq\n";
push @usage, "                    can be created from FASTQ files using preprocess.pl, a script for subtracting\n";
push @usage, "                    reads mapped to the host genome.\n";
push @usage, "  --fq2           Sequencing reads are FASTQ file. If this file is not provided and there isn't a\n";
push @usage, "                    file unmapped.2.fq in the directory specified by the user. The program will treat\n";
push @usage, "                    the input data as single end. Unmapped.2.fq can be created from FASTQ files\n";
push @usage, "                    using preprocess.pl, a script for subtracting reads mapped to the host genome.\n\n";
push @usage, "  --bamfile       The alignment bam file\n";
push @usage, "  --output        The directory to store software output\n";
push @usage, "  --virusDB       The virus database\n";
push @usage, "  --virus         Sequence file (in fasta format) of the virus. User can also provide virus's ID (or prefix\n";
push @usage, "                    of ID) in the virus database, e.g. gi_310698439_ref_NC_001526.2__Human_papillomavirus_type_16.\n";
push @usage, "                    This argument is mandatory unless user provides a file results-virus-top1.fa under\n";
push @usage, "                    the directory specified by the user.\n";
push @usage, "  --thread_no     Thread number, default is 1\n";
push @usage, "  --blastnIdxH    The index (prefix) used for blastn of human\n";
push @usage, "  --flank_reg     The flank region size used for sensitive mode of CREST, default is 4000\n";
push @usage, "  --mode          Sensitive or normal mode for virus integration detection, default value is normal\n";
push @usage, "  --SVDetect_dir  The directory where SVDetect located\n";

my $help;
#my $config_file;
my $output_dir;
my $read_file1;
my $read_file2;
my $virus_sequence;
my $mode = 'normal';
my $paired = 1;
my $bwa_bin = "bwa";
my $thread_no;
my $virus_database;
my $flank_region_size;
my $blastn_index_human;
my $alignment_file;
my $SVDetect_dir;

GetOptions
(
 'h|help|?'    => \$help,
# 'config=s'    => \$config_file,
 'output=s'    => \$output_dir,
 'virus=s'     => \$virus_sequence,
 'mode=s'      => \$mode,
 'fq1=s'       => \$read_file1,
 'fq2=s'       => \$read_file2,
 'bamfile=s'       => \$alignment_file,
 'thread_no=i' => \$thread_no,
 'virusDB=s'     => \$virus_database,
 'flank_reg=i' => \$flank_region_size,
 'blastnIdxH=s' => \$blastn_index_human,
 'SVDetect_dir=s' => \$SVDetect_dir,
);

$thread_no ||=1;
$flank_region_size ||=4000;

if ($help) {
   print @usage;
   exit(0);
}
#if (defined $config_file) {
#   if (!-e $config_file){
#   	print "\nThe configuration file $config_file does not exist!\n\n";
#        print @usage;
#	exit;
#   }
#}else{
#    print "Please provide a configuration file!\n\n";
#    print @usage;
#    exit;
#}

if (defined $output_dir) {
    if (!-e $output_dir){
    	print "\nThe output directory $output_dir does not exist!\n\n";
        print @usage;
	exit;
    }
}else{
    $output_dir = getcwd;
}
if (!-e "$output_dir/unmapped.1.fq" || !-s "$output_dir/unmapped.1.fq")
{
    if (defined $read_file1){
        if (!-e $read_file1 || !-s $read_file1){
            print "\nRead file1 $read_file1 does not exist!\n\n";
            print @usage;
            exit;
        }else{
            `rm "$output_dir/unmapped.1.fq"` if (-e "$output_dir/unmapped.1.fq");
            `ln -s $read_file1 $output_dir/unmapped.1.fq`;
        }
    }else{
    	print "\nPlease specify read file1 using '--fq1'!\n\n";
        print @usage;
	exit;
    }
}
if (!-e "$output_dir/unmapped.2.fq" || !-s "$output_dir/unmapped.2.fq")
{
    if (defined $read_file2){
        if (!-e $read_file2 || !-s $read_file2){
            print "\nRead file2 $read_file2 does not exist!\n\n";
            print @usage;
            exit;
        }else{
            `rm "$output_dir/unmapped.2.fq"` if (-e "$output_dir/unmapped.2.fq");
            `ln -s $read_file2 $output_dir/unmapped.2.fq`;
        }
    }else{
        $paired = 0;
    }
}



#my $config = new();
#$config->read($config_file);

#my $bwa_bin        = $config->get_value("bwa_bin");
my $script_dir     = "$FindBin::Bin";
#my $thread_no      = $config->get_value("thread_no");
#my $virus_database = $config->get_value("virus_database");


if (!-e "$output_dir/results-virus-top1.fa" || !-s "$output_dir/results-virus-top1.fa")
{
    if (defined $virus_sequence){
        if (!-e $virus_sequence || !-s $virus_sequence){

            print "$virus_sequence is not in the file system. Look up virus database...\n";
      	    PulloutTop1Virus($virus_sequence, $virus_database, "$output_dir/results-virus-top1.fa");

            if (!-e "$output_dir/results-virus-top1.fa" || !-s "$output_dir/results-virus-top1.fa"){
                print "\nError: $virus_sequence is not a file and it does not exist in virus database!\n\n";
                print @usage;
                exit;
            }
        }else{
            `rm "$output_dir/results-virus-top1.fa"` if (-e "$output_dir/results-virus-top1.fa");
            `ln -s $virus_sequence $output_dir/results-virus-top1.fa`;
        }
    }else{
    	print "\nPlease provide virus reference sequence using argument -v!\n\n";
        print @usage;
	exit;
    }
}
if(lc($mode) eq 'sensitive'){
    $mode = 'sensitive' ;
}else{
    $mode = 'normal' ;
}

#my $flank_region_size = 4000;
#if ($config->has_value("flank_region_size")) {
#    $flank_region_size = $config->get_value("flank_region_size");
#}

#if (scalar(@ARGV)<2){
#   print @usage;
#   exit(0);
#}
#
#my $config_file = $ARGV[0];
#my $output_dir  = $ARGV[1];
#
#
#die "\n$config_file does not exist!\n\n" if (!-e $config_file);
#die "\n$output_dir does not exist!\n\n" if (!-e $output_dir);
#die "\nCould not find unmapped.1.fq and unmapped.2.fq under the directory $output_dir!\n\n" if (!-e "$output_dir/unmapped.1.fq" || !-e "$output_dir/unmapped.2.fq");
#die "\nCould not find virus reference sequence, results-virus-top1.fa, under the directory $output_dir!\n\n" if (!-e "$output_dir/results-virus-top1.fa");
#
#my $mode = 'normal';
#if (scalar(@ARGV) > 2){
#   $mode = 'sensitive' if(lc($ARGV[2]) eq 'sensitive');
#}


my $start_dir = getcwd;
chdir($output_dir);


################################ Align unmapped reads to hg19+virus.fa #######################################

if (!-e "hg19.fa"){
#    my $blastn_index_human = $config->get_value("blastn_index_human");
    `ln -s $blastn_index_human.fa hg19.fa`;

}

$virus_sequence = "results-virus-top1.fa";

if (!-e 'hg19+virus.fa'){
    `cat hg19.fa $virus_sequence > hg19+virus.fa`;
}

if (!-e 'hg19+virus.fa.bwt'){
    system("$bwa_bin index hg19+virus.fa");
}

if (!-e 'alignment.sorted.bam'){

   my $RG = '@RG\tID:N\tSM:N\tLB:illumina\tPL:illumina';
   my $unmappedfq1 = 'unmapped.1.fq';
   my $unmappedfq2 = 'unmapped.2.fq';
   
    system("bwa mem -t $thread_no -M -R '$RG' hg19+virus.fa $unmappedfq1 $unmappedfq2 |samtools view -S -  -b -o alignment.bam");
    `samtools sort alignment.bam -o alignment.sorted.bam`;
    `samtools index  alignment.sorted.bam`;   

#   system("$bwa_bin aln -t $thread_no hg19+virus.fa $unmappedfq1 > unmapped.1.sai 2> /dev/null");
#
#   if ($paired){
#       system("$bwa_bin aln -t $thread_no hg19+virus.fa $unmappedfq2 > unmapped.2.sai 2> /dev/null");
#       system("$bwa_bin sampe -r '$RG' hg19+virus.fa unmapped.1.sai unmapped.2.sai $unmappedfq1 $unmappedfq2 > alignment.sam 2> /dev/null");
#   }else{
#       system("$bwa_bin samse -r '$RG' hg19+virus.fa unmapped.1.sai $unmappedfq1 > alignment.sam 2> /dev/null");
#   }
#
#   `samtools view -bS alignment.sam -o alignment.bam 2> /dev/null`;
#   `samtools sort   alignment.bam   alignment.sorted`;
#   `samtools index  alignment.sorted.bam`;
#   `rm alignment.sam alignment.bam`;
}


######################## Detect virus integration sites ###############################

my $ILIBs = GetIdir();
if ($paired){
    ## SVDetect
    if (!-e 'SVDetect' || !-e 'SVDetect/results-virus-loci.txt' || !-s 'SVDetect/results-virus-loci.txt'){
        print "Running SVDetect...\n";
        RunSVDetect();
    }
}

if (!-e 'crest' || !-e 'crest/results-virus-loci.txt' || !-s 'crest/results-virus-loci.txt'){
    print "Running CREST...\n";
    RunCREST();
}

if (!-e "$output_dir/crest/alignment.sorted.bam.cover" ||
    !-s "$output_dir/crest/alignment.sorted.bam.cover"){
    chdir $start_dir;
    die "\nWarning: CREST failed because there is no soft-clipped read!\n\n";
}

print "Consensus outputs...\n";
GetConsensusOutput();

`cp crest/results-virus-loci.txt results-virus-loci.txt`;

chdir $start_dir;


################################ Run SVDetect #######################################

sub RunSVDetect
{
    # Recuit aligned softcliped reads from the BAM file for SVDetect
#    if ($config->has_value("alignment_file") && !-e 'alignment+.bam') {
    if (defined $alignment_file) {
    	if (!-e 'alignment+.bam') {	
            my $RG = '@RG\tID:N\tSM:N\tLB:illumina\tPL:illumina';
            my $unmappedfq1 = 'unmapped.1+.fq';
            my $unmappedfq2 = 'unmapped.2+.fq';
            my $count = GetSoftClippedReads();

            if ($count > 0){
              `cat unmapped.1.fq >> unmapped.1+.fq`;
              `cat unmapped.2.fq >> unmapped.2+.fq`;

                system("bwa mem -t $thread_no -M -R '$RG' hg19+virus.fa $unmappedfq1 $unmappedfq2 |samtools view -S -  -b -o alignment+.bam");
#           system("$bwa_bin aln -t $thread_no hg19+virus.fa $unmappedfq1 > unmapped.1.sai 2> /dev/null");
#           system("$bwa_bin aln -t $thread_no hg19+virus.fa $unmappedfq2 > unmapped.2.sai 2> /dev/null");
#           system("$bwa_bin sampe -r '$RG' hg19+virus.fa unmapped.1.sai unmapped.2.sai $unmappedfq1 $unmappedfq2 > alignment+.sam 2> /dev/null");
#           `samtools view -bS alignment+.sam -o alignment+.bam`;
#           `rm alignment+.sam`;
            }
	    }
    }


    if (!-e 'SVDetect'){
       `mkdir SVDetect`;
    }

#    my $SVDetect_dir  = $config->get_value("SVDetect_dir");

    ## preprocess data
    if (!-e "SVDetect/alignment.sorted.ab.bam" && !-e 'SVDetect/alignment+.ab.bam'){
        if (-e 'alignment+.bam'){
            `perl $ILIBs $SVDetect_dir/scripts/BAM_preprocessingPairs.pl -t=1 -p=1 -o SVDetect alignment+.bam > SVDetect/preprocess.log 2>&1`;
        }else{
            `perl $ILIBs $SVDetect_dir/scripts/BAM_preprocessingPairs.pl -t=1 -p=1 -o SVDetect alignment.sorted.bam > SVDetect/preprocess.log 2>&1`;
        }
    }

    ## get lengths of reads
    my $read1_len = `sed '2q;d' unmapped.1.fq | awk '\{print length(\$0)\}'`;
    chomp $read1_len;
    my $read2_len = `sed '2q;d' unmapped.2.fq | awk '\{print length(\$0)\}'`;
    chomp $read2_len;

    # get lengths of chromosomes
    if (!-e 'SVDetect/hg19+virus.len'){
        my @chr_len = `samtools view -H alignment.sorted.bam | awk '/\@SQ/'`;
        open OUT, '>SVDetect/hg19+virus.len';
        my $idx = 1;
        foreach(@chr_len){
             s/[SN:|LN:]//g;

             print OUT $idx."\t".substr($_, 3);
             $idx++;
        }
        close OUT;
    }


    # write configuration file
    my $winsize;
    my $steplen;
    open IN, "<SVDetect/preprocess.log";
    while (<IN>){
        if (/mu length/){
            my @temp = split /,/;
            my $mu    = substr($temp[0], index($temp[0], '=')+1);
            my $sigma = substr($temp[1], index($temp[1], '=')+1);
            #$winsize  = 2*$mu+2.828*$sigma;
            #$steplen  = $winsize/4;
            $winsize  = $mu+2*$sigma;
            $steplen  = $winsize/4;
            last;
        }
    }
    close IN;

    open OUT, ">SVDetect/SVDetect.conf";
    print OUT "<general>\n";

    print OUT "input_format=bam\n";
    print OUT "sv_type=inter\n";
    print OUT "mates_orientation=FR\n";
    print OUT "read1_length=$read1_len\n";
    print OUT "read2_length=$read2_len\n";

    if (-e 'alignment+.bam'){
    	print OUT "mates_file=$output_dir/SVDetect/alignment+.ab.bam\n";
    }else{
    	print OUT "mates_file=$output_dir/SVDetect/alignment.sorted.ab.bam\n";
    }

    print OUT "cmap_file=$output_dir/SVDetect/hg19+virus.len\n";
    print OUT "num_threads=$thread_no\n";
    print OUT "output_dir=$output_dir/SVDetect\n";
    print OUT "tmp_dir=$output_dir/SVDetect\n";
    print OUT "</general>\n";

    print OUT "<detection>\n";
    print OUT "split_mate_file=1\n";
    print OUT "window_size=$winsize\n";
    print OUT "step_length=$steplen\n";
    print OUT "chromosomes=chrVirus\n";
    print OUT "</detection>\n";

    print OUT "<filtering>\n";
    print OUT "split_link_file=0\n";
    print OUT "nb_pairs_threshold=3\n";
    print OUT "strand_filtering=1\n";
    print OUT "</filtering>\n";
    close OUT;

    system("perl $ILIBs $SVDetect_dir/bin/SVDetect linking filtering links2SV -conf SVDetect/SVDetect.conf > /dev/null 2>&1");
    if (-e 'alignment+.bam'){
        if (!-e "$output_dir/SVDetect/alignment+.ab.bam.inter.links"){
            print "\nWarning: SVDetect did not work! Please make sure SVDetect is installed correctly!\n\n";
        }elsif(!-s "$output_dir/SVDetect/alignment+.ab.bam.inter.links"){
            print "\nWarning: SVDetect did find virus integration event!\n\n";
        }elsif(-e 'SVDetect/alignment+.ab.bam.inter.links.filtered.sv.txt'){
    	    `awk '{if(match(\$0,/chrVirus/) && \$13>=0.9)print}' SVDetect/alignment+.ab.bam.inter.links.filtered.sv.txt > SVDetect/results-virus-loci.txt`;
        }else{
            print "\nWarning: SVDetect did find virus integration event!\n\n";
        }
    }else{
        if (!-e "$output_dir/SVDetect/alignment.sorted.ab.bam.inter.links"){
            print "\nWarning: SVDetect did not work! Please make sure SVDetect is installed correctly!\n\n";
        }elsif(!-s "$output_dir/SVDetect/alignment.sorted.ab.bam.inter.links"){
            print "\nWarning: SVDetect did find virus integration event!\n\n";
        }elsif(-e 'SVDetect/alignment.sorted.ab.bam.inter.links.filtered.sv.txt'){
            `awk '{if(match(\$0,/chrVirus/) && \$13>=0.9)print}' SVDetect/alignment.sorted.ab.bam.inter.links.filtered.sv.txt > SVDetect/results-virus-loci.txt`;
        }else{
            print "\nWarning: SVDetect did find virus integration event!\n\n";
        }
    }
}

################## Get additional soft-clipped reads for SVDetect ######################

sub GetSoftClippedReads {
#    my $align_file = $config->get_value("alignment_file");
    my $align_file = $alignment_file;	
    my $awk_cmd    = '{if($6~/S/ && and($2,1292)==0)print $1}';
    my $sFILe      = "samtools view $align_file | awk '$awk_cmd' |";
    my %sReads     = ();

    open (IN, $sFILe);
    while (my $line = <IN>){
        chomp $line;
        $sReads{$line} = '';
    }
    close IN;

    my $count = 0;
    $sFILe = "samtools view $align_file |";
    open (IN, $sFILe);
    open (out1, '>unmapped.1+.fq');
    open (out2, '>unmapped.2+.fq');
    while (<IN>) {
    	my @data = split /\t/, $_;
        next if (!exists($sReads{$data[0]}));

        my $readend = ($data[1] >> 6) & 3;

        if (length($sReads{$data[0]})==0){
            $sReads{$data[0]} = "$readend\@$data[0]/$readend\n$data[9]\n+\n$data[10]\n";
        }else{
            my $prev_readend  = substr($sReads{$data[0]}, 0, 1);
            next if ($readend == $prev_readend);

            my @output = ();
            $output[$readend] = "\@$data[0]/$readend\n$data[9]\n+\n$data[10]\n";
            $output[$prev_readend] = substr($sReads{$data[0]}, 1);
            delete $sReads{$data[0]};

	    print out1 $output[1];
	    print out2 $output[2];
            $count++;
        }
    }
    close IN;
    close out1;
    close out2;
    return $count;
}


################################ Run CREST #######################################

sub RunCREST {

    my $CREST_dir      = "$script_dir/bin";
    my $cap3_bin       = "$CREST_dir/cap3";
    my $blat_bin       = "$CREST_dir/blat";
    my $faToTwoBit_bin = "$CREST_dir/faToTwoBit";

    if (!-e "crest"){
       `mkdir crest`;
    }

    ## get lengths of reads
    my $read1_len = `sed '2q;d' unmapped.1.fq | awk '\{print length(\$0)\}'`;
    chomp $read1_len;

    chdir($CREST_dir);

    ### Variant calling using CREST

    if (!-e "$output_dir/crest/alignment.sorted.bam.cover"){
       `perl extractSClip.pl -i $output_dir/alignment.sorted.bam -o $output_dir/crest/ --ref_genome $output_dir/hg19+virus.fa`;
    }
    if (!-e "$output_dir/crest/hg19+virus.2bit"){
       `$faToTwoBit_bin $output_dir/hg19+virus.fa $output_dir/crest/hg19+virus.2bit`;
    }
    if (!-e "$output_dir/crest/virus.2bit"){
       `$faToTwoBit_bin $output_dir/results-virus-top1.fa $output_dir/crest/virus.2bit`;
    }
    if (-e "$output_dir/crest/results.predSV.txt" && -s "$output_dir/crest/results.predSV.txt"){
        chdir $output_dir;
        return;
    }

    if ($paired){
    	if($mode eq 'normal'){
            #`perl $ILIBs CREST.pl -f $output_dir/crest/alignment.sorted.bam.cover -d $output_dir/alignment.sorted.bam --ref_genome $output_dir/hg19+virus.fa --boundary $output_dir/SVDetect/results-virus-loci.txt --flankregion 500                -t $output_dir/crest/hg19+virus.2bit -s $output_dir/crest/virus.2bit --cap3 $cap3_bin --blat $blat_bin -l $read1_len -o $output_dir/crest/ -p results --rmdup`;
            `perl $ILIBs CREST.pl -f $output_dir/crest/alignment.sorted.bam.cover -d $output_dir/alignment.sorted.bam --ref_genome $output_dir/hg19+virus.fa --boundary $output_dir/SVDetect/results-virus-loci.txt -t $output_dir/crest/hg19+virus.2bit -s $output_dir/crest/virus.2bit --cap3 $cap3_bin --blat $blat_bin -l $read1_len -o $output_dir/crest/ -p results --rmdup`;
        }else{
            #`perl CREST.pl -f $output_dir/crest/alignment.sorted.bam.cover -d $output_dir/alignment.sorted.bam --ref_genome $output_dir/hg19+virus.fa --boundary $output_dir/SVDetect/results-virus-loci.txt -t $output_dir/crest/hg19+virus.2bit -s $output_dir/crest/virus.2bit --cap3 $cap3_bin --blat $blat_bin -o $output_dir/crest/ -p results --rmdup --min_sclip_reads 1`;#  --sensitive
            `perl $ILIBs CREST.pl -f $output_dir/crest/alignment.sorted.bam.cover -d $output_dir/alignment.sorted.bam --ref_genome $output_dir/hg19+virus.fa --boundary $output_dir/SVDetect/results-virus-loci.txt --flankregion $flank_region_size -t $output_dir/crest/hg19+virus.2bit -s $output_dir/crest/virus.2bit --cap3 $cap3_bin --blat $blat_bin -l $read1_len -o $output_dir/crest/ -p results --rmdup  --sensitive`;
        }
    }else{
    	if($mode eq 'normal'){
            `perl $ILIBs CREST.pl -f $output_dir/crest/alignment.sorted.bam.cover -d $output_dir/alignment.sorted.bam --ref_genome $output_dir/hg19+virus.fa --nopaired -t $output_dir/crest/hg19+virus.2bit -s $output_dir/crest/virus.2bit --cap3 $cap3_bin --blat $blat_bin -l $read1_len -o $output_dir/crest/ -p results --rmdup`;
        }else{
            #`perl CREST.pl -f $output_dir/crest/alignment.sorted.bam.cover -d $output_dir/alignment.sorted.bam --ref_genome $output_dir/hg19+virus.fa --boundary $output_dir/SVDetect/results-virus-loci.txt -t $output_dir/crest/hg19+virus.2bit -s $output_dir/crest/virus.2bit --cap3 $cap3_bin --blat $blat_bin -o $output_dir/crest/ -p results --rmdup --min_sclip_reads 1`;#  --sensitive
            `perl $ILIBs CREST.pl -f $output_dir/crest/alignment.sorted.bam.cover -d $output_dir/alignment.sorted.bam --ref_genome $output_dir/hg19+virus.fa --nopaired -t $output_dir/crest/hg19+virus.2bit -s $output_dir/crest/virus.2bit --cap3 $cap3_bin --blat $blat_bin -l $read1_len -o $output_dir/crest/ -p results --rmdup  --sensitive`;
        }
    }

    chdir $output_dir;
}

################################ Get consensus output #######################################

sub GetConsensusOutput{

    ### Read SVDetect output

    my @SVDetectOutput = `awk \'{print \$4"\\t"\$5"\\t"\$3"\\t"\$7"\\t"\$8"\\t"\$6"\\t"\$9}\' SVDetect/results-virus-loci.txt`;
    my (@CHRs, @lowbound, @upbound, @span_num, @flag);

    foreach (@SVDetectOutput){
    	chomp;
        my @SD = split /\t/;

        my ($SD_chr, $SD_pos) = ($SD[0], $SD[1]);
        if ($SD_chr eq 'chrVirus'){
            $SD_chr = $SD[3];
            $SD_pos = $SD[4];
        }

        my @bound = split('-', $SD_pos);
        push @CHRs,     $SD_chr;
        push @lowbound, $bound[0];
        push @upbound,  $bound[1];
        push @span_num, $SD[6];
        push @flag,     0;
    }

    ### Output high confidence predictions of CREST

    open (OUT, '>crest/results-virus-loci.txt');

    my %outhash;
    if (-e 'crest/results.predSV.txt' && -s 'crest/results.predSV.txt'){
    	print OUT "Chromosome 1\tPosition 1\tStrand 1\tChromosome 2\tPosition 2\tStrand 2\t#Support reads (pair+softclip)\tConfidence\n";
        open (IN,  '<crest/results.predSV.txt');
        while (<IN>) {
            my @data = split /\t/;

            my ($data_chr, $data_pos, $split_n, $depth_n);
            if ($data[0] ne 'chrVirus'){
                $data_chr = $data[0];
                $data_pos = $data[1];
                $depth_n  = $data[9];
            }elsif ($data[4] ne 'chrVirus'){
                $data_chr = $data[4];
                $data_pos = $data[5];
            	$depth_n  = $data[10];
            }else{
            	next;
            }
            $split_n = $data[3]+$data[7];
            $data_chr = 'chr'.$data_chr if ($data_chr !~ /^chr/);

            next if (exists $outhash{"$data_chr$data_pos"});
            $outhash{"$data_chr$data_pos"} = 1;
            $outhash{$data_chr} = 1;

            my $span_n  = 0;
            for (my $i = 0; $i<scalar(@CHRs); $i++){
               if ($data_chr eq $CHRs[$i] && $data_pos>=$lowbound[$i] && $data_pos<=$upbound[$i]){
                  $span_n = $span_num[$i] if ($span_n < $span_num[$i]);
                  $flag[$i] = 1;
               }
            }
            $data[1] = $data[1]-1 if ($data[1] > 1);
            $data[5] = $data[5]-1 if ($data[5] > 1);
            print OUT "$data[0]\t$data[1]\t$data[2]\t$data[4]\t$data[5]\t$data[6]\t$span_n+$split_n\thigh\n";
        }
        close IN;
    }

    ### Output low confidence of CREST

    if (-e 'crest/results.predSV.txt-filtered' && -s 'crest/results.predSV.txt-filtered'){
    	if (!-e 'crest/results.predSV.txt' ||  !-s 'crest/results.predSV.txt'){
        print OUT "Chromosome 1\tPosition 1\tStrand 1\tChromosome 2\tPosition 2\tStrand 2\t#Support reads (pair+softclip)\tConfidence\n"; }

        open (IN,  '<crest/results.predSV.txt-filtered');
        while (<IN>){
            my @data = split /\t/;

            my ($data_chr, $data_pos, $split_n, $depth_n);
            if ($data[0] ne 'chrVirus'){
                $data_chr = $data[0];
                $data_pos = $data[1];
                $depth_n  = $data[9];
            }elsif ($data[4] ne 'chrVirus'){
                $data_chr = $data[4];
                $data_pos = $data[5];
            	$depth_n  = $data[10];
            }else{
            	next;
            }

            $split_n = $data[3]+$data[7];
            $data_chr = 'chr'.$data_chr if ($data_chr !~ /^chr/);

            next if (exists $outhash{"$data_chr$data_pos"});
            next if ($mode eq 'sensitive' && exists ($outhash{$data_chr}) && index($data_chr,'_')>=0);
            $outhash{"$data_chr$data_pos"} = 1;
            $outhash{$data_chr} = 1;

            my $span_n  = 0;
            my $toFilter = 0;
            for (my $i = 0; $i<scalar(@CHRs); $i++){
               if ($data_chr eq $CHRs[$i] && $data_pos>=$lowbound[$i] && $data_pos<=$upbound[$i]){
                  $span_n = $span_num[$i] if ($span_n < $span_num[$i]);
                  if ($flag[$i] == 0){
                  	$flag[$i] = 1;
                  }else{
                        $toFilter = 1;
                  }
               }
            }
            if ($toFilter == 0){
                $data[1] = $data[1]-1 if ($data[1] > 1);
                $data[5] = $data[5]-1 if ($data[5] > 1);
                print OUT "$data[0]\t$data[1]\t$data[2]\t$data[4]\t$data[5]\t$data[6]\t$span_n+$split_n\tlow\n";
            }
        }
        close IN;
    }

    ### Output consensus of CREST and SVDetect

    if (-e 'crest/results.predSV.txt2' && -s 'crest/results.predSV.txt2'){
        if ((!-e 'crest/results.predSV.txt' || !-s 'crest/results.predSV.txt') && (!-e 'crest/results.predSV.txt-filtered' || !-s 'crest/results.predSV.txt-filtered')){
            print OUT "Chromosome 1\tPosition 1\tStrand 1\tChromosome 2\tPosition 2\tStrand 2\t#Support reads (pair+softclip)\tConfidence\n";
        }
        open (IN,  '<crest/results.predSV.txt2');
        my %CrestOutput;
        while (<IN>){
            chomp;
            my @data = split /\t/;
	    $CrestOutput{$_} = $data[4];
        }
        close IN;

        for (my $i = 0; $i<scalar(@CHRs); $i++){
            next if ($flag[$i] == 1);
            my @SD = split /\t/, $SVDetectOutput[$i];

            my $prev_pos;
	    foreach my $key (sort {$CrestOutput{$b} <=> $CrestOutput{$a}} (keys(%CrestOutput))) {
                my @data = split /\t/, $key;
            	$data[2] = $data[2]-1 if ($data[2] > 1);

                $data[1] = 'chr'.$data[1] if ($data[1] !~ /^chr/);
                if ($data[1] eq $CHRs[$i] && $data[2]>=($lowbound[$i]-($mode eq 'sensitive'?$flank_region_size:500)) && $data[2]<=($upbound[$i]+($mode eq 'sensitive'?$flank_region_size:500))){
                     next if (exists $outhash{"$CHRs[$i]$data[2]"});
                     $outhash{"$CHRs[$i]$data[2]"} = 1;

                     if (!$prev_pos){
                         if ($SD[0] eq 'chrVirus'){
                             print OUT "$SD[0]\t$SD[1]\t\t$CHRs[$i]\t$data[2]\t$data[3]\t$SD[6]+$data[4]\tlow\n";
                         }else{
                             print OUT "$CHRs[$i]\t$data[2]\t$data[3]\t$SD[3]\t$SD[4]\t\t$SD[6]+$data[4]\tlow\n";
                         }
                     	 $prev_pos = $data[2];
                     }elsif(abs($prev_pos-$data[2])>10){
                         if ($SD[0] eq 'chrVirus'){
                             print OUT "$SD[0]\t$SD[1]\t\t$CHRs[$i]\t$data[2]\t$data[3]\t$SD[6]+$data[4]\tlow\n";
                         }else{
                             print OUT "$CHRs[$i]\t$data[2]\t$data[3]\t$SD[3]\t$SD[4]\t\t$SD[6]+$data[4]\tlow\n";
                         }
                         last;
                     }
                }
           }
        }
    }

    close OUT;

}