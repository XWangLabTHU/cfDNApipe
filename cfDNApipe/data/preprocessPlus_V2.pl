#!/usr/bin/perl -w

#modified from VirusFinder2


use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use Cwd;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin";
use Mosaic;

#my @usage;
#push @usage, "\nUsage:  preprocess.pl  <-c configuration file>  [options] \n\n";
#push @usage, "Options:\n";
#push @usage, "  -h, --help     Displays this information.\n";
#push @usage, "  -c, --config   Configuration file <required>.\n";
#push @usage, "  -o, --output   Full path of a directory to store results, default is current working directory.\n\n";
#push @usage, "Example:\n";
#push @usage, "  perl preprocess.pl -c config.txt -o /scratch/kingw/VirusFinder/simulation\n\n";
#push @usage, "Output summary:\n";
#push @usage, "  This program subtracts reads mapped to the host genome. It outputs 4 files for paired-end reads:\n";
#push @usage, "  (a) unmapped.1.fa, (b) unmapped.2.fa, (c) unmapped.1.fq, and (d) unmapped.2.fq. For single-end\n";
#push @usage, "  reads, it generates two files: (a) unmapped.1.fa, and (b) unmapped.1.fq\n\n";


my $help;
my $alignment_file;
my $output_dir;
my $fastq1;
my $fastq2;
my $paired = 1;
my $thread_no;

GetOptions
(
 'h|help|?'    => \$help,
 'bam=s'    => \$alignment_file,
 'output=s'    => \$output_dir,
 'fq1:s' => \$fastq1,
 'fq2:s' => \$fastq2,
 'thread_no=i' => \$thread_no,
 );

$thread_no ||=1;

my $usage = "Usage: $0 -bam <alignment_bam> -output <output_dir> -fq1 <unaligned_fastq1> -fq2 <unaligned_fastq2> -thread_no <thread_no>\n";
my $picard_bin = "picard SamToFastq ";

if ($help) {
   print $usage;
   exit(0);
}

if (!defined $alignment_file){
   print 'no bam';
   exit(0);
}

my $sample=basename($alignment_file,'.bam');

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

#if (defined $output_dir) {
#    if (!-e $output_dir){
#    	print "\nThe output directory $output_dir does not exist!\n\n";
#        print @usage;
#	exit;
#    }
#}else{
#    $output_dir = getcwd;
#}
#die "\nPlease provide:\n  (1) a configuration file, and \n  (2) a directory to store results.\n\n" if (scalar(@ARGV)<2);
#if (scalar(@ARGV)<2){
#  print "\nUsage:\t\tdetect_integration.pl <configuration file> <working directory>\n\n";
#  print "\t\tThis program requires: (1) a predefined configuration file; and, (2) full path of a directory to store results.\n\n";
#  print "\t\tThis program filters out reads mapped to the host genome. It outputs 4 files: (a) unmapped.1.fa, (b) unmapped.2.fa,\n";
#  print "\t\t(c) unmapped.1.fq, and, (d) unmapped.2.fq, which contains reads unaligned to the host genome.\n\n";
#  exit(0);
#}
#
#my $config_file = $ARGV[0];
#my $output_dir  = $ARGV[1];
#
#
#die "\n$config_file does not exist!\n\n" if (!-e $config_file);
#die "\n$output_dir does not exist!\n\n"  if (!-e $output_dir);
#

#my $config = new();
#$config->read($config_file);


#my $bowtie_index_human  = shift or die $usage;


#my $bowtie_bin         = $config->get_value("bowtie_bin");
#my $bowtie_index_human = $config->get_value("bowtie_index_human");
#my $thread_no          = $config->get_value("thread_no");



#if ($config->has_value("alignment_file")) {
#    PreprocessBam();
#}elsif ($config->has_value("fastq1")) {
#    PreprocessFastq();
#}else{
#    die "Please provide an alignment file (in BAM format) or one/two fastq files!\n";
#}

if (!-e "$output_dir/$sample.bam"){
    `ln -s $alignment_file $output_dir/$sample.bam`;
}

if (!-e "$output_dir/${sample}_unmapped.1.fa" || ($paired && !-e "$output_dir/${sample}_unmapped.2.fa")){

    print "Extract unaligned reads...\n";

    if ($paired){
        `samtools view -@ $thread_no -u -f 4  -F 264 $output_dir/${sample}.bam -o $output_dir/${sample}_L.unmapped.bam`;
        `samtools view -@ $thread_no -u -f 8  -F 260 $output_dir/${sample}.bam -o $output_dir/${sample}_R.unmapped.bam`;
        `samtools view -@ $thread_no -u -f 12 -F 256 $output_dir/${sample}.bam -o $output_dir/${sample}_B.unmapped.bam`;
        `samtools merge -@ $thread_no -u $output_dir/${sample}_LRB.unmapped.bam $output_dir/${sample}_[LRB].unmapped.bam`;
        `samtools sort -n -@ $thread_no $output_dir/${sample}_LRB.unmapped.bam -o $output_dir/${sample}_unmapped.bam`;

        `$picard_bin INCLUDE_NON_PF_READS=True VALIDATION_STRINGENCY=SILENT INPUT=$output_dir/${sample}_unmapped.bam FASTQ=$output_dir/${sample}_unmapped.1.fq SECOND_END_FASTQ=$output_dir/${sample}_unmapped.2.fq`;
        Fastq2Fasta("$output_dir/${sample}_unmapped.1.fq", "$output_dir/unmapped.1.fa");
        Fastq2Fasta("$output_dir/${sample}_unmapped.2.fq", "$output_dir/unmapped.2.fa");
    }else{
        `samtools view -@ $thread_no -uf 4 $output_dir/${sample}.bam > $output_dir/${sample}_unmapped.bam`;
        `$picard_bin INCLUDE_NON_PF_READS=True VALIDATION_STRINGENCY=SILENT INPUT=$output_dir/${sample}_unmapped.bam FASTQ=$output_dir/${sample}_unmapped.1.fq`;
        Fastq2Fasta("$output_dir/${sample}_unmapped.1.fq", "$output_dir/unmapped.1.fa");
    }
}

if (defined $fastq1){

	print "Extracting unaligned reads from bowtie2...";
	
	if($paired){
		Fastq2Fasta($fastq1, "$output_dir/${sample}_unmapped_tmp.1.fa");
		Fastq2Fasta($fastq2, "$output_dir/${sample}_unmapped_tmp.2.fa");
		`cat $output_dir/${sample}_unmapped_tmp.1.fa >> $output_dir/unmapped.1.fa`;
		`cat $output_dir/${sample}_unmapped_tmp.2.fa >> $output_dir/unmapped.2.fa`;
	}else{
		Fastq2Fasta($fastq1, "$output_dir/${sample}_unmapped_tmp.1.fa");
		`cat $output_dir/${sample}_unmapped_tmp.1.fa >> $output_dir/unmapped.1.fa`;
	}
}

#print "removing tmpFile...\n";
#`rm -rf $output_dir/*.bam $output_dir/*.fq $output_dir/*_tmp.*.fa`;
print "finish!\n";

#sub PreprocessBam {

#    my $alignment_file = $config->get_value("alignment_file");
#    die "Can't find BAM file $alignment_file!\n" if (!-e $alignment_file);
#    if (!-e "$output_dir/alignment.bam"){
#        `ln -s $alignment_file $output_dir/alignment.bam`;
#    }
#    $paired = `samtools view $alignment_file | head -10 | awk '{if(and(\$2,1)>0)print}' | wc -l`;
#    chomp $paired;
#}


#sub PreprocessFastq {
#
#    $fastq1 = $config->get_value("fastq1");
#    die "Can't find fastq file $fastq1!\n" if (!-e $fastq1);
#
#    if ($config->has_value("fastq2")){
#        $fastq2 = $config->get_value("fastq2");
#        die "Can't find fastq file $fastq2!\n" if (!-e $fastq2);
#    }else{
#        $paired = 0;
#    }
#
#    my $alignment_file = "$output_dir/alignment.bam";
#    return if (-e $alignment_file);
#
#    print "Do alignment using Bowtie2...\n";
#    if ($paired){
#       `$bowtie_bin -p $thread_no -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x $bowtie_index_human -1 $fastq1 -2 $fastq2 -S $output_dir/alignment.sam`;
#    }else{
#       `$bowtie_bin -p $thread_no -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x $bowtie_index_human -U $fastq1 -S $output_dir/alignment.sam`;
#    }
#
#    print "Convert SAM alignment file to BAM file...\n";
#    `samtools view -bS $output_dir/alignment.sam -o $alignment_file`;
#
#}
