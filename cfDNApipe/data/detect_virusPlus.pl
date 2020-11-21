#!/usr/bin/perl -w

#Modified from VirusFinder2

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use Mosaic;
use threads;

my $usethreads = 1;

my @usage;
push @usage, "Usage: perl detect_virusPlus.pl [options]\n\n";
push @usage, "Options:\n";
push @usage, "  --help          Usage information\n";
push @usage, "  --fa1           The unmapped.1.fa read file\n";
push @usage, "  --fa2           The unmapped.2.fa read file\n";
push @usage, "  --output        The directory to store software output\n";
push @usage, "  --virusDB       The virus database\n";
push @usage, "  --blastnIdxH    The index (prefix) used for blastn of human\n";
push @usage, "  --blastnIdxV    The index (prefix) used for blastn of viral \n";
push @usage, "  --thread_no     Thread number, default is 1\n";
push @usage, "  --minContig     The minmum contig length used for Trinity, default is 300\n";
push @usage, "  --chopL         The read length chopped for fasta file before used in blat against virus database,\n";
push @usage, "                       default is 25\n";
push @usage, "  --blastnThrd    The threshold evalue used for blastn of human, default is 0.05\n";
push @usage, "  --simThrd       The minimum identity used for  blat against virus database, default is 0.8\n\n";
push @usage, "  --minIdent      The minimum identity used in blat against virus database,\n";
push @usage, "                       default is 80\n";
push @usage, "  --dustcutoff          The dust cutoff value used for discarding low complexity contigs.\n\n";

my $help;
#my $config_file;
my $output_dir;
my $read_file1;
my $read_file2;
#my $dust_cutoff = 7;
my $dust_cutoff;
my $paired = 1;

#others in config file
my $thread_no;
my $min_contig_length;
my $blastn_evalue_thrd;
my $similarity_thrd;
my $chop_read_length;
my $minIdentity;
my $blat_bin   = "blat";
my $virus_database;
#my $bowtie_bin  = "bowtie2";
#my $bowtie_index_human;
my $trinity_script  = "Trinity";
my $blastn_bin  ="blastn";
my $blastn_index_human;
my $blastn_index_virus;


GetOptions
(
 'h|help|?'      => \$help,
# 'config=s'      => \$config_file,
 'output=s'      => \$output_dir,
 'fa1=s'         => \$read_file1,
 'fa2=s'         => \$read_file2,
 'dustcutoff=s'  => \$dust_cutoff,
 'virusDB=s'     => \$virus_database,
# 'bowtieIdxH=s' => \$bowtie_index_human,
 'blastnIdxH=s' => \$blastn_index_human,
 'blastnIdxV=s' => \$blastn_index_virus,
 'thread_no=i' => \$thread_no,
 'minContig=i' => \$min_contig_length,
 'blastnThrd=f' => \$blastn_evalue_thrd,
 'simThrd=f' => \$similarity_thrd,
 'chopL=i' => \$chop_read_length,
 'minIdent=i' => \$minIdentity,
);


#others in config file
$thread_no ||=1;
$min_contig_length ||= 300;
$blastn_evalue_thrd ||= 0.05;
$similarity_thrd ||= 0.8;
$chop_read_length ||= 25;
$minIdentity ||= 80;
$dust_cutoff ||= 7;

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

if (!-e "$output_dir/unmapped.1.fa" || !-s "$output_dir/unmapped.1.fa")
{
    if (defined $read_file1){
        if (!-e $read_file1 || !-s $read_file1){
            print "\nRead file1 $read_file1 does not exist!\n\n";
            print @usage;
            exit;
        }else{
            `rm $output_dir/unmapped.1.fa` if (-e "$output_dir/unmapped.1.fa");
            `ln -s $read_file1 $output_dir/unmapped.1.fa`;
        }
    }else{
    	print "\nPlease specify read file1 using '--fa1'!\n\n";
        print @usage;
	exit;
    }
}

if (!-e "$output_dir/unmapped.2.fa" || !-s "$output_dir/unmapped.2.fa")
{
    if (defined $read_file2){
        if (!-e $read_file2 || !-s $read_file2){
            print "\nRead file2 $read_file2 does not exist!\n\n";
            print @usage;
            exit;
        }else{
            `rm $output_dir/unmapped.2.fa` if (-e "$output_dir/unmapped.2.fa");
            `ln -s $read_file2 $output_dir/unmapped.2.fa`;
        }
    }else{
    	$paired = 0;
    }
}
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
#
#
#die "\nCould not find unmapped.1.fa and unmapped.2.fa under the directory $output_dir!\n\n" if (!-e "$output_dir/unmapped.1.fa" || !-e "$output_dir/unmapped.2.fa");


#my $config = new();
#$config->read($config_file);


# tunable parameters
#my $thread_no          = $config->get_value("thread_no");

#my $min_contig_length  = $config->get_value("min_contig_length");
#my $blastn_evalue_thrd = $config->get_value("blastn_evalue_thrd");
#my $similarity_thrd    = $config->get_value("similarity_thrd");
#my $chop_read_length   = $config->get_value("chop_read_length");
#my $minIdentity        = $config->get_value("minIdentity");

# executables
#$my $script_dir         = "$FindBin::Bin";
#my $blat_bin           = "$script_dir/bin/blat";
#my $blat_bin          = $config->get_value("blat_bin");
#my $virus_database     = $config->get_value("virus_database");
#my $bowtie_bin         = $config->get_value("bowtie_bin");
#my $bowtie_index_human = $config->get_value("bowtie_index_human");
#my $trinity_script     = $config->get_value("trinity_script");
#my $blastn_bin         = $config->get_value("blastn_bin");
#my $blastn_index_human = $config->get_value("blastn_index_human");
#my $blastn_index_virus = $config->get_value("blastn_index_virus");


my $start_dir = getcwd();
chdir($output_dir);


############################################ map reads to virus database #############################################


print "step 2.1 map reads to virus database...\n";

if (!-e 'chopped_unmapped.1.psl' || ($paired && !-e 'chopped_unmapped.2.psl')){
    if ($usethreads && $thread_no > 1){
	print "using 2 threads to map reads...\n";
        my @threads;
#        my $t1 = threads->new(sub{ blat_thread1($virus_database, $minIdentity, $chop_read_length, 'unmapped.1.fa') });
		my $t1 = threads->new(sub{ blat_thread1($virus_database, $minIdentity, 'unmapped.1.fa') });
        push(@threads,$t1);
        if ($paired){
#            my $t2 = threads->new(sub{ blat_thread1($virus_database, $minIdentity, $chop_read_length, 'unmapped.2.fa') });
			my $t2 = threads->new(sub{ blat_thread1($virus_database, $minIdentity, 'unmapped.2.fa') });
            push(@threads,$t2);
        }

        foreach (@threads) {
            $_->join;
        }
    }else{
#        blat_thread1($virus_database, $minIdentity, $chop_read_length, 'unmapped.1.fa');
		blat_thread1($virus_database, $minIdentity, 'unmapped.1.fa');
        if ($paired){
#            blat_thread1($virus_database, $minIdentity, $chop_read_length, 'unmapped.2.fa');
			blat_thread1($virus_database, $minIdentity, 'unmapped.2.fa');
        }
    }

    print "Find reads mapped to virus database...\n";
}

print "step 2.2 garner reads mapped to the virus database...\n";
if (!-e 'blat_out_candidate_singlelane.fa'){
   `cp chopped_unmapped.1.psl chopped_singlelane.psl`;
   `cp unmapped.1.fa singlelane_unmapped.fa`;

   if ($paired){
       `awk 'NR>5' chopped_unmapped.2.psl >> chopped_singlelane.psl`;
       `cat unmapped.2.fa >> singlelane_unmapped.fa`;
   }
   BlatOutCandidate('chopped_singlelane.psl',  'singlelane_unmapped.fa', 'blat_out_candidate_singlelane.fa');
#   `rm  chopped_singlelane.psl singlelane_unmapped.fa`;
}

############################################ de novo assembly #############################################


print "step 2.3 de novo assembly using Trinity...\n";

my $updated = 0;

if (!-e 'Trinity-init.fasta'){
    $updated = 1;

    my $ILIBs = GetIdir();
#    `perl $ILIBs $trinity_script --seqType fa  --JM 4G --single blat_out_candidate_singlelane.fa --min_contig_length $min_contig_length  --output trinity_output --CPU $thread_no --bfly_opts \"-V 10 --stderr\"`;
    `$trinity_script --seqType fa  --max_memory 10G --single blat_out_candidate_singlelane.fa --min_contig_length $min_contig_length  --output trinity_output --CPU $thread_no --verbose`;

    if (-e 'trinity_output/Trinity.fasta' && -s 'trinity_output/Trinity.fasta'){
    	print "Finished read assembly\n";
    }else{
        if (!-e 'trinity_output'){
            die "\nError: Trinity did not output contigs. Please make sure Trinity works.\n\n";
        }else{
            `rm -r trinity_output`;

            print "\nWarning: Trinity did not output contigs (min contig length=$min_contig_length)!\n\n";

            my $contig_len = $min_contig_length >> 1;
            my $read1_len  = `sed '2q;d' unmapped.1.fa | awk '\{print length(\$0)\}'`;
            chomp $read1_len;
            $contig_len = ($min_contig_length + $read1_len) >> 1 if ($contig_len <= $read1_len);

            print "\nRetry Trinity using a shorter contig length $contig_len...\n\n";

#            `$trinity_script --seqType fa  --JM 4G --single blat_out_candidate_singlelane.fa --min_contig_length $contig_len  --output trinity_output --CPU $thread_no --bfly_opts \"-V 10 --stderr\"`;
            `$trinity_script --seqType fa  --max_memory 10G --single blat_out_candidate_singlelane.fa --min_contig_length $contig_len  --output trinity_output --CPU $thread_no --verbose"`;

	    if (!-e 'trinity_output/Trinity.fasta' || !-s 'trinity_output/Trinity.fasta'){
            	`rm -r trinity_output`;
            	die "\nError: Trinity did not output contigs!\n\n";
            }else{
    		print "Finished read assembly\n";
            }
        }
    }

    `mv trinity_output/Trinity.fasta Trinity-init.fasta`;
#    `rm -r trinity_output`;
}

if ($updated || !-e 'Trinity.fasta'){
    $updated = 1;
    print "Discard low complexity contigs...\n";
    ReformatOutput('Trinity-init.fasta', 'Trinity.fasta', $dust_cutoff);
}

if ($updated || !-e 'clean_blastn.fa'){
    $updated = 1;
    print "blastn contigs against human genome...\n";
#    `$blastn_bin -query=Trinity.fasta -db=$blastn_index_human -evalue $blastn_evalue_thrd -outfmt 6 > human_contig.txt`;
    `$blastn_bin -query Trinity.fasta -db $blastn_index_human -evalue $blastn_evalue_thrd -outfmt 6 > human_contig.txt`;
	
    print "clean up blastn outputs\n";
    BlastnCleanup('human_contig.txt', 'Trinity.fasta', 'clean_blastn.fa', $similarity_thrd);   #pull out non-human contig into clean_blastn.fa
}

############################################ Pull out virus sequences #############################################

print "step 2.4 detect virus sequences...\n";
print "blastn trinity output against candidate viruses\n";
`cp clean_blastn.fa non_human_contig.fa`;

if ($updated || !-e 'singlelane.psl' || !-e 'non_human_contig_blastn.txt'){
    $updated = 1;
    `$blastn_bin -query non_human_contig.fa -db $blastn_index_virus -evalue $blastn_evalue_thrd -outfmt 6 > non_human_contig_blastn.txt`;


    print "step 2.5 map reads to contigs...\n\n";

    #`$blat_bin non_human_contig.fa -minIdentity=98 blat_out_candidate_singlelane.fa singlelane.psl`;

    if ($usethreads && $thread_no > 1){
        # split file
        my $line_num = `wc -l blat_out_candidate_singlelane.fa`;
        $line_num = substr($line_num, 0, length($line_num)-length('blat_out_candidate_singlelane.fa')-1);
        $line_num = (($line_num>>2) << 1) + 2;
        `split -l $line_num blat_out_candidate_singlelane.fa blat_out`;

        # mapping
        my @threads;
        my $t1 = threads->new(sub{ blat_thread2('non_human_contig.fa', 98, 'blat_outaa') });
        push(@threads,$t1);
        my $t2 = threads->new(sub{ blat_thread2('non_human_contig.fa', 98, 'blat_outab') });
        push(@threads,$t2);

        foreach (@threads) {
            $_->join;
        }

        `cp blat_outaa.psl singlelane.psl`;
        `awk 'NR>5' blat_outab.psl >> singlelane.psl`;
         `rm  blat_outa*`;
    }else{
        `$blat_bin non_human_contig.fa -minIdentity=98 blat_out_candidate_singlelane.fa singlelane.psl`;
    }
}


print "step 2.6 Pull out virus sequence...\n";

if ($updated || !-e 'results-virus.txt' || !-s 'results-virus.txt'){
   PulloutVirus('non_human_contig.fa', 'non_human_contig_blastn.txt', 'singlelane.psl');
}

if ($updated || !-e 'results-virus-top1.fa' || !-s 'results-virus-top1.fa'){
   PulloutTop1Virus('results-virus.txt', $virus_database, 'results-virus-top1.fa');
}

if ($updated || !-e 'results-virus-list.txt' || !-s 'results-virus-list.txt'){
   PulloutVirusList('non_human_contig.fa', 'non_human_contig_blastn.txt', 'singlelane.psl','human_contig.txt');
}

chdir $start_dir;

if (!-e "$output_dir/results-virus.txt" ||
    !-e "$output_dir/results-virus-top1.fa" ||
    !-s "$output_dir/results-virus-top1.fa"){
    print "\nFailed to detect viruses. \n\n";
    system "date";
    exit;
}

`ln -s $output_dir/results-virus.txt           $output_dir/virus.txt`;
`cp $output_dir/results-virus-list.txt      $output_dir/virus-list.txt`;
`ln -s $output_dir/results-contig.txt          $output_dir/contig.txt`;
`ln -s $output_dir/results-novel-contig-2.fa   $output_dir/novel-contig.fa` if (-e '$output_dir/results-novel-contig-2.fa');

if (-s "$output_dir/results-virus-top1.fa"){
    my $ret = `awk '{if(NR==2)print\$1}' $output_dir/virus.txt`;
    print "\nVirus detected, top one: $ret\n\n";
    print "Look at virus.txt, virus-list.txt, contig.txt and novel-contig.fa for details.\n\n\n";
}

sub blat_thread1 {
   my $blatRef = shift;
   my $cutoff  = shift;
#   my $minScore =shift;  #add
   my $faFile  = shift;

   my $chopped_faFile = 'chopped_'.$faFile;
   my $pslFile = $chopped_faFile;
   $pslFile =~ s/fa/psl/;

   if (!-e $chopped_faFile){
   	ChopReads($faFile, $chopped_faFile, $chop_read_length);
   }
#   `$blat_bin $blatRef -minIdentity=$cutoff -minScore=$minScore $chopped_faFile $pslFile`;
   `$blat_bin $blatRef -minIdentity=$cutoff -minScore=22 $chopped_faFile $pslFile`;
}


sub blat_thread2 {
   my $blatRef = shift;
   my $cutoff  = shift;
   my $faFile  = shift;

   my $pslFile = $faFile.'.psl';

   `$blat_bin $blatRef -minIdentity=$cutoff $faFile $pslFile`;
}
