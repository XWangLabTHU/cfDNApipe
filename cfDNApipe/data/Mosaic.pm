#!/usr/bin/perl -w

#from VirusFinder2

package Mosaic;

use diagnostics;
use strict;
use warnings;
use Bio::SeqIO;

$| = 1; #disable buffer

use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
our $VERSION     = 1.0;
our @ISA         = qw (Exporter);
our @EXPORT      = qw ($DEBUG_MODE new read has_value get_value average stdev Fastq2Fasta ChopReads ReformatOutput CalDustScore BlatOutCandidate BlastnCleanup PulloutVirus PulloutTop1Virus PulloutVirusList GetIdir GetChrPrefix GetInsertSize GetRefSubSeq FilterRefSubSeq WriteRefSubSeq ExtractMappedReads);
our @EXPORT_OK   = qw ($DEBUG_MODE new read has_value get_value average stdev Fastq2Fasta ChopReads ReformatOutput CalDustScore BlatOutCandidate BlastnCleanup PulloutVirus PulloutTop1Virus PulloutVirusList GetIdir GetChrPrefix GetInsertSize GetRefSubSeq FilterRefSubSeq WriteRefSubSeq ExtractMappedReads);
our %EXPORT_TAGS = (DEFAULT => [qw(&new &read &has_value &get_value &average &stdev &Fastq2Fasta &ChopReads &ReformatOutput &CalDustScore &BlatOutCandidate &BlastnCleanup &PulloutVirus &PulloutTop1Virus &PulloutVirusList &GetIdir &GetChrPrefix &GetInsertSize &GetRefSubSeq &FilterRefSubSeq &WriteRefSubSeq &ExtractMappedReads)]);



sub new {
  my $self = {};
  bless($self);
  return $self;
}


sub read {
  my $self = shift;
  my $config_file = shift;
  my %config_values;

  open CFG, $config_file or die "Error: Unable to open $config_file\n";
  while (<CFG>){
      chomp;
      /^\s*([^=\s]+)\s*=\s*(.*)$/;

      my $key = $1;
      my $value = $2;

      next if not defined $key;
      next if not defined $value;

      $config_values{$key} = $value;
  }
  close CFG;

  foreach my $key (keys %config_values){
      while ($config_values{$key} =~ /\$\(([^)]+)\)/){
	  my $other_key = $1;

	  if (not defined $config_values{$other_key}){
	      die "Error: no value for $other_key in config file $config_file\n";
	  }

	  $config_values{$key} =~ s/\$\($other_key\)/$config_values{$other_key}/;
      }
  }

  $self->{"config_values"} = \%config_values;
  $self->{"config_file"} = $config_file;
}


sub has_value {
  my $self = shift;
  my $key = shift;

  my $config_values = $self->{"config_values"};
  my $config_file = $self->{"config_file"};

  defined $config_values and defined $config_file or die "Error: config not read\n";

  return defined $config_values->{$key};
}


sub get_value {
  my $self = shift;
  my $key = shift;

  my $config_values = $self->{"config_values"};
  my $config_file = $self->{"config_file"};

  defined $config_values and defined $config_file or die "Error: config not read\n";

  return $config_values->{$key};
}

sub average{
    my($data) = @_;
    if (not @$data) {
       	die("Empty array\n");
    }
    my $total = 0;
    foreach (@$data) {
    	$total += $_;
    }
    my $average = $total/@$data;
    return $average;
}

sub stdev{
    my($data) = @_;
    if(@$data == 1){
    	return 0;
    }
    my $average = &average($data);
    my $sqtotal = 0;
    foreach(@$data) {
    	$sqtotal += ($average-$_) ** 2;
    }
    my $std = ($sqtotal/(@$data-1)) ** 0.5;
    return $std;
}

sub Fastq2Fasta {
    my $infile  = shift @_;
    my $outfile = shift @_;

    my $InFileStr;
    if ($infile =~ /\.gz$/){
        $InFileStr = "gunzip -c $infile |";
    }else{
        $InFileStr = "<$infile";
    }

    open (IN, $InFileStr);
    open (OUT, ">$outfile");
    my $count = 0;
    while (my $line=<IN>) {
        $count ++;
        if ($count%1000000 eq 0) {
            print "$count\n";
        }

        chomp $line;

        $line=~s/\s+/\_/g;
        print OUT ">$line\n";
        $line = <IN>;
        print OUT $line;
        $line = <IN>;
        $line = <IN>;
    }

    close IN;
    close OUT;
}


sub ChopReads {
    my $fasta       = shift @_;
    my $out         = shift @_;
    my $chop_length = shift @_;

    open (FASTA, "<$fasta");
    open (OUT, ">$out");

    my ($count, %namesfasta, $title);

    while (my $line = <FASTA>) {
        chomp $line;
        $count++;
        if ($count%1000000 == 0){
            print "$count\n";
        }
        if ($count%2 == 1) {
            $title = $line;
        }else{
            my @sequence = split (//, $line);
            my $totalbases = $#sequence + 1;
            my $number = int($totalbases/$chop_length);

            my $i;
            for ($i=1; $i<=$number; $i++){
                my $start = ($i-1)*$chop_length;
                my $end = $i*$chop_length - 1;
                my $seq =join ("", @sequence[$start..$end]);
                print OUT "$title$i\n$seq\n";
            }

            if (($totalbases - $number*$chop_length)>=20) {
                my $start = $#sequence - $chop_length + 1;
                my $end = $#sequence;
                my $seq = join ("", @sequence[$start..$end]);
                print OUT "$title$i\n$seq\n";
            }
        }
    }

    close (FASTA);
    close (OUT);
}

sub ReformatOutput {
    my $infile         = shift @_;
    my $outfile        = shift @_;
    my $low_com_cutoff = shift @_;

    die "" if (-z $infile);
    open (IN,  "<$infile");
    open (OUT, ">$outfile");

    my ($line, $name, $seq);

    while ($line=<IN>) {
        chomp $line;

        if ($line=~/^\>/) {
             if ($name){
                if (CalDustScore($seq) <= $low_com_cutoff){
                   print OUT "$name\n$seq\n";
                }
             }
             $name = $line;
             $seq  = '';
        }else {
             if ($seq){
                 $seq .= "$line";
             }else{
                 $seq  = $line;
             }
        }
    }
    if (CalDustScore($seq) <= $low_com_cutoff){
    	print OUT "$name\n$seq\n";
    }

    close IN;
    close OUT;
}

sub CalDustScore{
    #Check for sequence complexity. Adapted from PRINSEQ
    my $seqn  = shift @_;
    my $length = length($seqn);

    my $WINDOWSIZE = 64;
    my $WINDOWSTEP = 32;
    my $WORDSIZE = 3;
    my @WINDOWSIZEARRAY = (0..61);
    my $POINTFIVE = 1/2;

    my ($rest,$steps,@vals,$str,$num,$bynum,%counts);
    if($length <= $WINDOWSIZE) {
        $rest = $length;
        $steps = 0;
    } else {
        $steps = int(($length - $WINDOWSIZE) / $WINDOWSTEP) + 1;
        $rest = $length - $steps * $WINDOWSTEP;
        unless($rest > $WINDOWSTEP) {
            $rest += $WINDOWSTEP;
            $steps--;
        }
    }
    $num = $WINDOWSIZE-2;
    $bynum = 1/$num;
    $num--;

    my $dustscore;
    foreach my $i (0..$steps-1) {
        $str = substr($seqn,($i * $WINDOWSTEP),$WINDOWSIZE);
        %counts = ();
        foreach my $i (@WINDOWSIZEARRAY) {
            $counts{substr($str,$i,3)}++;
        }
        $dustscore = 0;
        foreach(values %counts) {
            $dustscore += ($_ * ($_ - 1) * $POINTFIVE);
        }
        push(@vals,($dustscore * $bynum));
    }
    #last step
    if($rest > 5) {
        $str = substr($seqn,($steps * $WINDOWSTEP),$rest);
        %counts = ();
        $num = $rest-2;
        foreach my $i (0..($num - 1)) {
            $counts{substr($str,$i,3)}++;
        }
        $dustscore = 0;
        foreach(values %counts) {
            $dustscore += ($_ * ($_ - 1) * $POINTFIVE);
        }
        push(@vals,(($dustscore / ($num-1)) * (($WINDOWSIZE - 2) / $num)));
    } else {
        push(@vals,31); #to assign a maximum score based on the scaling factor 100/31
    }
    my $mean = 0;
    $mean += $_ for @vals;
    $mean /= @vals;

    return int($mean * 100 / 31);
}

sub BlatOutCandidate {
    my $infile1 = shift @_;
    my $infile2 = shift @_;
    my $outfile = shift @_;

    my %flag;
    my $chopname = 0;
    if ($infile1 =~ /chopped/){
        $chopname = 1;
    }

    open (IN, "<$infile1");
    my $line;
    for (my $i=1;$i<=5;$i++) {
      $line = <IN>;
    }

    while ($line=<IN>) {
        chomp $line;
        my @data = split /\t/, $line;
        my $name = $data[9];
        if ($chopname){
             chop $name;
        }
        $name =~s/\/[0-9]+$//;
        $flag{$name} = 1;
    }
    close IN;


    open (IN, "<$infile2");
    open (OUT, ">$outfile");

    while (my $line=<IN>) {
        chomp $line;
        my $name = $line;
        $name =~s/\/[0-9]+$//;
        $name =~s/^\>//;

        if ($flag{$name}) {
            print OUT "$line\n";
            $line = <IN>;
            print OUT $line;
        } else {
            $line = <IN>;
        }
    }

    close IN;
    close OUT;
}


sub BlastnCleanup {
    my $input       = shift @_;
    my $fasta_input = shift @_;
    my $output      = shift @_;
    my $similarity_thrd = shift @_;

    my %flag = ();
    my %seqlen = ();
    my $line;

    open (IN, "<$fasta_input");
    while ($line=<IN>){
      my @temp = split (/ /, $line);
      my $length = substr($temp[1], 4);
      $seqlen{substr($temp[0],1)} = $length;
      $line = <IN>;
    }
    close IN;

    open (IN, "<$input");
    while (<IN>){
      chomp;
      my @data = split /\t/;
      my $numbermatch = $data[3];
      my $evalue = $data[10];
      my $bit_score = $data[-1];

      die "Wrong data format\n!" if (!($seqlen{$data[0]}));

      my $length = $seqlen{$data[0]};
      my $percent = $bit_score/$length;

      if ($percent >= $similarity_thrd){
          $flag{$data[0]} = 1;
      }
    }

    close IN;

    open (IN, "<$fasta_input");
    open (OUT, ">$output");
    while ($line=<IN>) {
      chomp $line;
      my @temp = split (/ /, $line);

      if (!exists $flag{substr($temp[0],1)}) {
        print OUT "$line\n";
        $line = <IN>;
        print OUT $line;
      }
      else {
        $line = <IN>;
      }
    }

    close IN;
    close OUT;
}

sub PulloutVirus {

    my $contig_fa_file     = shift @_;
    my $contig_blastn_file = shift @_;
    my $psl_file           = shift @_;

    my $line;
    my $name;
    my %seq = ();
    my %count = ();
    my %seqlen = ();
    my %contig_virus  = ();
    my %contig_virus2 = ();
    #my %virus_read = ();

    open (IN, "<$contig_fa_file");
    while ($line=<IN>) {
      chomp $line;

      my @temp = split (/ /, $line);
      my $length = substr($temp[1], 4);
      $seqlen{substr($temp[0], 1)} = $length;

      $line = <IN>;
      chomp $line;
      $seq{substr($temp[0], 1)} = $line;
    }
    close IN;


    my @contig_virus_lines;
    open (IN, "<$contig_blastn_file");

    while ($line=<IN>) {
      chomp $line;
      my @data   = split /\t/, $line;
      my $length = $seqlen{$data[0]};

      next if ($data[11]/$length < 0.75 || $data[2] < 85);

      push @contig_virus_lines, $line;

      if (!(exists $contig_virus{$data[0]}{$data[1]}) || $contig_virus{$data[0]}{$data[1]} < $data[3]){
          $contig_virus{$data[0]}{$data[1]}  = $data[3];
          $contig_virus2{$data[0]}{$data[1]} = $data[2];
      }
    }
    close IN;


    open (IN, "<$psl_file");
    for (my $i=1;$i<=5;$i++) {
      $line = <IN>;
    }
    while ($line=<IN>) {
      chomp $line;
      my @data = split /\t/, $line;
      next if (scalar(@data) < 14);

      $count{$data[13]} ++;

      #for my $v (keys %{$contig_virus{$data[13]}} ) {
      #    $virus_read{$v}{$data[9]} = 1;
      #}
    }
    close IN;

    my $output_prefix = 'results';
    open (OUT,  ">$output_prefix-contig.txt");
    open (OUT2, ">$output_prefix-contig.fa");
    print OUT "Contig name\t#Reads\tVirus\tE-value\tBit score\tSequence\n";

    my %flag = ();
    my $max_depth = 0;
    foreach $line (@contig_virus_lines) {
      my @data = split /\t/, $line;

      my $pline = "$data[0]\t$count{$data[0]}\t$data[1]\t$data[10]\t$data[11]\t$seq{$data[0]}";

      if (!exists($flag{$pline})) {
          print OUT  "$pline\n";
          $flag{$pline} = 1;
      }
      if (!exists($flag{$data[0]})) {
          print OUT2 "\>$data[0]\n$seq{$data[0]}\n";
          $flag{$data[0]} = 1;
      }
      if (exists $count{$data[0]}){
          my $depth = $count{$data[0]}/$seqlen{$data[0]};
          if ($depth > $max_depth){
             $max_depth = $depth;
          }
      }
    }

    close OUT;
    close OUT2;

    my %virus_contig = ();
    my %virus_contig_len = ();
    my %virus_contig_len_max = ();
    my %virus_contig_max_len = ();
    my %contig_read_num = ();
    my %virus_read_depth = ();

    for my $c (keys %contig_virus) {
        for my $v (keys %{$contig_virus{$c}} ) {
            if (exists $virus_contig{$v}){
                $virus_contig{$v} = $virus_contig{$v}.', '.$c;
                $virus_contig_len{$v} = $virus_contig_len{$v}.', '.$contig_virus{$c}{$v}.'/'.$contig_virus2{$c}{$v};
                if ($virus_contig_len_max{$v} < $contig_virus{$c}{$v}*$contig_virus2{$c}{$v}){
                    $virus_contig_len_max{$v} = $contig_virus{$c}{$v}*$contig_virus2{$c}{$v};
                }
                $virus_contig_max_len{$v} = $virus_contig_max_len{$v}.', '.$seqlen{$c};
                $contig_read_num{$v} = $contig_read_num{$v}.', '.$count{$c};
            }else{
                $virus_contig{$v} = $c;
                $virus_contig_len{$v} = $contig_virus{$c}{$v}.'/'.$contig_virus2{$c}{$v};
                $virus_contig_len_max{$v} = $contig_virus{$c}{$v}*$contig_virus2{$c}{$v};
                $virus_contig_max_len{$v} = $seqlen{$c};
                $contig_read_num{$v} = $count{$c};
            }
	    if (exists $virus_read_depth{$v}){
            	$virus_read_depth{$v} = $count{$c}/$seqlen{$c} if ($virus_read_depth{$v} < $count{$c}/$seqlen{$c});
            }else{
            	$virus_read_depth{$v} = $count{$c}/$seqlen{$c};
            }
        }
    }

    my $top1 = '';
    open (OUT, ">$output_prefix-virus.txt");
    print OUT "Virus name\tContigs\tContig length (bp)\tMapped length/rate of contigs\t#Reads fallen on contigs\n";
    for my $v (sort {$virus_contig_len_max{$b}<=>$virus_contig_len_max{$a}} keys %virus_contig_len_max) {
    	#next if ($virus_read_depth{$v} < 0.2*$max_depth);
        print OUT "$v\t$virus_contig{$v}\t$virus_contig_max_len{$v}\t$virus_contig_len{$v}\t$contig_read_num{$v}\n";
        if (length($top1) == 0){
            $top1 = $v;
        }
    }
    close OUT;

}

sub PulloutTop1Virus {
    my $virus_id_file  = shift @_;
    my $blat_reference = shift @_;
    my $out_file_name  = shift @_;

    my $top1;
    if (!-e $virus_id_file){
        $top1 = $virus_id_file;
    }else{
        open (IN, $virus_id_file);
        my $line=<IN>;
        $line=<IN>;
        if (!$line){
            `touch $out_file_name`;
            return;
        }

        my @data = split /\t/, $line;
        $top1 = $data[0];
        close IN;
    }

    open (IN, $blat_reference);
    open (OUT, '>'.$out_file_name.'-temp');

    my $flag = 0;
    while (my $line=<IN>) {
      if ($flag == 1){
          if (substr($line,0,1) eq '>'){
              last;
          }else{
              print OUT $line;
          }
      }else{
          if ($line =~ /^\>\Q$top1\E/){
              $flag = 1;
              print OUT ">chrVirus\n";
          }
      }
    }
    close IN;
    close OUT;

    my $in  = Bio::SeqIO->new(-file => $out_file_name.'-temp', -format => 'Fasta');
    my $out = Bio::SeqIO->new(-file => '>'.$out_file_name,     -format => 'Fasta');
    while ( my $seq = $in->next_seq() ) {$out->write_seq($seq); }
    `rm $out_file_name-temp`;
}


sub PulloutVirusList {

    my $contig_fa_file     = shift @_;
    my $contig_blastn_file = shift @_;
    my $psl_file           = shift @_;
    my $human_contigs      = shift @_;

    my (%contig_virus, %contig_len, %contig_identity, %contig_evalue, %contig_reads, %human_contigs);

    my $line;
    open (IN, "<$contig_blastn_file");
    while ($line=<IN>) {
      chomp $line;
      my @data   = split /\t/, $line;

      if (!(exists $contig_virus{$data[0]}) || $contig_evalue{$data[0]} < $data[2]*$data[3]){
          $contig_virus{$data[0]}=$data[1];
          $contig_evalue{$data[0]}= $data[2]*$data[3];
          $contig_identity{$data[0]}=$data[2];
      }
    }
    close IN;

    open (IN, "<$human_contigs");
    while ($line=<IN>) {
      my @data = split /\t/, $line;
      $human_contigs{$data[0]}=1
    }
    close IN;


    my $contig_seq = '';
    my $contig_seq2 = '';
    open (IN, "<$contig_fa_file");
    while ($line=<IN>) {
      chomp $line;

      my @temp = split (/ /, $line);
      my $length = substr($temp[1], 4);
      $contig_len{substr($temp[0], 1)} = $length;

      $length = $line;
      $line = <IN>;
      if (!(exists $contig_virus{substr($temp[0], 1)})){
	  $contig_seq  .= $length."\n".$line;
	  $contig_seq2 .= $length."\n".$line if (!exists $human_contigs{substr($temp[0], 1)});
      }
    }
    close IN;

    my $output_prefix = 'results';
    if (length($contig_seq) > 0){
        open (OUT,  ">$output_prefix-novel-contig-1.fa");
        print OUT $contig_seq;
        close OUT;
    }
    if (length($contig_seq2) > 0){
        open (OUT,  ">$output_prefix-novel-contig-2.fa");
        print OUT $contig_seq2;
        close OUT;
    }

    open (IN, "<$psl_file");
    for (my $i=1;$i<=5;$i++) {
      $line = <IN>;
    }
    while ($line=<IN>) {
      chomp $line;
      my @data = split /\t/, $line;
      next if (scalar(@data) < 14);

      $contig_reads{$data[13]} ++;
    }
    close IN;

    my %virus;
    open (OUT, ">$output_prefix-virus-list.txt");
    print OUT "Virus name\tContigs\tContig length (bp)\tIdentities\(\%\)\t#Reads fallen on contigs\n";
    for my $c (sort {$contig_evalue{$b}<=>$contig_evalue{$a}} keys %contig_evalue) {

        if (!exists $virus{$contig_virus{$c}}){
            print OUT "$contig_virus{$c}\t$c\t$contig_len{$c}\t$contig_identity{$c}\t$contig_reads{$c}\n";
            $virus{$contig_virus{$c}} = 1;
        }
    }
    close OUT;

}

sub GetIdir{
    my $PerlEnv  = `perl -e 'print join(\" \",\@INC)'`;
    $PerlEnv .= ' ';

    my $iDir = "";
    foreach(@INC) {
        if (index($PerlEnv, $_.' ') == -1) {
            $iDir .= "-I $_ ";
        }
    }

    return $iDir;
}



############### Check the 'chr' prefix of the reference in bam file ####################

sub GetChrPrefix {
    my $BamFile = shift;

    my $awk_cmd  = '{if(and($2,0x2)>0)print $3}';
    my $ret = `samtools view $BamFile | awk '$awk_cmd' | head -1`;

    if (substr($ret, 0, 3) eq 'chr'){
	return 'chr';
    }else{
        return '';
    }
}

################## Get insert size distribution ######################

sub GetInsertSize{

    my $BamFile = shift;

    my $q=20;  #read quality threshold
    my $c=4;   #times of std from mean

    my @insert_stat;
    my $recordcounter=0;
    my $expected_max=50000;

    open(BAM,"samtools view -h $BamFile |") || die "unable to open $BamFile\n";
    while(<BAM>){
      chomp;
        next if(/^\@RG/);
        next if(/^\@/);
        last if($recordcounter>$expected_max);
        $recordcounter++;

        my @s=split /\t+/;
        next if(@s<10);

        my ($flag,$qual,$dist)=@s[1,4,8];
        next if ($qual<=$q);  #skip low quality mapped reads
        next unless(($flag & 0x0002 || $flag=~/P/) && $dist>=0);

        push @insert_stat, $dist;
    }
    close(BAM);

    my $mean = average(\@insert_stat);
    my $std=stdev(\@insert_stat);

    my @insertsize;
    foreach my $x(@insert_stat){
      next if($x>$mean+5*$std);
      push @insertsize, $x;
    }

    $mean=average(\@insertsize);
    $std=stdev(\@insertsize);

    my ($stdm,$stdp)=(0,0);
    my ($nstdm,$nstdp)=(0,0);
    foreach my $x(@insertsize){
      if($x>$mean){
        $stdp+=($x-$mean)**2;
        $nstdp++;
      }
      else{
        $stdm+=($x-$mean)**2;
        $nstdm++;
      }
    }
    $stdm=sqrt($stdm/($nstdm-1));
    $stdp=sqrt($stdp/($nstdp-1));

    my $upper=$mean+$c*$stdp;
    my $lower=$mean-$c*$stdm;
    $lower=0 if(defined $lower && $lower<0);

    return (int($lower),int($upper),int($mean),$std);
}


################## Get sub sequences that harbor virus integrations ######################

sub GetRefSubSeq{
    my $SVDetect_file   = shift;
    my $flank_region_sz = shift;

    my @boundaries    = `awk \'{print \$4"\\t"\$5"\\t"\$3"\\t"\$7"\\t"\$8"\\t"\$6"\\t"\$9}\' $SVDetect_file`;

    my %ChrHash;
    foreach (@boundaries){
        my @data = split /\t/;

        my $chr = $data[0];
        $chr = $data[3] if ($chr eq 'chrVirus');
        #$chr = substr($chr, 3) if (lc(substr($chr, 0, 3)) eq 'chr');
        $ChrHash{$chr}=1;
    }

    my (@chrs, @lowbounds, @upbounds);
    #for my $i (sort {$a<=>$b} keys %ChrHash){
	#my $chr = "chr$i";
    for my $chr (sort {$a cmp $b} keys %ChrHash){

        my %myhash;
        foreach (@boundaries){
            chomp;
            my @data = split /\t/;

            my ($tmp_chr, $tmp_pos) = ($data[0], $data[1]);
            if ($tmp_chr eq 'chrVirus'){
                $tmp_chr = $data[3];
                $tmp_pos = $data[4];
            }
            next if ($tmp_chr ne $chr);

            my @pos = split('-', $tmp_pos);
            $pos[0] = $pos[0]-$flank_region_sz;
            $pos[0] = 1 if ($pos[0]<1);
            $pos[1] = $pos[1]+$flank_region_sz;
            $myhash{$pos[0]}=$pos[1];
        }

        my (@low_pos, @up_pos);
        for my $low (sort {$a<=>$b} keys %myhash){
            my $up = $myhash{$low};
            my $size = scalar(@low_pos);
            if ($size > 0 && $low<=$up_pos[$size-1]){
                $up_pos[$size-1] = $up if ($up_pos[$size-1] < $up);
            }else{
                push @low_pos, $low;
                push @up_pos,  $up;
            }
        }

        for (my $j=0; $j<scalar(@low_pos); $j++){
            push @chrs, $chr;
            push @lowbounds, $low_pos[$j];
            push @upbounds,  $up_pos[$j];
        }
    }

    return (\@chrs, \@lowbounds, \@upbounds);
}


################## Get unexpored regions that harbor virus integrations ######################

sub FilterRefSubSeq {
    my $chrs 	  = shift;
    my $lowbounds = shift;
    my $upbounds  = shift;
    my $vLocifile = shift;
    my $sclipfile = shift;

    ## Read integration loci detected previsouly
    my @vLoci;
    if (-e $vLocifile && -s $vLocifile){
        open (IN, $vLocifile);
        @vLoci = <IN>;
        close(IN);
    }

    ## Check and remove the overlap between the two detections
    my (@chrs2, @lowbounds2, @upbounds2);

    for (my $i=0; $i<scalar(@$chrs); $i++){
        my $flag = 0;
        for (my $j=1; $j<scalar(@vLoci); $j++){
            my @data = split(/\t/, $vLoci[$j]);

            my ($chr, $pos);
            if ($data[0] ne 'chrVirus'){
                $chr = $data[0];
                $pos = $data[1];
            }else{
                $chr = $data[3];
                $pos = $data[4];
            }
            if ($chr eq @$chrs[$i] &&
            	$pos>=(@$lowbounds[$i]) &&
                $pos<=(@$upbounds[$i])){
            	$flag = 1;
                last;
            }
        }
        if ($flag == 0){
            push @chrs2,      @$chrs[$i];
            push @lowbounds2, @$lowbounds[$i];
            push @upbounds2,  @$upbounds[$i];
        }
    }

#    if (!$sclipfile || scalar(@chrs2) == 0){
    	return (\@chrs2,\@lowbounds2,\@upbounds2);
#    }else{
        ## Check soft-clipped reads
#        my (@chrs3, @lowbounds3, @upbounds3);
#        for (my $i=0; $i<scalar(@chrs2); $i++){
#            my $tmp_chr = substr($chrs2[$i], 3);
#            my $awk_cmd  = '{if(($1=="'.$chrs2[$i].'" || $1=="'.$tmp_chr.'") && $2>='.$lowbounds2[$i].' && $2<='.$upbounds2[$i].')print}';
#            my $soft_num = `awk '$awk_cmd' $sclipfile | wc -l`;
#
#            if ($soft_num >= 2){
#                push @chrs3,      $chrs2[$i];
#                push @lowbounds3, $lowbounds2[$i];
#                push @upbounds3,  $upbounds2[$i];
#            }
#        }
#
#        return (\@chrs3,\@lowbounds3,\@upbounds3);
#    }
}


################## Write regions that harbor virus integrations ######################

sub WriteRefSubSeq{
    my $chrs 	  = shift;
    my $lowbounds = shift;
    my $upbounds  = shift;
    my $ref_seq   = shift;
    my $OutFile   = shift;

    `rm $OutFile` if (-e $OutFile);

    my $ret = `grep '>' $ref_seq | head -1`;

    for (my $j=0; $j<scalar(@$chrs); $j++){
        if ($ret =~ /chr/){
            print "\t@$chrs[$j]:@$lowbounds[$j]-@$upbounds[$j]\n";
            `samtools faidx $ref_seq @$chrs[$j]:@$lowbounds[$j]-@$upbounds[$j] >> $OutFile`;
        }else{
            my $tmp_chr = substr(@$chrs[$j],3);
            print "\t$tmp_chr:@$lowbounds[$j]-@$upbounds[$j]\n";
            `samtools faidx $ref_seq $tmp_chr:@$lowbounds[$j]-@$upbounds[$j] >> $OutFile`;
        }
    }

    `sed 's/[:-]/_/g' $OutFile > $OutFile.2`;
    `mv $OutFile.2 $OutFile`;
}


############### Get additional reads that mapped to the reference####################

sub ExtractMappedReads {

    my $BamFile   = shift;
    my $fastq1    = shift;
    my $fastq2    = shift;
    my $chrs 	  = shift;
    my $lowbounds = shift;
    my $upbounds  = shift;


    my $BamRefPrefix = GetChrPrefix($BamFile);

    my $RefPrefix = 'chr';
    foreach (@$chrs){
        next if (/^chr/);
        $RefPrefix = '';
        last;
    }

    if ($RefPrefix ne $BamRefPrefix){
        for (my $j=0; $j<scalar(@$chrs); $j++){
            if ($BamRefPrefix eq ''){
            	@$chrs[$j] = substr(@$chrs[$j], 3);
            }else{
            	@$chrs[$j] = 'chr'.@$chrs[$j];
            }
        }
    }

    my $prev_chr='';
    my $awk_cmd = '{if(';
    for (my $j=0; $j<scalar(@$chrs); $j++){
    	#next if ($prev_chr eq @$chrs[$j]);

        if ($prev_chr eq ''){
            $awk_cmd .=     '($3=="'.@$chrs[$j].'" && $4>='.@$lowbounds[$j].' && $4<='.@$upbounds[$j].')';
        }else{
            $awk_cmd .= ' || ($3=="'.@$chrs[$j].'" && $4>='.@$lowbounds[$j].' && $4<='.@$upbounds[$j].')';
        }
        $prev_chr = @$chrs[$j];
    }
    my $awk_cmd1 = $awk_cmd . ')print}';
    $awk_cmd .= ')print $1}';

    my $sFILe   = "samtools view -f 2 $BamFile | awk '$awk_cmd' |";
    my %sReads  = ();
    open (IN, $sFILe);
    while (<IN>){
        chomp;
        $sReads{$_} = '';
    }
    close IN;


    my $count = 0;
    $sFILe = "samtools view -f 2 $BamFile | awk '$awk_cmd1' |";
    open (IN, $sFILe);
    open (out1, ">$fastq1");
    open (out2, ">$fastq2");
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




1;