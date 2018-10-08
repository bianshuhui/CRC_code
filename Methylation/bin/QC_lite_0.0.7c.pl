#!/usr/bin/perl
use Inline C;
use warnings;
use strict;
use Getopt::Long;
use File::Basename qw/ basename /;
use List::Util qw/ sum sum0 /;

####
#
#	Description: A simple DP-baseed adapter trimming and
#	             quality control script for pair-ended fastq files.
#	Main Function: 
#	(1) trim 3' end low quality bases
#	(2) remove Illumina adapters (DP algorithm paraphrased from cutadapt in c)
#	(3) brief QC parameter synopsis
#	(4) brief flaw analysis
#	Author: Xylem
#	Last modified on 2015-08-30
#	
###

# R path
my $R_PATH = "R";
# Constants          Values
my $INT_MIXED       = 0;
my $INT_THREE_PRIME = 1;
my $INT_FIVE_PRIME  = 2;


#--------------------------help and options guide--------------------------#
my $usage = <<"END_USAGE" ;
Usage:
	perl $0 [options] R1.fastq.gz R2.fastq.gz
Options:
--quality      INT    Quality trimming cutoff from 3' end. Default: 10.
--stringency   INT    Minimal overlap with adapter to be trimmed. Default: 1.
--length       INT    Minimal length of a kept read. Default: 50.
--lowqual_max  FLOAT  Upper limit of percentage of low quality bases in a kept 
                      read. Range from 0.0 (most stringent) to 1.0 (no abandoning
                      because of low quality). Default: 1.0.
--left_R1      INT    Trim INT bp from left of Read 1 before quality and adapter 
                      trimming. Default: 0.
--left_R2      INT    For Read 2. see above. Default: 0.
--right_R1     INT    Trim INT bp from right of Read 1 before quality and adapter
                      trimming. Default: 0.
--right_R2     INT    For Read 2. Default: 0.

--a1           STR    adapter sequence that possibly occurs in Read 1.
--a2           STR    adapter sequence that possibly occurs in Read 2.
                      Default for --a1 and --a2: "AGATCGGAAGAGC"
--recession_a1 INT    Further inwards recession from adapter position in Read 1.
                      Default: 0.
--recession_a2 INT    Further inwards recession from adapter position in Read 2.
                      Default: 0.
--left_primers STR    oligos that possibly occurs at 5' end of Read 1 or 2.
                      Multiple candidate sequences are separated by comma.
--right_primers STR   oligos that possibly occurs at 3' end of Read 1 or 2.
                      Multiple candidate sequences are separated by comma.

--trimx        INT    further INT bp from 3' end after qualtity and adapter
                      trimming.
 
--symmetry            symmetry read length in a pair. Default: off.
--graph               QC results visualization. Default: off.
--outdir       DIR    the output dir of clean data fastq files.

Preset mode:
--scBS|pbat           abbr. of '--symmetry --left_R1 9 --left_R2 9 --trimx 1'.

--scRNA|tang          abbr. of '
                      --left_primers "ATATGGATCCGGCGCGCCGTCGACTTTTTTTTTTTTTTTTTTTTTTTT,\
                      ATATCTCGAGGGCGCGCCGGATCCTTTTTTTTTTTTTTTTTTTTTTTT" 
                      --right_primers "AAAAAAAAAAAAAAAAAAAAAAAAGTCGACGGCGCGCCGGATCCATAT,\
                      AAAAAAAAAAAAAAAAAAAAAAAAGGATCCGGCGCGCCCTCGAGATAT,"
                      --stringency 3'

If fastq files were not provided definitely, following options will help to
locate fastqs files. 
--indir        DIR    the input dir of the samples. Provided definite fastq files,
--sample	   DIR    sample. If provided, output dir will be changed to 
                            <outdir>/<sample>


END_USAGE

my (
    $indir, $outdir, $sample, 
    $qual_cutoff, $min_len,
    $min_adapter_overlap, $lowqual_ratio,
    $preleft_1,  $preleft_2,
    $preright_1, $preright_2,
    $adapter_1,  $adapter_2, 
    $adapter_recession_1, $adapter_recession_2,
    $lp_str, $rp_str,
    $trimx,
    $symmetry, $pbat, $tang,
    $gzip, $not_gzip, $graph, $help,
); 
GetOptions(
    "outdir=s"        => \$outdir,
    "indir=s"         => \$indir,
	"sample=s"        => \$sample,
    "quality=i"       => \$qual_cutoff,
    "lowqual_max=f"   => \$qual_cutoff,
    "length=i"        => \$min_len,
    "left_R1=i"       => \$preleft_1,
    "left_R2=i"       => \$preleft_2,
    "right_R1=i"      => \$preright_1,
    "right_R2=i"      => \$preright_2,
    "a1=s"            => \$adapter_1,
    "a2=s"            => \$adapter_2,
    "recession_a1=i"  => \$adapter_recession_1,
    "recession_a2=i"  => \$adapter_recession_2,
    "left_primers=s"  => \$lp_str,
    "right_primers=s" => \$rp_str,
    "stringency=i"    => \$min_adapter_overlap,
    "trimx=i"         => \$trimx,
    "symmetry"        => \$symmetry,
    "pbat|scBS"       => \$pbat,
    "tang_protocol|scRNA"   => \$tang,
    "graph"           => \$graph,
    "not_gzip"        => \$not_gzip,
    "help"            => \$help,
);
die $usage if $help;
#
#--------------------------help and options guide--------------------------#

#-default parameters value
$outdir       ||= ".";
$preleft_1    ||= 0;
$preleft_2    ||= 0;
$preright_1   ||= 0;
$preright_2   ||= 0;
$adapter_1    ||= "AGATCGGAAGAGC";
$adapter_2    ||= "AGATCGGAAGAGC";
$lp_str       ||= "";
$rp_str       ||= "";
$adapter_recession_1 ||= 0;
$adapter_recession_2 ||= 0;
$symmetry     ||= 0;
$trimx        ||= 0;
$graph        ||= 0;
$gzip         ||= 1;
$qual_cutoff  ||= 10;
$min_len      ||= 50;
$min_adapter_overlap ||= 1;
$lowqual_ratio ||= 1.0;

# text output instead of gzip compressed files
if ($not_gzip) {
    $gzip  = 0;
}

# scBS-Seq or pbat
if ($pbat) {
    $symmetry  = 1;
    $preleft_1 = 9;
    $preleft_2 = 9;
    $trimx     = 1;
    $min_adapter_overlap = 1;
}

# Tang protocol
if ($tang) {
    $lp_str = "ATATGGATCCGGCGCGCCGTCGACTTTTTTTTTTTTTTTTTTTTTTT,"
             ."ATATCTCGAGGGCGCGCCGGATCCTTTTTTTTTTTTTTTTTTTTTTT,";
    $rp_str = "AAAAAAAAAAAAAAAAAAAAAAAGTCGACGGCGCGCCGGATCCATAT,"
             ."AAAAAAAAAAAAAAAAAAAAAAAGGATCCGGCGCGCCCTCGAGATAT,";
    $min_adapter_overlap = 3;
}

#-trim_slash 
$outdir = trim_slash($outdir);


#--------------------------input files----------------------------------  #

my ($in_fq1, $in_fq2) = @ARGV;
my ($bn_1, $bn_2);
my ($out_fq1, $out_fq2);

if ( $in_fq1 && $in_fq2 ) {
    if ( (! -e $in_fq1) || (! -e $in_fq2) ) {
        die "fastqs files not found\n";
    }
    ($bn_1, $bn_2) = map { basename($_) } ($in_fq1, $in_fq2);
    if ($bn_1 =~ /(.*)\.(fastq|fq)\.gz$/) {
        $out_fq1 = "$outdir/$1"."_val_1.fq";
    }
    if ($bn_2 =~ /(.*)\.(fastq|fq)\.gz$/) {
        $out_fq2 = "$outdir/$1"."_val_2.fq";
    }
}
elsif ( $indir && $sample ) {
    $indir  = trim_slash($indir);
    $sample = trim_slash($sample);
    my @fastqs;
    opendir DH, "$indir/$sample" or die $!;
    while (my $file = readdir DH) {
        if ($file =~ /(fastq|fq).gz$/) {
            push @fastqs, $file;
        }
    }
    closedir DH;
    ($bn_1, $bn_2) = sort @fastqs;
    ($in_fq1, $in_fq2) = map { "$indir/$sample/".$_ } ($bn_1, $bn_2);

    if ( (! -e $in_fq1) || (! -e $in_fq2) ) {
        die "fastqs files not found\n";
    }

    $outdir = "$outdir/$sample";
    ($out_fq1, $out_fq2) = map { "$outdir/$sample.".$_.".clean.fq" } qw/ R1 R2 /;
}
else {
    die "Fastq files not available\n".$usage;
}
    
system "mkdir", "-p", "$outdir" unless(-d "$outdir"); 

print STDERR "Input fastqs: $in_fq1,$in_fq2\n";
print STDERR "Output fastqs: $out_fq1,$out_fq2\n";

#-create the Log file
my $log_prefix = ($sample) ? $sample : $bn_1;
open LOG,">$outdir/$log_prefix.paired.QC.log" or die $!;
my $now = scalar localtime;
print LOG "$log_prefix begin at: $now\n";



my (
    $paired, $qual_base,
    $error_rate, 
    $qual_thres,
);
#-default parameters value
$paired      ||= 1;
$qual_base   ||= 33;
$error_rate  ||= 0.1;
$qual_thres  ||= "20,30";

my $adapter_len_1 = length $adapter_1;
$adapter_len_1 += $preleft_2 if ($symmetry);
my $adapter_len_2 = length $adapter_2;
$adapter_len_2 += $preleft_1 if ($symmetry);


my @thres = split /\,/, $qual_thres;
my @left_primers  = split /\,/, $lp_str;
my @right_primers = split /\,/, $rp_str;

my @left_primer_lens  = map {length $_} @left_primers;
my @right_primer_lens = map {length $_} @right_primers;

# declairation and initialization
my ($total_reads, $total_bases, $clean_reads, $clean_bases,) = (0)x4; 
my ($total_pairs, $total_bases_1, $total_bases_2, 
    $clean_pairs, $clean_bases_1, $clean_bases_2,) = (0)x6; 

my (%seen_input_readlen_1, %seen_input_readlen_2);
my (%seen_clean_readlen_1, %seen_clean_readlen_2);

## base compositon (AoH)
my (@post_base_dist_1, @post_base_dist_2,);

## quality distribtion (AoH)
my (@post_qual_dist_1, @post_qual_dist_2);

## flaws anaylsis
my %flaw_of = (
    "three_lowqual_1" => 0,
    "three_lowqual_2" => 0,
    "three_adapter_1" => 0,
    "three_adapter_2" => 0,
);
my %discard_because_of = (
    "three_lowqual"   => 0,    
    "general_lowqual" => 0,    
    "three_adapter"   => 0,    
);
my @adapter_1_counts;
my @adapter_2_counts;

#____________________________________processing begin_______________________________________________#


#____________________________________for the default paired end sequencing______________________________________#

if ($paired) {

	#-open the input files
	open IN_1, "gzip -dc $in_fq1 |" or die $!;
	open IN_2, "gzip -dc $in_fq2 |" or die $!;

	#-open the output files
	open OUT_1, "> $out_fq1" or die $!;
	open OUT_2, "> $out_fq2" or die $!;
	

	#---------------------------read in the reads information-----------------------#

	while (1) {
		
		#-get the reads and corresponding information in each 4 lines
		my $line1_1 = <IN_1>;
		my $line1_2 = <IN_1>;
		my $line1_3 = <IN_1>;
		my $line1_4 = <IN_1>;

		my $line2_1 = <IN_2>;
		my $line2_2 = <IN_2>;
		my $line2_3 = <IN_2>;
		my $line2_4 = <IN_2>;

		# check the end of the file
		last unless (defined($line1_1) and defined($line2_1));
		chomp ($line1_1,$line1_2,$line1_3,$line1_4,$line2_1,$line2_2,$line2_3,$line2_4);
		$total_pairs += 1;
        
        
        # parse array    
        my $seq_1 = uc $line1_2;
        my $seq_2 = uc $line2_2;
        my @qual_1 = split //, $line1_4;
        my @qual_2 = split //, $line2_4;
        my @int_qual_1 = map { ord($_) - $qual_base } @qual_1;
        my @int_qual_2 = map { ord($_) - $qual_base } @qual_2;

        my $seqlen_1 = length $seq_1;
        my $seqlen_2 = length $seq_2;
        $total_bases_1 += $seqlen_1;
        $total_bases_2 += $seqlen_2;
        
        $seen_input_readlen_1{$seqlen_1} = 1;
        $seen_input_readlen_2{$seqlen_2} = 1;

        #-------------------------------------------------------------
        # data grooming begins
        # preliminary trim
        my ($ad1_leading, $ad2_leading);

        if ($preleft_1 || $preright_1) {
            $seqlen_1 = $seqlen_1 - $preright_1 - $preleft_1;
            die "ill-defined left or right trim parameters" if ($seqlen_1 <=0);
            $ad2_leading = rev_comp(substr($seq_1,0,$preleft_1));
            $seq_1 = substr $seq_1, $preleft_1, $seqlen_1;
            @int_qual_1 = splice @int_qual_1, $preleft_1, $seqlen_1;
        }
        if ($preleft_2 || $preright_2) {
            $seqlen_2 = $seqlen_2 - $preright_2 - $preleft_2;
            die "ill-defined left or right trim parameters" if ($seqlen_2 <=0);
            $ad1_leading = rev_comp(substr($seq_2,0,$preleft_2));
            $seq_2 = substr $seq_2, $preleft_2, $seqlen_2;
            @int_qual_2 = splice @int_qual_2, $preleft_2, $seqlen_2;
        }
            

		# trim the low quality 3' end low quality base
		my $qual_removal_start_1 = label3lowqual(\@int_qual_1, $qual_cutoff); 
		my $qual_removal_start_2 = label3lowqual(\@int_qual_2, $qual_cutoff); 

        if ($qual_removal_start_1 < $#int_qual_1) {
            ++$flaw_of{"three_lowqual_1"};
            $seqlen_1 = $qual_removal_start_1 + 1;
        }
        if ($qual_removal_start_2 < $#int_qual_2) {
            ++$flaw_of{"three_lowqual_2"};
            $seqlen_2 = $qual_removal_start_2 + 1;
        }

        if (($seqlen_1 < $min_len) || ($seqlen_2 < $min_len)) {
            ++$discard_because_of{"three_lowqual"};
            next;
        }
    
        if ($symmetry) {
            ($seqlen_1 > $seqlen_2) ? ($seqlen_1 = $seqlen_2) : ($seqlen_2 = $seqlen_1);
        }

        # discard reads with more than 50% of low-quality bases
        if ($lowqual_ratio < 1) {
            my $lowqual_base_1 = scalar grep {$_ < $qual_cutoff} @int_qual_1[0..($seqlen_1-1)];
            my $lowqual_base_2 = scalar grep {$_ < $qual_cutoff} @int_qual_2[0..($seqlen_2-1)];
            if ($lowqual_base_1 > $seqlen_1 * $lowqual_ratio || $lowqual_base_2 > $seqlen_2 * $lowqual_ratio) {
                ++$discard_because_of{"general_lowqual"};
                next;
            }
        }

        # trim illumina adapter if needed
        my $ad_1 = ($preleft_2 && $symmetry) ? $ad1_leading.$adapter_1 : $adapter_1;
        my $ad_2 = ($preleft_1 && $symmetry) ? $ad2_leading.$adapter_2 : $adapter_2;

        my ($ad_start_1, $ad_stop_1, $ad_cost_1, $ad_match_1) 
            = locate_single_adapter($ad_1, $seq_1, $adapter_len_1, $seqlen_1,
                                    $error_rate, $INT_THREE_PRIME, $min_adapter_overlap);
        my ($ad_start_2, $ad_stop_2, $ad_cost_2, $ad_match_2) 
            = locate_single_adapter($ad_2, $seq_2, $adapter_len_2, $seqlen_2,
                                    $error_rate, $INT_THREE_PRIME, $min_adapter_overlap);

        if ($ad_start_1 != -1) {
            ++$flaw_of{"three_adapter_1"};
            ++$adapter_1_counts[$ad_start_1];
            $seqlen_1 = $ad_start_1 - $adapter_recession_1;
            $seqlen_1 = 0 if ($seqlen_1 < 0);
        }
        if ($ad_start_2 != -1) {
            ++$flaw_of{"three_adapter_2"};
            ++$adapter_2_counts[$ad_start_2];
            $seqlen_2 = $ad_start_2 - $adapter_recession_2;
            $seqlen_2 = 0 if ($seqlen_2 < 0);
        }
        
        if ( ($seqlen_1 < $min_len) || ($seqlen_2 < $min_len) ) {
            ++$discard_because_of{"three_adapter"};
            next;
        }

        if ($symmetry) {
            ($seqlen_1 > $seqlen_2) ? ($seqlen_1 = $seqlen_2) : ($seqlen_2 = $seqlen_1);
        }
      
        
        my ($begin_1, $begin_2) = (0)x2;

        if (scalar @right_primers) {
            #my ($rp_start_1, $rp_stop_1, $rp_cost_1,$rp_match_1) = (-1,-1,1000,0);
            #my ($rp_start_2, $rp_stop_2, $rp_cost_2,$rp_match_2) = (-1,-1,1000,0);
            for my $i (0..$#right_primers) {
                my $rp = $right_primers[$i];
                my $rp_len = $right_primer_lens[$i];
                my ($rp_start_1, $rp_stop_1, $rp_cost_1,$rp_match_1)
                    = locate_single_adapter($rp, $seq_1, $rp_len, $seqlen_1,
                                            $error_rate, $INT_THREE_PRIME, $min_adapter_overlap);
                if ($rp_start_1 != -1) {
                    ++$flaw_of{ "rp".$i."_1" };
                    $seqlen_1 = $rp_start_1;
                    $seqlen_1 = 0 if ($seqlen_1 < 0);
                    last;
                }
            }
            for my $i (0..$#right_primers) {
                my $rp = $right_primers[$i];
                my $rp_len = $right_primer_lens[$i];
                my ($rp_start_2, $rp_stop_2, $rp_cost_2,$rp_match_2)
                    = locate_single_adapter($rp, $seq_2, $rp_len, $seqlen_2,
                                            $error_rate, $INT_THREE_PRIME, $min_adapter_overlap);
                if ($rp_start_2 != -1) {
                    ++$flaw_of{ "rp".$i."_2" };
                    $seqlen_2 = $rp_start_2;
                    $seqlen_2 = 0 if ($seqlen_2 < 0);
                    last;
                }
            }
            if ( ($seqlen_1 < $min_len) || ($seqlen_2 < $min_len) ) {
                ++$discard_because_of{"right_primers"};
                next;
            }

            if ($symmetry) {
                ($seqlen_1 > $seqlen_2) ? ($seqlen_1 = $seqlen_2) : ($seqlen_2 = $seqlen_1);
            }
        }
        
        if (scalar @left_primers) {
            #my ($lp_start_1, $lp_stop_1, $lp_cost_1,$lp_match_1) = (-1,-1,1000,0);
            #my ($lp_start_2, $lp_stop_2, $lp_cost_2,$lp_match_2) = (-1,-1,1000,0);
            for my $i (0..$#left_primers) {
                my $lp = $left_primers[$i];
                my $lp_len = $left_primer_lens[$i];
                my ($lp_start_1, $lp_stop_1, $lp_cost_1,$lp_match_1)
                    = locate_single_adapter($lp, $seq_1, $lp_len, $seqlen_1,
                                            $error_rate, $INT_FIVE_PRIME, $min_adapter_overlap);
                if ($lp_stop_1 != -1) {
                    ++$flaw_of{ "lp".$i."_1" };
                    $begin_1 = $lp_stop_1 + 1;
                    $seqlen_1 -= $begin_1;
                    $seqlen_1 = 0 if ($seqlen_1 < 0);
                    last;
                }
            }
            for my $i (0..$#left_primers) {
                my $lp = $left_primers[$i];
                my $lp_len = $right_primer_lens[$i];
                my ($lp_start_2, $lp_stop_2, $lp_cost_2,$lp_match_2)
                    = locate_single_adapter($lp, $seq_2, $lp_len, $seqlen_2,
                                            $error_rate, $INT_FIVE_PRIME, $min_adapter_overlap);
                if ($lp_stop_2 != -1) {
                    ++$flaw_of{ "lp".$i."_2" };
                    $begin_2 = $lp_stop_2 + 1;
                    $seqlen_2 -= $begin_2;
                    $seqlen_2 = 0 if ($seqlen_2 < 0);
                    last;
                }
            }
            if ( ($seqlen_1 < $min_len) || ($seqlen_2 < $min_len) ) {
                ++$discard_because_of{"left_primers"};
                next;
            }

            if ($symmetry) {
                ($seqlen_1 > $seqlen_2) ? ($seqlen_1 = $seqlen_2) : ($seqlen_2 = $seqlen_1);
            }
        }

        
 
        # final 
        if ($trimx) {
            $seqlen_1 -= $trimx;
            $seqlen_2 -= $trimx;
            next if ($seqlen_1 <= 0 || $seqlen_2 <= 0);
        }
 
        $seq_1 = substr $seq_1, $begin_1, $seqlen_1;
        @int_qual_1 = splice @int_qual_1, $begin_1, $seqlen_1;
        $seq_2 = substr $seq_2, $begin_2, $seqlen_2;
        @int_qual_2 = splice @int_qual_2, $begin_2, $seqlen_2;

        #------------------------ end of data grooming --------------------------
        # post base stat
        base_stat($seq_1, $seqlen_1, \@post_base_dist_1);
        base_stat($seq_2, $seqlen_2, \@post_base_dist_2);
        # post quality stat
        qual_stat(\@int_qual_1, $seqlen_1, \@post_qual_dist_1);
        qual_stat(\@int_qual_2, $seqlen_2, \@post_qual_dist_2);
        
        $seen_clean_readlen_1{$seqlen_1} = 1;
        $seen_clean_readlen_2{$seqlen_2} = 1;

        $clean_pairs += 1;
        $clean_bases_1 += $seqlen_1;
        $clean_bases_2 += $seqlen_2;
        
        $line1_2 = $seq_1;
        $line1_4 = join "", (map {chr($_ + $qual_base)} @int_qual_1);
        $line2_2 = $seq_2;
        $line2_4 = join "", (map {chr($_ + $qual_base)} @int_qual_2);
        
        
		#out put the remanent reads
		print OUT_1 "$line1_1\n$line1_2\n$line1_3\n$line1_4\n";
        print OUT_2 "$line2_1\n$line2_2\n$line2_3\n$line2_4\n";
	}
	
	
	#-close the file handle
	close IN_1;close IN_2;close OUT_1;close OUT_2;
    if ($gzip) {
        unlink "$out_fq1.gz" if (-e "$out_fq1.gz");
        unlink "$out_fq2.gz" if (-e "$out_fq2.gz");
        system "gzip", "$out_fq1";
        system "gzip", "$out_fq2";
    }
	#---------------------------------read in done----------------------------#
	
    #------------------------------- BEGIN Synopsis  ----------------------------#
    
    pe_synopsis();
    pe_base_composition_and_quality();
    pe_adapter_distribution();
}
#____________________________________for the default pair end sequencing done______________________________________#





#________________________________________________Subrutines begin___________________________________________________#

# summerise total/clean reads/bases, clean ratio, post- GC content, Q20/30
sub pe_synopsis {
    # ----------------------------- input reads --------------------------------------------
    # input read length
    my @input_lens_1 = sort {$b <=> $a} keys %seen_input_readlen_1;
    my @input_lens_2 = sort {$b <=> $a} keys %seen_input_readlen_2;
    die "read no sequence\n" if (@input_lens_1 == 0);
    my ($inputlen_max_1, $inputlen_min_1) = @input_lens_1[0,-1];
    my ($inputlen_max_2, $inputlen_min_2) = @input_lens_2[0,-1];
    print LOG "Input read 1 length:\t[$inputlen_min_1,$inputlen_max_1]\n";
    print LOG "Input read 2 length:\t[$inputlen_min_2,$inputlen_max_2]\n";

    # total reads
    $total_reads = 2 * $total_pairs;
    # total bases
    $total_bases = $total_bases_1 + $total_bases_2;
    print LOG "Total pairs:\t$total_pairs\n";
    print LOG "Total reads:\t$total_reads\n";
    print LOG "Total bases:\t$total_bases\n";

    
    # ------------------------------- clean reads --------------------------------------------
    # clean read length
    my @clean_lens_1 = sort {$b <=> $a} keys %seen_clean_readlen_1;
    my @clean_lens_2 = sort {$b <=> $a} keys %seen_clean_readlen_2;
    my ($cleanlen_max_1, $cleanlen_min_1) = @clean_lens_1[0,-1];
    my ($cleanlen_max_2, $cleanlen_min_2) = @clean_lens_2[0,-1];
    print LOG "Clean read 1 length:\t[$cleanlen_min_1,$cleanlen_max_1]\n";
    print LOG "Clean read 2 length:\t[$cleanlen_min_2,$cleanlen_max_2]\n";

    # clean reads
    $clean_reads = 2 * $clean_pairs;
    # clean bases
    $clean_bases = $clean_bases_1 + $clean_bases_2;
    # clean ratio
    my $clean_ratio = $clean_pairs / $total_pairs;
    print LOG "Clean pairs:\t$clean_pairs\n";
    print LOG "Clean reads:\t$clean_reads\n";
    print LOG "Clean bases:\t$clean_bases\n";
    print LOG "Clean rate:\t$clean_ratio\n";
    
    # ------------------------------ clean base and quality stat -------------------------
    open QUALLOG,">$outdir/$log_prefix.seq.quality.txt"      or die $!;
    open BASELOG,">$outdir/$log_prefix.base.composition.txt" or die $!;

    # pre-QC quality and GC content
    my ($post_1_GC_bases, $post_2_GC_bases,) = (0)x2;
    my (%qualcount_1, %qualcount_2);

    # position iteration
    for my $i (0..($cleanlen_max_1-1)) {
        my $pos = $i + 1;
        # base
        my $href_basecount = $post_base_dist_1[$i];
        my $base_sum = sum0( values %{$href_basecount} );
        for my $base ('A','T','C','G') {
            if (! exists ${$href_basecount}{$base}) {
                $href_basecount->{$base} = 0;
            }
            my $basecount = $href_basecount->{$base} ;
            my $content = $basecount / $base_sum;
            print BASELOG "Read1\t$pos\t$base\t$content\n";
        }
        $post_1_GC_bases += sum0( @{$href_basecount}{ qw/ G C / } );
        
        # qual
        my $href_qualcount = $post_qual_dist_1[$i];
        my ($qual_mean, $qual_sd);
        weight_hash_stat($href_qualcount, \$qual_mean, \$qual_sd);
        print QUALLOG "Read1\t$pos\t$qual_mean\t$qual_sd\n";
    
        while( my ($qual,$count) = each %{$href_qualcount}) {
            $qualcount_1{$qual} += $count;
        }
    }
    for my $i (0..($cleanlen_max_2-1)) {
        my $pos = $i + 1;
        # base
        my $href_basecount = $post_base_dist_2[$i];
        my $base_sum = sum0( values %{$href_basecount} );
        for my $base ('A','T','C','G') {
            if (! exists ${$href_basecount}{$base}) {
                $href_basecount->{$base} = 0;
            }
            my $basecount = $href_basecount->{$base} ;
            my $content = $basecount / $base_sum;
            print BASELOG "Read2\t$pos\t$base\t$content\n";
        }
        $post_2_GC_bases += sum0( @{$href_basecount}{ qw/ G C / } );
        
        # qual
        my $href_qualcount = $post_qual_dist_2[$i];
        my ($qual_mean, $qual_sd);
        weight_hash_stat($href_qualcount, \$qual_mean, \$qual_sd);
        print QUALLOG "Read2\t$pos\t$qual_mean\t$qual_sd\n";
    
        while( my ($qual,$count) = each %{$href_qualcount}) {
            $qualcount_2{$qual} += $count;
        }
    }
    
    # Q20/Q30 or something like that
    for my $threshold (@thres) {
        my @above_quals_1 = grep { $_ >= $threshold } (keys %qualcount_1);
        my $above_count_1 = sum(@qualcount_1{ @above_quals_1 });
        my $above_rate_1  = $above_count_1 / $clean_bases_1;
        my @above_quals_2 = grep { $_ >= $threshold } (keys %qualcount_2);
        my $above_count_2 = sum(@qualcount_2{ @above_quals_2 });
        my $above_rate_2  = $above_count_2 / $clean_bases_2;
        printf LOG "Q%d:\t%.2f,%.2f\n", $threshold, $above_rate_1, $above_rate_2;
    }
    
    my $post_1_GC_content = $post_1_GC_bases / $clean_bases_1;
    my $post_2_GC_content = $post_2_GC_bases / $clean_bases_2;
    printf LOG "GC content:\t%.3f,%.3f\n", $post_1_GC_content, $post_2_GC_content;
    
    print LOG "\n";
    print LOG "Read 1 with low quality 3' end:\t$flaw_of{three_lowqual_1}\n";
    print LOG "Read 2 with low quality 3' end:\t$flaw_of{three_lowqual_2}\n";
    print LOG "Read 1 with 3' end adapter contamination:\t$flaw_of{three_adapter_1}\n";
    print LOG "Read 2 with 3' end adapter contamination:\t$flaw_of{three_adapter_2}\n";
    print LOG "Discarded pairs because of 3' end low quality:\t$discard_because_of{three_lowqual}\n";
    print LOG "Discarded pairs because of general low quality:\t$discard_because_of{general_lowqual}\n";
    print LOG "Discarded pairs because of adapter contamination:\t$discard_because_of{three_adapter}\n";
    print LOG "\n";

    if (scalar @left_primers ) {
        for my $i (0..$#left_primers) {
            if (exists $flaw_of{"lp".$i."_1"}) {
                my $primer_count = $flaw_of{"lp".$i."_1"};
                my $one_based_number = $i + 1;
                print LOG "Read 1 with 5' primer $one_based_number:\t$primer_count\n";
            }
            if (exists $flaw_of{"lp".$i."_2"}) {
                my $primer_count = $flaw_of{"lp".$i."_2"};
                my $one_based_number = $i + 1;
                print LOG "Read 2 with 5' primer $one_based_number:\t$primer_count\n";
            }
        }
        if (exists $discard_because_of{"left_primers"}) {
            print LOG "Discarded pairs because of 5' primers contamination:\t$discard_because_of{left_primers}\n";
        }
    }
    if (scalar @right_primers ) {
        for my $i (0..$#right_primers) {
            if (exists $flaw_of{"rp".$i."_1"}) {
                my $primer_count = $flaw_of{"rp".$i."_1"};
                my $one_based_number = $i + 1;
                print LOG "Read 1 with 3' primer $one_based_number:\t$primer_count\n";
            }
            if (exists $flaw_of{"rp".$i."_2"}) {
                my $primer_count = $flaw_of{"rp".$i."_2"};
                my $one_based_number = $i + 1;
                print LOG "Read 2 with 3' primer $one_based_number:\t$primer_count\n";
            }
        }
        if (exists $discard_because_of{"right_primers"}) {
            print LOG "Discarded pairs because of 3' primers contamination:\t$discard_because_of{right_primers}\n";
        }
    }

    my $end_time = scalar localtime;
    print LOG "$log_prefix end at: $end_time\n";
    close LOG;
    close QUALLOG;
    close BASELOG;
    # ------------------------------------- end of log ---------------------------------------
}

# base composition and quality analysis
sub pe_base_composition_and_quality {
    
   

    if ($graph) {
        # base composition graph
        my $base_comp_plot_script = << "END_BASE_GRAPH_SCRIPT";
library(grid)
library(ggplot2)
setwd("$outdir")
base <- read.delim("$log_prefix.base.composition.txt",header=F)
colnames(base) <- c("read","base","nucleotide","freq")
ggplot(base,aes(base,freq,color=nucleotide)) + geom_line() + facet_grid(~read) + theme_bw() + theme(panel.margin=unit(5,"lines"),legend.position=c(.5,.5)) + labs(title="$log_prefix",x="Base",y="Frequency") + coord_cartesian(ylim=c(0,1))
ggsave("$log_prefix.base.composition.pdf",plot=last_plot(), width=16,height=5,units="in",dpi=300)
quit(save="no")
END_BASE_GRAPH_SCRIPT
        open BASEGRAPH, "| $R_PATH --vanilla --slave" or die $!;
        print BASEGRAPH $base_comp_plot_script;
        close BASEGRAPH;

        # sequence quality graph
        my $qual_plot_script = << "END_QUAL_GRAPH_SCRIPT";
library(grid)
library(ggplot2)
setwd("$outdir")
qual <- read.delim("$log_prefix.seq.quality.txt",header=F)
colnames(qual) <- c("read","base","qual","sd")
ggplot(qual,aes(base,qual)) + geom_line(color="blue") + facet_grid(~read) + theme_bw() + labs(title="$log_prefix",x="Base in read (bp)",y="Quality") + coord_cartesian(ylim=c(0,40))
ggsave("$log_prefix.seq.quality.pdf",plot=last_plot(), width=16,height=5,units="in",dpi=300)
quit(save="no")
END_QUAL_GRAPH_SCRIPT
        open QUALGRAPH, "|$R_PATH --vanilla --slave" or die $!;
        print QUALGRAPH $qual_plot_script;
        close QUALGRAPH;
        unlink("$outdir/Rplots.pdf");
    }
}


sub pe_adapter_distribution {
    
    # clean read length
    my @clean_lens_1 = sort {$b <=> $a} keys %seen_clean_readlen_1;
    my @clean_lens_2 = sort {$b <=> $a} keys %seen_clean_readlen_2;
    my ($cleanlen_max_1, $cleanlen_min_1) = @clean_lens_1[1,-1];
    my ($cleanlen_max_2, $cleanlen_min_2) = @clean_lens_2[1,-1];

    open ADAPTDIST,">$outdir/$log_prefix.adapter.content.txt" or die $!;
   
    # read 1 
    for my $i (0..($cleanlen_max_1-1)) {
        my $pos = $i + 1;
        my ($count, $percentage) = (0)x2;
        if ( $adapter_1_counts[$i] ) {
            $count = $adapter_1_counts[$i];
            $percentage = $count / $total_pairs;
        }
        print ADAPTDIST "Read1\t$pos\t$count\t$percentage\n";
    }
    
    # read 2
    for my $i (0..($cleanlen_max_2-1)) {
        my $pos = $i + 1;
        my ($count, $percentage) = (0)x2;
        if ( $adapter_2_counts[$i] ) {
            $count = $adapter_2_counts[$i];
            $percentage = $count / $total_pairs;
        }
        print ADAPTDIST "Read2\t$pos\t$count\t$percentage\n";
    }
    close ADAPTDIST;

    if ($graph) {
        # adapter content distribution graph
        my $adapter_dist_graph_script = << "END_ADAPTER_DIST_GRAPH_SCRIPT";
library(grid)
library(ggplot2)
setwd("$outdir")
ad_dist <- read.delim("$log_prefix.adapter.content.txt",header=F)
colnames(ad_dist) <- c("read","base","count","percentage")
ggplot(ad_dist,aes(base,percentage,color=read)) + geom_line() + facet_grid(~read) + theme_bw() + theme(panel.margin=unit(5,"lines"),legend.position=c(.5,.5)) + labs(title="$log_prefix",x="Base",y="%Adapter") + coord_cartesian(ylim=c(0,1))
ggsave("$log_prefix.adapter.content.pdf",plot=last_plot(), width=16,height=5,units="in",dpi=300)
quit(save="no")
END_ADAPTER_DIST_GRAPH_SCRIPT
        open ADAPTGRAPH, "| $R_PATH --vanilla --slave" or die $!;
        print ADAPTGRAPH $adapter_dist_graph_script;
        close ADAPTGRAPH;
        unlink("$outdir/Rplots.pdf");
    }

}

#-dir trimming
sub trim_slash {
	my ($dir) = @_;
	chop $dir if ($dir =~ /\/$/);
	return $dir;
}


sub rev_comp {
    my ($seq) = @_;
    my $rc = reverse $seq;
    $rc =~ tr/ATCGN/TAGCN/;
    return $rc;
}


# remove low qual base at 3' end
sub label3lowqual {
    my ($aref, $thres) = @_;
    my @reduction = map {$_ - $thres} @{$aref};
    my $partial_sum = 0;
    for my $i (reverse (0..$#reduction)) {
        $partial_sum += $reduction[$i];
        return $i if ($partial_sum > 0);
    }
    return $#reduction;
}


# collect base composition and deposit in a assigned hash
sub base_stat {
    my ($seq, $len, $AoHref_to_basecount) = @_;
    for my $i (0..($len-1)) {
        my $base = substr $seq, $i, 1;
        $AoHref_to_basecount->[$i]{$base} += 1;
    }
}

# collect sum of base quality
sub qual_stat {
    my ($aref_to_qual, $len, $AoHref_to_qualcount) = @_;
    for my $i (0..($len-1)) {
        my $basequal = $aref_to_qual->[$i];
        $AoHref_to_qualcount->[$i]{$basequal} += 1;
    }
}

# 
sub weight_hash_stat {
    my ($href, $p_mean, $p_sd) = @_;
    my $sum = 0;
    my $n   = 0;
    my $square_sum = 0;
    while (my ($obs, $weight) = each %{$href}) {
        $n   += $weight;
        $sum += $obs * $weight;
        $square_sum += $obs * $obs * $weight;
    }
    $$p_mean = $sum / $n;
    $$p_sd   = sqrt(($square_sum  - $sum * $sum / $n) / ($n - 1)) if ($n >1);
}

#########################   END OF PERL PART    #########################
#########################################################################


__DATA__
__C__

#define DELETION_COST 1
#define INSERTION_COST 1
#define MISMATCH_COST 1
#define MATCH_COST 0

#define INT_SEMI_GLOB 0
#define INT_THREE_PRIME 1
#define INT_FIVE_PRIME 2

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))


void locate_single_adapter (const char *s1, const char *s2, int m, int n, double error_rate_thres, int adapter_type_int,  int min_overlap) {

	/*
    # DP Matrix:
                 read
                s2 (j)
              ----------> n
             |
     s1 (i)  |
     adapter |
             V
             m
	*/

    int can_start_in_s1 = 0;
    int can_start_in_s2 = 0;
    int can_stop_in_s1 = 0;
    int can_stop_in_s2 = 0;

    if (adapter_type_int == INT_THREE_PRIME) {
        can_start_in_s2 = 1;
        can_stop_in_s1  = 1;
        can_stop_in_s2  = 1;
    }
    else if (adapter_type_int == INT_FIVE_PRIME) {
        can_start_in_s1 = 1;
        can_start_in_s2 = 1;
        can_stop_in_s2  = 1;
    }
    else if (adapter_type_int == INT_SEMI_GLOB) {
        can_start_in_s2 = 1;
        can_start_in_s1 = 1;
        can_stop_in_s1  = 1;
        can_stop_in_s2  = 1;
    }
    else {
        // return negative default
        Inline_Stack_Vars;

        Inline_Stack_Reset;
        Inline_Stack_Push(sv_2mortal(newSViv(-1)));
        Inline_Stack_Push(sv_2mortal(newSViv(-1)));
        Inline_Stack_Push(sv_2mortal(newSViv(m+n)));
        Inline_Stack_Push(sv_2mortal(newSViv(0)));
        Inline_Stack_Done;
        Inline_Stack_Return(4);
    }

    int i, j;

    /* maximum no. of errors */
    int k = (int)(m * error_rate_thres);
    
    // Determine largest and smallest column
    int max_n = n;
    int min_n = 0;

    if (can_start_in_s2 == 0) {
        // costs can only get worse after column m 
        max_n = min(n, m + k);
    }
    if (can_stop_in_s2 == 0) { 
        // for suffix case
        min_n = max(0, n - m - k);
    }
    
    /* fill the auxilary columns array: @Cost, @Match, @Origin
    Four cases:
    not startin1, not startin2: c(i,j) = max(i,j); origin(i, j) = 0
        startin1, not startin2: c(i,j) = j       ; origin(i, j) = min(0, j - i)
    not startin1,     startin2: c(i,j) = i       ; origin(i, j) =
        startin1,     startin2: c(i,j) = min(i,j) 
    */

    int cost[m+1];
    int match[m+1];
    int origin[m+1];

    if ((can_start_in_s1==0) && (can_start_in_s2==0)) {
        for (i=0; i<=m; i++) {
            match[i] = 0;
            cost[i]  = max(i, min_n);
            origin[i] = 0;
        }
    }
    else if ((can_start_in_s1) && (can_start_in_s2==0)) {
        for (i=0; i<=m; i++) {
            match[i] = 0;
            cost[i] = min_n;
            origin[i] = min(0, min_n - i);
        }
    }
    else if ((can_start_in_s1==0) && (can_start_in_s2)) {
        for (i=0; i<=m; i++) {
            match[i] = 0;
            cost[i] = i;
            origin[i] = max(0, min_n - i);
        }
    }
    else {
        for (i=0; i<=m; i++) {
            match[i] = 0;
            cost[i] = min(i, min_n);
            origin[i] = min_n - i;
        }
    }

    // alignment variable
    int best_cost   = m + n;
    int best_match  = 0;
    int best_origin = 0;
    int best_stop_s1 = m;
    int best_stop_s2 = n;

	// break cutoff
    int last = (can_start_in_s1) ? m : k + 1;


	// iterate over columns
	for (j = min_n + 1; j <= max_n; ++j) {
        // diagnal parameters
        int diag_prev_cost   = cost[0];
        int diag_prev_match  = match[0];
        int diag_prev_origin = origin[0];

		if (can_start_in_s2) {
            origin[0] = j;
        }
        else {
            cost[0] = j;
        }
        // printf("%d\t%d\t%d\t%d\n",j,cost[0],match[0],origin[0]);
		for (i = 1; i <= last; ++i) {
			int cur_dist = (s1[i-1] == 'N')     ? 0 :
			               (s2[j-1] == 'N')     ? 0 :
			               (s1[i-1] == s2[j-1]) ? 0 : 1;
			int cost_candid_diag = diag_prev_cost + cur_dist;
			int cost_candid_del  = cost[i]   + DELETION_COST;
			int cost_candid_ins  = cost[i-1] + INSERTION_COST;
            
            // left parameters
            int helper_cost   = cost[i];
            int helper_match  = match[i];
            int helper_origin = origin[i];
            
			if (cost_candid_diag <= cost_candid_del && cost_candid_diag <= cost_candid_ins) {
				// MATCH or MISMATCH
				cost[i]   = cost_candid_diag;
				origin[i] = diag_prev_origin;
				match[i]  = (cur_dist==0) ? (diag_prev_match + 1) : diag_prev_match;
			} 
            else if (cost_candid_ins <= cost_candid_del) {
				// INSERTION
				cost[i]   = cost_candid_ins;
				origin[i] = origin[i-1];
				match[i]  = match[i-1];
			} 
            else {
				// DELETION
				cost[i] = cost_candid_del;
			}
            
            // printf("%d\t%d\t%d\t%d\n",j,cost[i],match[i],origin[i]);
			// remember current cell for next iteration
            diag_prev_cost   = helper_cost;
            diag_prev_match  = helper_match;
            diag_prev_origin = helper_origin;

		}
		// column finished
		//
        while ((last >= 0) && (cost[last] > k)) {
            --last;
        }
        
        if (last < m) {
            last += 1;
        }
        else if (can_stop_in_s2) {
            int length = m + min(origin[m],0);
            int final_cost  = cost[m];
            int final_match = match[m];
            if ( (length >= min_overlap) 
                && (final_cost <= length * error_rate_thres) 
                && (final_match > best_match || (final_match==best_match && final_cost<best_cost)) 
            ) {
                best_cost    = final_cost;
                best_match   = final_match;
                best_origin  = origin[m];
                best_stop_s1 = m;
                best_stop_s2 = j;
                // exact match
                if (final_cost == 0 && final_match == m) {
                    break;
                }
            }
        }
	}   // END_DP_colunm_iteration
    
    if (max_n == n) {
        int first_i = (can_stop_in_s1) ? 0 : m;
        for (i=first_i; i<=m; i++) {
            int length = i + min(origin[i],0);
            int final_cost  = cost[i];
            int final_match = match[i];
            if ( (length >= min_overlap) 
                && (final_cost <= length * error_rate_thres) 
                && (final_match > best_match || (final_match==best_match && final_cost<best_cost)) 
            ) {
                best_cost   = final_cost;
                best_match  = final_match;
                best_origin = origin[i];
                best_stop_s1 = i;
                best_stop_s2 = n;
            }
        }
    }
    
    if (best_cost == (m + n)) {
        // return negative default
        Inline_Stack_Vars;

        Inline_Stack_Reset;
        Inline_Stack_Push(sv_2mortal(newSViv(-1)));
        Inline_Stack_Push(sv_2mortal(newSViv(-1)));
        Inline_Stack_Push(sv_2mortal(newSViv(m+n)));
        Inline_Stack_Push(sv_2mortal(newSViv(0)));
        Inline_Stack_Done;
        Inline_Stack_Return(4);
    }
    
	int start_s1, start_s2;
	if (best_origin >= 0) {
		start_s1 = 0;
		start_s2 = best_origin;
	} else {
		start_s1 = -best_origin;
		start_s2 = 0;
	}

    if (best_stop_s1 - start_s1 <= 0) {
        // empty alignment
        // return negative default
        Inline_Stack_Vars;

        Inline_Stack_Reset;
        Inline_Stack_Push(sv_2mortal(newSViv(-1)));
        Inline_Stack_Push(sv_2mortal(newSViv(-1)));
        Inline_Stack_Push(sv_2mortal(newSViv(m+n)));
        Inline_Stack_Push(sv_2mortal(newSViv(0)));
        Inline_Stack_Done;
        Inline_Stack_Return(4);
    }
    else {
        // return negative default
        Inline_Stack_Vars;

        Inline_Stack_Reset;
        Inline_Stack_Push(sv_2mortal(newSViv(start_s2)));
        Inline_Stack_Push(sv_2mortal(newSViv(best_stop_s2)));
        Inline_Stack_Push(sv_2mortal(newSViv(best_cost)));
        Inline_Stack_Push(sv_2mortal(newSViv(best_match)));
        Inline_Stack_Done;
        Inline_Stack_Return(4);
    }
}



