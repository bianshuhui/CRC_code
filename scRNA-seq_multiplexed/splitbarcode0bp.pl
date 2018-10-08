#!/usr/bin/perl
use strict;
use warnings;
no strict 'refs';

use PerlIO::gzip;
my $usage = "perl $0 IN1.gz IN2.gz IN3.gz IN4.gz ... OUT1.gz OUT2.gz...OUTx.gz
<IN1>:fq
<IN2>:(start pos 0-base /base number before matching with regular expression)5
<IN3>:GGTCTT
<IN4>:TCTGGT
...
<OUT1>:fq
<OUT2>:fq
...
<OUTx>:fq(has no id)";


my ($n,$file,$i,$out,$j,$line1,$line2,$line3,$line4,$temp,$m,$error,$lengthb);
die $usage unless @ARGV >=3;
$n = @ARGV;
$file = ($n - 3)/2;

open IN1,"gzip -dc $ARGV[0] |" or die $!;
for($i=$file + 3 - 1;$i<=$n - 1;$i++)
        {
         $out = "OUT".$i;
         open $out,"| gzip -c > $ARGV[$i].2.fq.gz" or die $!;
         }
while (<IN1>)
        {
         $j = 0;
         $line1 = $_;
         $line2 = <IN1>;
         $line3 = <IN1>;
         $line4 = <IN1>;
         for($i=3 - 1;$i<=$file + 2 - 1;$i++)
                {
                 $lengthb = length($ARGV[$i]);
                 $m = 0;
                 $error = 0;
                 while ($m <= $lengthb - 1)
                        {
                         if (substr($ARGV[$i],$m,1) ne substr($line2,$ARGV[1]+$m,1))
                                {
                                 $error ++;
                                }
                         $m ++;
                        }
                 if ($error < 1)
                        {
                         $temp = $i + $file;
                         $out = "OUT".$temp;
                         print $out "$line1$line2$line3$line4";
                         $j = $file + 1;
                        }
                 $j ++;
                 if ($j eq $file)
                        {
                         $temp = $n - 1;
                         $out = "OUT".$temp;
                         print $out "$line1$line2$line3$line4";
                        }
                }

        }
for($i=$file + 3 - 1;$i<=$n - 1;$i++)
        {
         $out = "OUT".$i;
         close $out;
        }
