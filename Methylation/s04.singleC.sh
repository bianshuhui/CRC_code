meth_name=$1
ref=$2

main_dir=$3

bam_dir=$main_dir/01.bam
singleC_dir=$main_dir/02.singleC

samtools_exe= #/software_dir/samtools-0.1.18/samtools
database= #/database_dir
perl_exe= #/usr/bin/perl
pile2singleC_pl= #/bin_dir/singleC_metLevel.hg19.pl

bam_sam_dir=$bam_dir/$meth_name/genome
singleC_sam_dir=$singleC_dir/$meth_name

if [ ! -d $singleC_sam_dir ]
        then mkdir -p $singleC_sam_dir
fi

$samtools_exe mpileup -O -f $database/${ref}.fa                     \
    $bam_sam_dir/$meth_name.sort.rmdup.bam                          \
    > $singleC_sam_dir/$meth_name.pileup                         && \

$perl_exe $pile2singleC_pl                                          \
    $singleC_sam_dir/$meth_name.pileup                              \
    > $singleC_sam_dir/$meth_name.singleC                        && \



rm $singleC_sam_dir/$meth_name.pileup
gzip -f $singleC_sam_dir/$meth_name.singleC