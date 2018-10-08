in_fq1=$1
in_fq2=$2
meth_name=$3
main_dir=$4

trim_dir=$main_dir/00.1.trim_data


perl_exe=/usr/bin/perl
bin_dir=/bin_dir

trim_sam_dir=$trim_dir/$meth_name

$perl_exe $bin_dir/QC_lite_0.0.7c.pl                                        \
    --quality 20 --length 50 --left_R1 9 --left_R2 9 --trimx 1              \
    --symmetry  --graph                                                     \
    --outdir $trim_sam_dir                                                  \
    $in_fq1 $in_fq2
