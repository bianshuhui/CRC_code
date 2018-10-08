meth_name=$1
trim_fq1=$2
trim_fq2=$3
ref=$4
main_dir=$5

perl_exe=/usr/bin/perl
bismark_exe=/software_dir/bismark_v0.7.6/bismark
bowtie_dir=/software_dir/bowtie-1.1.1
changeID_pl=/bin_dir/ChangeReadID.pl

trim_fq_dir=$main_dir/00.1.trim_data
bam_dir=$main_dir/01.bam

database=/date/bianshuhui/hg19/DNA_seq
samtools_exe=/datd/bianshuhui/software/samtools-0.1.18/samtools

trim_sam_dir=$trim_fq_dir/$meth_name
bam_sam_dir=$bam_dir/$meth_name/genome

sam=$bam_sam_dir/${trim_fq1}_bismark_pe.sam
unmap_fq1=$bam_sam_dir/${trim_fq1}_unmapped_reads_1.txt
unmap_fq2=$bam_sam_dir/${trim_fq2}_unmapped_reads_2.txt

sam1_unmap=$bam_sam_dir/unmap1/${trim_fq1}_unmapped_reads_1.txt_bismark
sam2_unmap=$bam_sam_dir/unmap2/${trim_fq2}_unmapped_reads_2.txt_bismark


if [ ! -d $bam_sam_dir ]
        then mkdir -p $bam_sam_dir
fi



$bismark_exe     --fastq    --non_directional   --unmapped                  \
    --phred33-quals --path_to_bowtie $bowtie_dir                            \
    --output_dir $bam_sam_dir --temp_dir $bam_sam_dir $database             \
    -1 $trim_sam_dir/$trim_fq1 -2 $trim_sam_dir/$trim_fq2                && \
    $samtools_exe view -u -b -S -t $database/${ref}.fa $sam                |\
    $samtools_exe sort -m 200000000 - $sam.sort


$bismark_exe     --fastq    --non_directional   --unmapped                  \
    --phred33-quals --path_to_bowtie $bowtie_dir                            \
    --output_dir $bam_sam_dir/unmap1 --temp_dir $bam_sam_dir/unmap1         \
    $database   $unmap_fq1                                               && \
    $samtools_exe view -uSb -t $database/${ref}.fa $sam1_unmap.sam         |\
    $samtools_exe sort -m 200000000 - $sam1_unmap.sort


$bismark_exe     --fastq    --non_directional   --unmapped                  \
    --phred33-quals --path_to_bowtie $bowtie_dir                            \
    --output_dir $bam_sam_dir/unmap2 --temp_dir $bam_sam_dir/unmap2         \
    $database   $unmap_fq2                                               && \
    $samtools_exe view -uSb -t $database/${ref}.fa $sam2_unmap.sam         |\
    $samtools_exe sort -m 200000000 - $sam2_unmap.sort


$perl_exe $changeID_pl $sam.sort.bam       $sam.sort.ReID.bam            && \
$samtools_exe rmdup    $sam.sort.ReID.bam                                   \
                       $sam.sort.ReID.rmdup.bam


$perl_exe $changeID_pl $sam1_unmap.sort.bam $sam1_unmap.sort.ReID.bam    && \
$samtools_exe rmdup -s $sam1_unmap.sort.ReID.bam                            \
                       $sam1_unmap.sort.ReID.rmdup.bam


$perl_exe $changeID_pl $sam2_unmap.sort.bam $sam2_unmap.sort.ReID.bam    && \
$samtools_exe rmdup -s $sam2_unmap.sort.ReID.bam                            \
                       $sam2_unmap.sort.ReID.rmdup.bam


# merge bam
$samtools_exe merge -f $bam_sam_dir/$meth_name.rmdup.bam                    \
                       $sam.sort.ReID.rmdup.bam                             \
                       $sam1_unmap.sort.ReID.rmdup.bam                      \
                       $sam2_unmap.sort.ReID.rmdup.bam


$samtools_exe sort -m 200000000 $bam_sam_dir/$meth_name.rmdup.bam           \
    $bam_sam_dir/$meth_name.sort.rmdup

$samtools_exe index $bam_sam_dir/$meth_name.sort.rmdup.bam

##rm  $sam1_unmap.sort.ReID.bam $sam2_unmap.sort.ReID.bam                     \
##    $sam1_unmap.sort.bam      $sam2_unmap.sort.bam                          \
##    $sam  $sam.sort.ReID.bam  $unmap_fq1   $unmap_fq2                       \
##    ${sam1_unmap/bismark/unmapped_reads.txt}                                \
##    ${sam2_unmap/bismark/unmapped_reads.txt}                                \
##    $sam1_unmap.sam   $sam2_unmap.sam                                       \
##    $sam1_unmap.sort.ReID.rmdup.bam $sam2_unmap.sort.ReID.rmdup.bam

rm $sam $unmap_fq1 $unmap_fq2 $sam.sort.ReID.bam $sam.rmdup.bam $sam.sort.ReID.rmdup.bam                                                \
$sam1_unmap.sam ${sam1_unmap/bismark/unmapped_reads.txt} $sam1_unmap.sort.bam $sam1_unmap.sort.ReID.bam $sam1_unmap.sort.ReID.rmdup.bam \
$sam2_unmap.sam ${sam2_unmap/bismark/unmapped_reads.txt} $sam2_unmap.sort.bam $sam2_unmap.sort.ReID.bam $sam2_unmap.sort.ReID.rmdup.bam


[ -e $bam_sam_dir/$meth_name.sort.rmdup.bam ] && rm $bam_sam_dir/$meth_name.rmdup.bam