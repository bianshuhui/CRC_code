meth_name=$1
trim_fq1_list=$2
trim_fq2_list=$3

main_dir=$4

bam_dir= #$main_dir/01.bam
samtools_exe= #/software_dir/samtools-0.1.18/samtools
bam_sam_dir= #$bam_dir/$meth_name/genome

pe_bam_list=''
sam1_unmap_list=''
sam2_unmap_list=''

for trim_fq1 in `echo $trim_fq1_list | cut -d',' -f1- --output-delimiter " "`
do
sam=$bam_sam_dir/${trim_fq1}_bismark_pe.sam
pe_bam_list="${pe_bam_list} $sam.sort.ReID.rmdup.bam"

sam1_unmap=$bam_sam_dir/unmap1/${trim_fq1}_unmapped_reads_1.txt_bismark
sam1_unmap_list="${sam1_unmap_list} $sam1_unmap.sort.ReID.rmdup.bam"
done

for trim_fq2 in `echo $trim_fq2_list | cut -d',' -f1- --output-delimiter " "`
do
sam2_unmap=$bam_sam_dir/unmap2/${trim_fq2}_unmapped_reads_2.txt_bismark
sam2_unmap_list="${sam2_unmap_list} $sam2_unmap.sort.ReID.rmdup.bam"
done

$samtools_exe merge -f $bam_sam_dir/$meth_name.rmdup.bam                    \
                       $pe_bam_list                                         \
                       $sam1_unmap_list                                     \
                       $sam2_unmap_list                                  && \


$samtools_exe sort -m 200000000  $bam_sam_dir/$meth_name.rmdup.bam          \
    $bam_sam_dir/$meth_name.sort.rmdup                                   && \

$samtools_exe index $bam_sam_dir/$meth_name.sort.rmdup.bam


[ -e $bam_sam_dir/$meth_name.sort.rmdup.bam ] && rm $bam_sam_dir/$meth_name.rmdup.bam