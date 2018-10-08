meth_name=$1
trim_fq1=$2
trim_fq2=$3
ref=$4
main_dir=$5

bismark_exe= #/software_dir/bismark_v0.7.6/bismark
bowtie_dir= #/software_dir/bowtie-1.1.1
changeID_pl= #/bin_dir/ChangeReadID.pl

trim_fq_dir=$main_dir/00.1.trim_data
bam_dir=$main_dir/01.bam

database=/date/bianshuhui/hg19/lambda_DNA/bis_conv_lambda_DNA
samtools_exe=/datd/bianshuhui/software/samtools-0.1.18/samtools

trim_sam_dir=$trim_fq_dir/$meth_name
bam_sam_dir=$bam_dir/$meth_name/lambda

if [ ! -d $bam_sam_dir ]
        then mkdir -p $bam_sam_dir
fi

sam=$bam_sam_dir/${trim_fq1}_bismark_pe.sam
unmap_fq1=$bam_sam_dir/${trim_fq1}_unmapped_reads_1.txt
unmap_fq2=$bam_sam_dir/${trim_fq2}_unmapped_reads_2.txt

sam1_unmap=$bam_sam_dir/unmap1/${trim_fq1}_unmapped_reads_1.txt_bismark
sam2_unmap=$bam_sam_dir/unmap2/${trim_fq2}_unmapped_reads_2.txt_bismark


$bismark_exe     --fastq    --non_directional                               \
    --phred33-quals --path_to_bowtie $bowtie_dir                            \
    --output_dir $bam_sam_dir --temp_dir $bam_sam_dir $database             \
    -1 $trim_sam_dir/$trim_fq1 -2 $trim_sam_dir/$trim_fq2
    
    $samtools_exe view -u -b -S -t $database/${ref}.fa $sam                |\
    $samtools_exe sort -m 200000000 - $sam.sort
