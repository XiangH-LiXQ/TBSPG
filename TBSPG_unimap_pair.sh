print_help () 
{
echo "This script automatically trim the reads, run the uniquely mapping, sort, and  deduplicate the alignment map, then count the number of reads (readsno.txt) and coverage (cov.txt), using entire or half length reads."
echo "-h help"
echo "-p 'E' = entire reads analysis, 'H' = half reads analysis"
echo "-i input untrimmed read1 file.fastq"
echo "-j input untrimmed read2 file.fastq"
echo "-r input reference genome file.fasta"
echo "-o output file DIRECTORY"
echo "Example: sh TBSPG_unimap_pair.sh   -p E/H   -i READ1.fastq   -j READ2.fastq   -r REF_GENOME.fasta   -o OUT"
}

if [ $# -lt 5 ]; then
print_help
exit 1
fi

while getopts "h:p:i:j:o:r:" opt
do
case "$opt" in
h) print_help; exit 1
;;
p) half=$OPTARG
;;
i) in=$OPTARG
;;
j) inn=$OPTARG
;;
r) ref=$OPTARG
;;
o) outputdir=$OPTARG
;;
\?) print_help; exit 1
;;
	esac
done


#### Formating of reference genome. If the genome has been formated, please close them. ####
bwa index -a bwtsw $ref
samtools faidx $ref
file=${ref}
java -jar programs/picard-tools-1.101/CreateSequenceDictionary.jar R=$ref O=${file%.*}.dict


mkdir ${outputdir}


#### Anlysis of entire reads ####
if [ $half = E ]; then

## trimmomatic
java -jar programs/Trimmomatic-0.30/trimmomatic-0.30.jar PE -phred33 $in $inn ${outputdir}/entire_reads1_trimmed.fq.gz ${outputdir}/entire_reads1_unpaired.fq.gz ${outputdir}/entire_reads2_trimmed.fq.gz ${outputdir}/entire_reads2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
gunzip ${outputdir}/entire_reads1_trimmed.fq.gz
gunzip ${outputdir}/entire_reads2_trimmed.fq.gz
rm -r ${outputdir}/entire_reads1_unpaired.fq.gz ${outputdir}/entire_reads2_unpaired.fq.gz

## mapping
bwa aln $ref ${outputdir}/entire_reads1_trimmed.fq >${outputdir}/entire_reads1_trimmed.sai
bwa aln $ref ${outputdir}/entire_reads2_trimmed.fq >${outputdir}/entire_reads2_trimmed.sai
bwa sampe -r "@RG\tID:IDa\tSM:SM\tPL:Illumina" -P $ref ${outputdir}/entire_reads1_trimmed.sai ${outputdir}/entire_reads2_trimmed.sai ${outputdir}/entire_reads1_trimmed.fq ${outputdir}/entire_reads2_trimmed.fq >${outputdir}/entire_out.sam
rm -r ${outputdir}/entire_reads1_trimmed.fq ${outputdir}/entire_reads1_trimmed.sai ${outputdir}/entire_reads2_trimmed.fq ${outputdir}/entire_reads2_trimmed.sai

## unimap
samtools view -h -F 4 -q 1 -S ${outputdir}/entire_out.sam >${outputdir}/entire_multimap.sam
more ${outputdir}/entire_multimap.sam |perl -ne 'chomp;my @ID=split/\t/;if(($ID[0]=~/\@SQ|\@RG|\@PG/)||($ID[12] eq "XT:A:U")){print "$_\n"}' >${outputdir}/entire_unimap.sam
samtools view -b -S ${outputdir}/entire_unimap.sam >${outputdir}/entire_unimap.bam
rm -r ${outputdir}/entire_multimap.sam ${outputdir}/entire_out.sam ${outputdir}/entire_unimap.sam

## sort bam
java -jar programs/picard-tools-1.101/SortSam.jar I=${outputdir}/entire_unimap.bam O=${outputdir}/entire_unimap_sort.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT
samtools index ${outputdir}/entire_unimap_sort.bam
rm -r ${outputdir}/entire_unimap.bam

## dedup bam
java -jar programs/picard-tools-1.101/MarkDuplicates.jar ASSUME_SORTED=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT INPUT=${outputdir}/entire_unimap_sort.bam OUTPUT=${outputdir}/entire_unimap_sort_dedup.bam METRICS_FILE=${outputdir}/entire_unimap_sort_dedup.metrics
rm -r ${outputdir}/entire_unimap_sort.bam ${outputdir}/entire_unimap_sort.bam.bai

## count reads no
samtools view ${outputdir}/entire_unimap_sort_dedup.bam >${outputdir}/entire_unimap_sort_dedup.sam
perl programs/bin/calc_reads.pl ${outputdir}/entire_unimap_sort_dedup.sam >${outputdir}/entire_unimap_sort_dedup_readsno.txt
rm -r ${outputdir}/entire_unimap_sort_dedup.sam

## count coverage
java -jar programs/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar -T DepthOfCoverage --validation_strictness SILENT -omitIntervals -omitLocusTable -ct 0 -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -ct 25 -ct 30 -allowPotentiallyMisencodedQuals -R $ref -I ${outputdir}/entire_unimap_sort_dedup.bam -o ${outputdir}/entire_unimap_sort_dedup_cov
perl programs/bin/sum_cov.pl ${outputdir}/entire_unimap_sort_dedup_cov >${outputdir}/entire_unimap_sort_dedup_cov_sum.txt

## integrate the results of reads_no and coverage
cat ${outputdir}/entire_unimap_sort_dedup_readsno.txt ${outputdir}/entire_unimap_sort_dedup_cov_sum.txt |perl programs/bin/write_resluts.pl | perl programs/bin/calc_cov.pl >${outputdir}/entire_unimap_sort_dedup_readsno_cov.xls


#### Anlysis of half reads ####
elif [ $half = H ]; then

## get half-length reads
perl programs/bin/get_half.pl $in >${outputdir}/half_reads1.fq
perl programs/bin/get_half.pl $inn >${outputdir}/half_reads2.fq

## trimmomatic
java -jar programs/Trimmomatic-0.30/trimmomatic-0.30.jar PE -phred33 ${outputdir}/half_reads1.fq ${outputdir}/half_reads2.fq ${outputdir}/half_reads1_trimmed.fq.gz ${outputdir}/half_reads1_unpaired.fq.gz ${outputdir}/half_reads2_trimmed.fq.gz ${outputdir}/half_reads2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
gunzip ${outputdir}/half_reads1_trimmed.fq.gz
gunzip ${outputdir}/half_reads2_trimmed.fq.gz
rm -r ${outputdir}/half_reads1_unpaired.fq.gz ${outputdir}/half_reads2_unpaired.fq.gz

## mapping
bwa aln $ref ${outputdir}/half_reads1_trimmed.fq >${outputdir}/half_reads1_trimmed.sai
bwa aln $ref ${outputdir}/half_reads2_trimmed.fq >${outputdir}/half_reads2_trimmed.sai
bwa sampe -r "@RG\tID:IDa\tSM:SM\tPL:Illumina" -P $ref ${outputdir}/half_reads1_trimmed.sai ${outputdir}/half_reads2_trimmed.sai ${outputdir}/half_reads1_trimmed.fq ${outputdir}/half_reads2_trimmed.fq >${outputdir}/half_out.sam
rm -r ${outputdir}/half_reads1_trimmed.fq ${outputdir}/half_reads1_trimmed.sai ${outputdir}/half_reads2_trimmed.fq ${outputdir}/half_reads2_trimmed.sai

## unimap
samtools view -h -F 4 -q 1 -S ${outputdir}/half_out.sam >${outputdir}/half_multimap.sam
more ${outputdir}/half_multimap.sam |perl -ne 'chomp;my @ID=split/\t/;if(($ID[0]=~/\@SQ|\@RG|\@PG/)||($ID[12] eq "XT:A:U")){print "$_\n"}' >${outputdir}/half_unimap.sam
samtools view -b -S ${outputdir}/half_unimap.sam >${outputdir}/half_unimap.bam
rm -r ${outputdir}/half_multimap.sam ${outputdir}/half_out.sam ${outputdir}/half_unimap.sam

## sort bam
java -jar programs/picard-tools-1.101/SortSam.jar I=${outputdir}/half_unimap.bam O=${outputdir}/half_unimap_sort.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT
samtools index ${outputdir}/half_unimap_sort.bam
rm -r ${outputdir}/half_unimap.bam

## dedup bam
java -jar programs/picard-tools-1.101/MarkDuplicates.jar ASSUME_SORTED=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT INPUT=${outputdir}/half_unimap_sort.bam OUTPUT=${outputdir}/half_unimap_sort_dedup.bam METRICS_FILE=${outputdir}/half_unimap_sort_dedup.metrics
rm -r ${outputdir}/half_unimap_sort.bam ${outputdir}/half_unimap_sort.bam.bai

## count reads no
samtools view ${outputdir}/half_unimap_sort_dedup.bam >${outputdir}/half_unimap_sort_dedup.sam
perl programs/bin/calc_reads.pl ${outputdir}/half_unimap_sort_dedup.sam >${outputdir}/half_unimap_sort_dedup_readsno.txt
rm -r ${outputdir}/half_unimap_sort_dedup.sam

## count coverage
java -jar programs/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar -T DepthOfCoverage --validation_strictness SILENT -omitIntervals -omitLocusTable -ct 0 -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -ct 25 -ct 30 -allowPotentiallyMisencodedQuals -R $ref -I ${outputdir}/half_unimap_sort_dedup.bam -o ${outputdir}/half_unimap_sort_dedup_cov
perl programs/bin/sum_cov.pl ${outputdir}/half_unimap_sort_dedup_cov >${outputdir}/half_unimap_sort_dedup_cov_sum.txt

## integrate the results of reads_no and coverage
cat ${outputdir}/half_unimap_sort_dedup_readsno.txt ${outputdir}/half_unimap_sort_dedup_cov_sum.txt |perl programs/bin/write_resluts.pl | perl programs/bin/calc_cov.pl >${outputdir}/half_unimap_sort_dedup_readsno_cov.xls


fi;
