print_help () 
{
echo "This script automatically trim the reads, run the multi-mapping, remove the shared reads among Mitochondria, Chloroplast and Chromosomes, and sort, deduplicate the alignment map, then count the number of reads (readsno.txt) and coverage (cov.txt), using entire or half length single reads."
echo "-h help"
echo "-p 'E' = entire reads analysis, 'H' = half reads analysis"
echo "-i input untrimmed read file.fastq"
echo "-r input reference genome file.fasta"
echo "-l input a file with the name of reference Chromosome"
echo "-m input a file with the name of reference Mitochondria"
echo "-n input a file with the name of reference Chloroplast"
echo "-o output file DIRECTORY"
echo "Example: sh TBSPG_multimap_noshare_single.sh   -p E/H   -i READ.fastq   -l Chromosome_name  -m  Mitochondria_name    -n  Chloroplast_name   -r REF_GENOME.fasta   -o OUT"
}

if [ $# -lt 7 ]; then
print_help
exit 1
fi

while getopts "h:p:i:o:r:l:n:m:" opt
do
case "$opt" in
h) print_help; exit 1
;;
p) half=$OPTARG
;;
i) in=$OPTARG
;;
r) ref=$OPTARG
;;
o) outputdir=$OPTARG
;;
l) chr=$OPTARG
;;
m) mit=$OPTARG
;;
n) chl=$OPTARG
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
java -jar programs/Trimmomatic-0.30/trimmomatic-0.30.jar SE -phred33 $in ${outputdir}/entire_reads_trimmed.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
gunzip ${outputdir}/entire_reads_trimmed.fq.gz

## mapping
bwa aln $ref ${outputdir}/entire_reads_trimmed.fq >${outputdir}/entire_reads_trimmed.sai
bwa samse -r "@RG\tID:IDa\tSM:SM\tPL:Illumina" $ref ${outputdir}/entire_reads_trimmed.sai ${outputdir}/entire_reads_trimmed.fq >${outputdir}/entire_out.sam
rm -r ${outputdir}/entire_reads_trimmed.fq ${outputdir}/entire_reads_trimmed.sai

## multimap
samtools view -h -F 4 -S ${outputdir}/entire_out.sam >${outputdir}/entire_multimap.sam
rm -r ${outputdir}/entire_out.sam

## remove share
perl programs/bin/remove_share.pl $chr $mit $chl ${outputdir}/entire_multimap.sam ${outputdir}/entire_multimap_shared.sam ${outputdir}/entire_multimap_noshare.sam
samtools view -b -S ${outputdir}/entire_multimap_noshare.sam >${outputdir}/entire_multimap_noshare.bam
rm -r ${outputdir}/entire_multimap_noshare.sam ${outputdir}/entire_multimap.sam

## sort bam
java -jar programs/picard-tools-1.101/SortSam.jar I=${outputdir}/entire_multimap_noshare.bam O=${outputdir}/entire_multimap_noshare_sort.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT
samtools index ${outputdir}/entire_multimap_noshare_sort.bam
rm -r ${outputdir}/entire_multimap_noshare.bam

## dedup bam
java -jar programs/picard-tools-1.101/MarkDuplicates.jar ASSUME_SORTED=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT INPUT=${outputdir}/entire_multimap_noshare_sort.bam OUTPUT=${outputdir}/entire_multimap_noshare_sort_dedup.bam METRICS_FILE=${outputdir}/entire_multimap_noshare_sort_dedup.metrics
rm -r ${outputdir}/entire_multimap_noshare_sort.bam ${outputdir}/entire_multimap_noshare_sort.bam.bai

## count reads no
samtools view ${outputdir}/entire_multimap_noshare_sort_dedup.bam >${outputdir}/entire_multimap_noshare_sort_dedup.sam
perl programs/bin/calc_reads.pl ${outputdir}/entire_multimap_noshare_sort_dedup.sam >${outputdir}/entire_multimap_noshare_sort_dedup_readsno.txt
rm -r ${outputdir}/entire_multimap_noshare_sort_dedup.sam

## count coverage
java -jar programs/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar -T DepthOfCoverage --validation_strictness SILENT -omitIntervals -omitLocusTable -ct 0 -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -ct 25 -ct 30 -allowPotentiallyMisencodedQuals -R $ref -I ${outputdir}/entire_multimap_noshare_sort_dedup.bam -o ${outputdir}/entire_multimap_noshare_sort_dedup_cov
perl programs/bin/sum_cov.pl ${outputdir}/entire_multimap_noshare_sort_dedup_cov >${outputdir}/entire_multimap_noshare_sort_dedup_cov_sum.txt

## integrate the results of reads_no and coverage
cat ${outputdir}/entire_multimap_noshare_sort_dedup_readsno.txt ${outputdir}/entire_multimap_noshare_sort_dedup_cov_sum.txt |perl programs/bin/write_resluts.pl | perl programs/bin/calc_cov.pl >${outputdir}/entire_multimap_noshare_sort_dedup_readsno_cov.xls


#### Anlysis of half reads ####
elif [ $half = H ]; then

## get half-length reads
perl programs/bin/get_half.pl $in >${outputdir}/half_reads.fq

## trimmomatic
java -jar programs/Trimmomatic-0.30/trimmomatic-0.30.jar SE -phred33 ${outputdir}/half_reads.fq ${outputdir}/half_reads_trimmed.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
gunzip ${outputdir}/half_reads_trimmed.fq.gz

## mapping
bwa aln $ref ${outputdir}/half_reads_trimmed.fq >${outputdir}/half_reads_trimmed.sai
bwa samse -r "@RG\tID:IDa\tSM:SM\tPL:Illumina" $ref ${outputdir}/half_reads_trimmed.sai ${outputdir}/half_reads_trimmed.fq >${outputdir}/half_out.sam
rm -r ${outputdir}/half_reads_trimmed.fq ${outputdir}/half_reads_trimmed.sai

## multimap
samtools view -h -F 4 -S ${outputdir}/half_out.sam >${outputdir}/half_multimap.sam
rm -r ${outputdir}/half_out.sam

## remove share
perl programs/bin/remove_share.pl $chr $mit $chl ${outputdir}/half_multimap.sam ${outputdir}/half_multimap_shared.sam ${outputdir}/half_multimap_noshare.sam
samtools view -b -S ${outputdir}/half_multimap_noshare.sam >${outputdir}/half_multimap_noshare.bam
rm -r ${outputdir}/half_multimap_noshare.sam ${outputdir}/half_multimap.sam

## sort bam
java -jar programs/picard-tools-1.101/SortSam.jar I=${outputdir}/half_multimap_noshare.bam O=${outputdir}/half_multimap_noshare_sort.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT
samtools index ${outputdir}/half_multimap_noshare_sort.bam
rm -r ${outputdir}/half_multimap_noshare.bam

## dedup bam
java -jar programs/picard-tools-1.101/MarkDuplicates.jar ASSUME_SORTED=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT INPUT=${outputdir}/half_multimap_noshare_sort.bam OUTPUT=${outputdir}/half_multimap_noshare_sort_dedup.bam METRICS_FILE=${outputdir}/half_multimap_noshare_sort_dedup.metrics
rm -r ${outputdir}/half_multimap_noshare_sort.bam ${outputdir}/half_multimap_noshare_sort.bam.bai

## count reads no
samtools view ${outputdir}/half_multimap_noshare_sort_dedup.bam >${outputdir}/half_multimap_noshare_sort_dedup.sam
perl programs/bin/calc_reads.pl ${outputdir}/half_multimap_noshare_sort_dedup.sam >${outputdir}/half_multimap_noshare_sort_dedup_readsno.txt
rm -r ${outputdir}/half_multimap_noshare_sort_dedup.sam

## count coverage
java -jar programs/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar -T DepthOfCoverage --validation_strictness SILENT -omitIntervals -omitLocusTable -ct 0 -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -ct 25 -ct 30 -allowPotentiallyMisencodedQuals -R $ref -I ${outputdir}/half_multimap_noshare_sort_dedup.bam -o ${outputdir}/half_multimap_noshare_sort_dedup_cov
perl programs/bin/sum_cov.pl ${outputdir}/half_multimap_noshare_sort_dedup_cov >${outputdir}/half_multimap_noshare_sort_dedup_cov_sum.txt

## integrate the results of reads_no and coverage
cat ${outputdir}/half_multimap_noshare_sort_dedup_readsno.txt ${outputdir}/half_multimap_noshare_sort_dedup_cov_sum.txt |perl programs/bin/write_resluts.pl | perl programs/bin/calc_cov.pl >${outputdir}/half_multimap_noshare_sort_dedup_readsno_cov.xls


fi;
