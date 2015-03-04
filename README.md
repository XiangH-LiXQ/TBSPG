# TBSPG

TBSPG automatically trim the reads, run the multi/noshare/uni mapping, sort, and  deduplicate the alignment map, then count the number of reads (readsno.txt) and coverage (cov.txt), using entire or half length reads.

- The TBSPG pipeline supplementary package is available on https://github.com/XiangH-LiXQ/TBSPG.git, which contains six shell scripts and two folders.
- Six shell scripts written for the TBSPG pipelines are:
-- TBSPG_multimap_pair.sh:  A pipeline for conducting multi-mapping (multimap) of paired-end reads without removing the reads shared between compartments (nuclear, chloroplast, and mitochondrial).
-- TBSPG_multimap_noshare_pair.sh:  A pipeline for multi-mapping (multimap) of paired-end reads with removal of the reads shared between compartments.
-- TBSPG_unimap_pair.sh:  A pipeline for unique mapping of paired-end reads.
-- TBSPG_multimap_single.sh:  A pipeline for multi-mapping (multimap) of single reads without removal of the reads shared between compartments.
-- TBSPG_multimap_noshare_single.sh:  A pipeline for multi-mapping of single reads with removal of the reads shared between compartments.
-- TBSPG_unimap_single.sh:  A pipeline for unique mapping for single reads.
- programs: a folder containing all PERL scripts in the bin folder, and five publicly available third-party programs, where BWA and Samtools need to be installed, and GATK, Picard and Trimmomatic need the Java platform.
- example_input: a folder containing the example input files, i.e., Reads1.fq, Reads2.fq. (Note that the Reads1.fq file contains 1000 reads extracted from the paired-end R1 sequence file from Illumina HighSeq2000 sequencing of the potato root DNA D3356. The Reads2.fq file contains 1000 reads extracted from the paired-end R2 sequence file from the same sequencing of D3356.)
