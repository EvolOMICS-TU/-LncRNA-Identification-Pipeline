# LncRNA-identification-pipeline
This is under development pipeline for analysis of LncRNA from RNA-Seq data by team EVoLOMICS.

The files are as follows:
1. Alpha_one_PE(Identification) and Alpha_one_SE(Identification) are LncRNA identification scripts from RNA-Seq data. 
These bash scripts accept raw fastq pair-end or single-end files as input and generate a fasta file containing putative novel LncRNA sequences.
2. Alpha_two_PE(Identification + Quantification) and Alpha_two_SE(Identification + Quantification) are combined identification and quantification
scripts that run the identification pipeline as well as quantification pipelines on transcripts. They generate raw count data of transcripts as 
output file.
3. Quantification(Salmon) is a transcript quantification script that accepts a transcriptome and fastq reads to do a pseudoalignment mapping to 
the transcriptome to generate transcript counts. This script accepts raw fastq files, processes them through fastp and runs salmon to generate raw
read counts of transcripts.
4. tximport_DESeq.R is an R script that takes the output of salmon (quant.sf) and runs it through an offset(tximport) and then processes it through
DESeq2 to generate differentially expressed transcripts which are then filtered to get DELncRNAs.
5. 200ntfilter.pl is a perl script that works on FASTA files to filter out sequences having more than 200 nucleotides.
