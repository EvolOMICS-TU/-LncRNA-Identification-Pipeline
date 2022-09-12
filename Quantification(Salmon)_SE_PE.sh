#!/bin/bash

#Quantification by salmon

#Building index 

mkdir -p /home/cluster/nath/quant_test
mkdir -p /home/cluster/nath/quant_test/fastp
mkdir -p /home/cluster/nath/quant_test/fastqc

transcripts=/home/cluster/nath/pipeline/transcriptome/*.fa
fastp=/home/cluster/nath/TEST/results/fastp
base=/home/cluster/nath/quant_pneumonia
deseq=/home/cluster/nath/pipeline/deseq/*.R
fastp=/home/cluster/nath/quant_test/fastp
fastq=/home/cluster/nath/quant_test/fastq
reads=/home/cluster/nath/pneumonia_dataset/
threads="50"

cd $reads

find $reads -name "*.fastq.gz" | sort | paste - - | while read A B

do
a=`basename ${A} | sed 's/.sra_1/_1/' | awk -F "." '{print $1}'`
b=`basename ${B} | sed 's/.sra_2/_2/' | awk -F "." '{print $1}'`

echo ""
echo [`date +"%Y-%m-%d %H:%M:%S"`] "TRIMMING FILES"

fastp --thread=16 --length_required=10 --qualified_quality_phred=32 --in1=${A} --in2=${B} --out1=$a\_trimmed.fastq.gz --out2=$b\_trimmed.fastq.gz --json=$a.json --html=$a.html

mv -v $a\_trimmed.fastq.gz $b\_trimmed.fastq.gz $fastp

echo ""
echo "--------------------------------$a AND $b TRIMMED--------------------------------"

done
mv *.html "$fastp"  
mv *.json "$fastp" 

echo [`date +"%Y-%m-%d %H:%M:%S"`] "REPORTS GENERATED"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING FASTQC"

cd $fastp
fastqc *.gz -t $threads
mv *.html "$fastqc"  
mv *.zip "$fastqc" 

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QC REPORTS GENERATED"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QUALITY CHECK AND PROCESSING DONE"

cd $base

salmon index -t $transcripts -i index -k 31

#Quantification

find $fastp -name "*_trimmed.fastq.gz" | sort | paste - - | while read A B

do

a=`basename ${A} | awk -F "." '{print $1}' | awk -F "_" '{print $1}'`

salmon quant -i index -l A -1 ${A} -2 ${B} --validateMappings -o $a

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QUANTIFICATION DONE!!!!"

: << 'END'

#Tximport and DESeq

cp $deseq $base

R -f tximport_deseq.R

mv -v res.csv bak_res.csv
echo -n "", > res.csv; cat bak_res.csv >> res.csv #fixes the left shift of column names
rm bak_res.csv
 
END


