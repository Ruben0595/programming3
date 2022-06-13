#!/bin/bash

export OUTPUT=/students/2021-2022/master/Ruben_DSLS/output
export FILE1=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R1_001_BC24EVACXX.filt.fastq
export FILE2=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R2_001_BC24EVACXX.filt.fastq

mkdir -p  ${OUTPUT}
mkdir -p output
#seq 29 2 31 | parallel -j60 "velveth $OUTPUT/{} {} -longPaired -fastq $FILE1 $FILE2 && velvetg $OUTPUT/{}" | (python3 N50.py)

#cat $OUTPUT/21/contigs.fa | python3 N50.py


#echo "$( cat ${OUTPUT}/23/contigs.fa)"
#parallel echo "$( cat ${OUTPUT}/{} /contigs.fa)" ::: {21..23..2}  | python3 N50.py
#seq 21 2 23 | ${OUTPUT}{} | python3 N50.py
seq 25 2 31 | parallel -j16 'velveth $OUTPUT/{} {} -longPaired -fastq $FILE1 $FILE2 && velvetg $OUTPUT/{} && cat $OUTPUT/{}/contigs.fa | (python3 N50.py && echo -e Kmer_size:{}; ) >> output/output.csv'
# Python script to clean the data
python3 fileclean.py
# remove redundant files
rm -rf $OUTPUT