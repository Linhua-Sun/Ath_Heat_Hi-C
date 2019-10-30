#!/bin/bash


source activate 

conda activate /data1/linhua/software/anaconda2/envs/hic


ID=$(basename ${1} .h5)

hicConvertFormat -m ${ID}.h5 --inputFormat h5 --outputFormat cool --outFileName ${ID}.cool
cooler dump -t pixels --join ${ID}.cool > ${ID}_pixels_join.tsv
cooler dump -t pixels ${ID}.cool > ${ID}_pixels.tsv
paste ${ID}_pixels_join.tsv ${ID}_pixels.tsv |awk '{OFS="\t"; print $1,$8,$4,$9,$7}' > ${ID}.HIC
rm -rf ${ID}.cool ${ID}_pixels_join.tsv ${ID}_pixels.tsv

gzip ${ID}.HIC
