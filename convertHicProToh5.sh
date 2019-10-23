#!/bin/bash

conda activate Peaks

BED1K="/data1/linhua/QIANLAB/PROJECT/hic/Cold/RES_1000.bed"

ID=$(basename ${1} .matrix)

hicConvertFormat -m ${1} --bedFileHicpro ${BED1K}  --inputFormat hicpro --outputFormat h5 -o ${ID}.h5

for i in 2 10 20 40 50 100; do echo "hicMergeMatrixBins -m ${ID}.h5 --outFileName ${ID}_${i}kb.h5 --numBins ${i} &" ;done > BatchMergeIntoLargerBins_${ID}.sh

sh BatchMergeIntoLargerBins_${ID}.sh
