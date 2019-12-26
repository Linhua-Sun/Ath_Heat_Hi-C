#!/bin/bash

# set several traps
set -e
set -u
set -o pipefail

#------------------------------------------------------------------------------------
#
# Aim: single-end RNA-seq data analysis pipeline (short reads, have a lot of reads with length < 50)
# Add STAR mapping | Add 
# Author: Linhua Sun (Linhua_Sun@pku.edu.cn)
# Last modified: 2017-09-23 | 2018-10-07 | 2018-11-30 | 2019-11-10
# Usage: sh script.sh CL100012044_L02_34_1.fq.gz CL100012044_L02_34_1 $(pwd)/OUT ## example Use aboso path
#
#------------------------------------------------------------------------------------

# parameters required to change (based to different machines' limitation)

THREADS="8"

#------------------------------------------------------------------------------------
##Locations require to change (based to different computers) && Assign the I/O variables

# Example fastq file:

FqFile=${1}
ID=${2}
TOTAL_OUTPUT=${3}

BIN="/data1/linhua/QIANLAB/PROJECT/hybrids/Mapping/BIN"

OUTPUT="${TOTAL_OUTPUT}/${ID}" ## Each sample should have a uniq Dir for store files

LOG="${OUTPUT}/LOG"
QC="${OUTPUT}/Quality_Control"

#------------------------------------------------------------------------------------
##mkdir new directories

if [ ! -d ${TOTAL_OUTPUT} ]
                then mkdir -p ${TOTAL_OUTPUT}
fi

if [ ! -d ${OUTPUT} ]
                then mkdir -p ${OUTPUT}
fi

if [ ! -d ${LOG} ]
                then mkdir -p ${LOG}
fi

if [ ! -d ${QC} ]
                then mkdir -p ${QC}
fi


##Quality control on raw/clean sequencing data from fastp!

echo "Use fastp remove low quality reads and remove adapters auto!" >> ${TOTAL_OUTPUT}/${ID}_general_log.txt

fastp -w ${THREADS} --length_required 30 \
	-i ${1} \
	-o ${OUTPUT}/${ID}_fastp.fq.gz \
	--json ${OUTPUT}/${ID}_fastp.json \
	--html ${OUTPUT}/${ID}_fastp.html > ${OUTPUT}/${ID}_fastp.log 2>&1

## Mapping by STAR (Extremely fast (also does splice alignment, requires at least 30 Gb memory))
## Sorted
## BAM Output
## wig file (RPM Norm)
## Intron Length

#------------------------------------------------------------------------------------
##Mapping by STAR

echo "Mapping the reads to reference genome undier the guide annotation!" >> ${TOTAL_OUTPUT}/${ID}_general_log.txt

STAR  --runMode alignReads \
    --genomeDir /data1/linhua/QIANLAB/TAIR10/TAIR10_Araport11_STAR_INDEX/STAR_Araport11  \
    --runThreadN ${THREADS}   \
    --readFilesIn ${OUTPUT}/${ID}_fastp.fq.gz \
    --outFileNamePrefix ${OUTPUT}/${ID}_STAR \
    --readFilesCommand gunzip -c \
    --outReadsUnmapped Fastx \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI NM MD AS XS \
    --outSAMstrandField intronMotif \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNmax 2  \
    --alignIntronMin 35  \
    --alignIntronMax 2000 \
    --outWigType wiggle \
    --outWigNorm RPM \
    --outWigStrand Unstranded > ${OUTPUT}/${ID}_STAR.log 2>&1

samtools index -@ ${THREADS} ${OUTPUT}/${ID}_STARAligned.sortedByCoord.out.bam &

SIZE="/data1/linhua/software/BIN/chrom.sizes"

#wigToBigWig ${OUTPUT}/${ID}_STARSignal.Unique.str1.out.wig         ${SIZE} ${OUTPUT}/${ID}_STARSignal.Unique.str1.out.bw &
#wigToBigWig ${OUTPUT}/${ID}_STARSignal.UniqueMultiple.str1.out.wig ${SIZE} ${OUTPUT}/${ID}_STARSignal.UniqueMultiple.str1.out.bw &


echo "The pipeline for SE50 RNA-seq have finished!" >> ${TOTAL_OUTPUT}/${ID}_general_log.txt

#------------------------------------------------------------------------------------
##get uniq mapped bam and bigwig file || Qualilty control on bams

echo "Convert bam file to bigwig for visulization and quantification by !" >> ${TOTAL_OUTPUT}/${ID}_general_log.txt

## Uniq map

BAM="${OUTPUT}/${ID}_STARAligned.sortedByCoord.out.bam"

samtools view -@ ${THREADS} -h ${BAM} |grep -E '^@|NH:i:1' |samtools view -@ ${THREADS} -Sb - > ${OUTPUT}/${ID}_STAR_Uniq.bam

sambamba index -t ${THREADS} ${OUTPUT}/${ID}_STAR_Uniq.bam

## convert to wigwig
bam2wig.py -s ${CHROM} -i ${OUTPUT}/${ID}_STAR_Uniq.bam  -o ${OUTPUT}/${ID}_STAR_Uniq -u

rm -rf *.wig
echo "The pipeline for SE50 RNA-seq have finished!" >> ${TOTAL_OUTPUT}/${ID}_general_log.txt