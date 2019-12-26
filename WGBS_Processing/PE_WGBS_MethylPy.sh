#!/bin/bash
# set several traps
set -e
set -u
set -o pipefail

#------------------------------------------------------------------------------------
#
# Author: Linhua Sun (Linhua_Sun@pku.edu.cn)
# Last modified: 2018-04-16 | 2019-02-27
# Add audo QC
# Usage: sh script.sh 1h-2_clean_r1.fq.gz 1h-2_clean_r2.fq.gz 1h-2_clean ./
# Usage: sh script.sh samplename.R1.fq.gz samplename.R2.fq.gz SAMPLENAME OUTPUTDIR
#
#------------------------------------------------------------------------------------

# parameters required to change (based to different machines' limitation)

THREADS="20"

#------------------------------------------------------------------------------------
##Locations require to change (based to different computers) && Assign the I/O variables

FastqFile1=${1} 	# Anhui-Qianshan_united_R1.fq.gz
FastqFile2=${2} 	# Anhui-Qianshan_united_R2.fq.gz
SAMPLE_NAME=${3} 	# Anhui-Qianshan_united
TOTAL_OUTPUT=${4} 	# /data4/maliang/Linhua_Epi_Project/Methy/China_Methy/MethyOut

# Examples
# Anhui-Qianshan_united_R1.fq.gz
# Anhui-Qianshan_united_R2.fq.gz
# Anhui-Qianshan_united
# /data4/maliang/Linhua_Epi_Project/Methy/China_Methy/MethyOut

OUTPUT="${TOTAL_OUTPUT}/${SAMPLE_NAME}" ## Each sample should have a uniq Dir for store files
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

#------------------------------------------------------------------------------------ % Need to rewrite, not good enough %
##Quality control on raw/clean sequencing data from the company by fastqc and trim_galore

echo "Calculate the basic fastq info!" >> ${TOTAL_OUTPUT}/${SAMPLE_NAME}_general_log.txt

seqkit stats ${FastqFile1} > ${LOG}/${SAMPLE_NAME}_1_seq_stats.txt &
seqkit stats ${FastqFile2} > ${LOG}/${SAMPLE_NAME}_2_seq_stats.txt &

echo "Fastqc proform quality control on clean fastq file!" >> ${TOTAL_OUTPUT}/${SAMPLE_NAME}_general_log.txt

fastqc -t ${THREADS} -o ${QC} ${FastqFile1} ${FastqFile2} &

echo "Auto trim fastq file!" >> ${TOTAL_OUTPUT}/${SAMPLE_NAME}_general_log.txt

## trim_galore --length 70 --paired  --output_dir ${OUTPUT}

fastp --in1 ${FastqFile1} --in2 ${FastqFile2} --out1 ${OUTPUT}/${SAMPLE_NAME}_R1_val_1.fq.gz --out2 ${OUTPUT}/${SAMPLE_NAME}_R2_val_2.fq.gz -j ${OUTPUT}/${SAMPLE_NAME}.json -h ${OUTPUT}/${SAMPLE_NAME}.html --thread 8 --length_required 70

CFq_gz1=${OUTPUT}/${SAMPLE_NAME}_R1_val_1.fq.gz
CFq_gz2=${OUTPUT}/${SAMPLE_NAME}_R2_val_2.fq.gz

echo "Fastqc proform quality control on trimmmed fastq file!" >> ${TOTAL_OUTPUT}/${SAMPLE_NAME}_general_log.txt
fastqc -t ${THREADS} -o ${QC} ${CFq_gz1} ${CFq_gz2} &

#------------------------------------------------------------------------------------
## Clean Fastq files && Mapping by bismark

echo "methylpy mapping!" >> ${TOTAL_OUTPUT}/${SAMPLE_NAME}_general_log.txt

seqkit stats ${CFq_gz1} > ${LOG}/${SAMPLE_NAME}_ult_1_seq_stats.txt &
seqkit stats ${CFq_gz2} > ${LOG}/${SAMPLE_NAME}_ult_2_seq_stats.txt &

Genome="/data1/linhua/QIANLAB/PROJECT/ddccm/methylpy/example_data/genome"

## Chromosome name start with chr.
## Index Arabidopsis ----
#	methylpy build-reference \
#	--input-files genome/tair10.fa \
#	--output-prefix genome/tair10 \
#	--aligner bowtie2 \
#	--num-procs 10

methylpy paired-end-pipeline \
	--read1-files ${CFq_gz1} \
	--read2-files ${CFq_gz2} \
	--sample ${SAMPLE_NAME} \
	--forward-ref ${Genome}/tair10_f \
	--reverse-ref ${Genome}/tair10_r \
	--ref-fasta   ${Genome}/tair10.fa \
	--path-to-output ${OUTPUT} \
	--num-procs 30 \
	--sort-mem 50G \
	--remove-clonal True \
	--path-to-picard /data1/linhua/software/ \
	--unmethylated-control chrC: \
	--binom-test True \
	> ${LOG}/${SAMPLE_NAME}_methylpy_MAP.log 2>&1

for i in CG CHG CHH CNN
do
	echo "methylpy allc-to-bigwig ${i}!" >> ${TOTAL_OUTPUT}/${SAMPLE_NAME}_general_log.txt
	methylpy allc-to-bigwig \
		--mc-type ${i} \
		--allc-file ${OUTPUT}/allc_${SAMPLE_NAME}.tsv.gz \
		--output-file ${OUTPUT}/m${i}_${SAMPLE_NAME}.bw \
		--ref-fasta ${Genome}/tair10.fa \
		--bin-size 1 \
		--min-site-cov 3 \
		--add-chr-prefix True &
done

echo "MethylPy start!" >> ${TOTAL_OUTPUT}/${SAMPLE_NAME}_general_log.txt 





## QC Started ----

echo "WGBS QC start!" >> ${TOTAL_OUTPUT}/${SAMPLE_NAME}_general_log.txt
#######################################################################################################
cd ${OUTPUT}

## Generate bigwig and DSS format files
allcToDSSBigWig.R allc_${SAMPLE_NAME}.tsv.gz 5

samtools index -@ ${THREADS} ${OUTPUT}/${SAMPLE_NAME}_processed_reads_no_clonal.bam
mosdepth --threads ${THREADS} ${SAMPLE_NAME} ${OUTPUT}/${SAMPLE_NAME}_processed_reads_no_clonal.bam
samtools idxstats ${OUTPUT}/${SAMPLE_NAME}_processed_reads_no_clonal.bam > ${SAMPLE_NAME}_idxstats.txt

## Quality control extract

echo "Collecting all stat info" >> ${TOTAL_OUTPUT}/${SAMPLE_NAME}_general_log.txt
PyBIN="/data1/linhua/QIANLAB/PROJECT/BS/BIN"

cat ${LOG}/*seq_stats*|cat|grep -v "file" > ${SAMPLE_NAME}_fastq_stats.txt

python ${PyBIN}/get_fastq_info_from_seqkit.py ${SAMPLE_NAME}_fastq_stats.txt Clean > ${TOTAL_OUTPUT}/${SAMPLE_NAME}_QC.log


python ${PyBIN}/get_bam_info_from_idxstats.py ${SAMPLE_NAME}_idxstats.txt Clean >>${TOTAL_OUTPUT}/${SAMPLE_NAME}_QC.log

## Extract infos from methylpy

python ${PyBIN}/get_info_from_methylpy_PE.py  ${LOG}/${SAMPLE_NAME}_methylpy_MAP.log >> ${TOTAL_OUTPUT}/${SAMPLE_NAME}_QC.log

Rscript ${PyBIN}/get_info_from_mosdepth.R ${SAMPLE_NAME}.mosdepth.global.dist.txt  >>${TOTAL_OUTPUT}/${SAMPLE_NAME}_QC.log

#######################################################################################################
echo "WGBS QC finished!" >> ${TOTAL_OUTPUT}/${SAMPLE_NAME}_general_log.txt
