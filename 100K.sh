#!/bin/bash

## Aim:
## Date:
## Linhua Sun

conda activate hic

# 100Kb analysis

## 1. hicNormalize

hicNormalize -m  Control_100K.h5 Control_Rep1_100K.h5 Control_Rep2_100K.h5 Heat_100K.h5 Heat_Rep1_100K.h5 Heat_Rep2_100K.h5  --normalize smallest -o  Control_100K_Sm.h5 Control_Rep1_100K_Sm.h5 Control_Rep2_100K_Sm.h5 Heat_100K_Sm.h5 Heat_Rep1_100K_Sm.h5 Heat_Rep2_100K_Sm.h5

## 2. hicCorrectMatrix
hicCorrectMatrix correct --correctionMethod KR -m Control_100K_Sm.h5 -o Control_100K_SmKr.h5
hicCorrectMatrix correct --correctionMethod KR -m Control_Rep1_100K_Sm.h5 -o Control_Rep1_100K_SmKr.h5
hicCorrectMatrix correct --correctionMethod KR -m Control_Rep2_100K_Sm.h5 -o Control_Rep2_100K_SmKr.h5
hicCorrectMatrix correct --correctionMethod KR -m Heat_100K_Sm.h5 -o Heat_100K_SmKr.h5
hicCorrectMatrix correct --correctionMethod KR -m Heat_Rep1_100K_Sm.h5 -o Heat_Rep1_100K_SmKr.h5
hicCorrectMatrix correct --correctionMethod KR -m Heat_Rep2_100K_Sm.h5 -o Heat_Rep2_100K_SmKr.h5

## 3. hicAdjustMatrix by mask

Peri="Peri.bed"
hicAdjustMatrix --regions ${Peri} --action mask -m Control_100K_SmKr.h5  --o Control_100K_SmKrArms.h5
hicAdjustMatrix --regions ${Peri} --action mask -m Control_Rep1_100K_SmKr.h5  --o Control_Rep1_100K_SmKrArms.h5
hicAdjustMatrix --regions ${Peri} --action mask -m Control_Rep2_100K_SmKr.h5  --o Control_Rep2_100K_SmKrArms.h5
hicAdjustMatrix --regions ${Peri} --action mask -m Heat_100K_SmKr.h5  --o Heat_100K_SmKrArms.h5
hicAdjustMatrix --regions ${Peri} --action mask -m Heat_Rep1_100K_SmKr.h5  --o Heat_Rep1_100K_SmKrArms.h5
hicAdjustMatrix --regions ${Peri} --action mask -m Heat_Rep2_100K_SmKr.h5  --o Heat_Rep2_100K_SmKrArms.h5
Arms="HandArms.bed"
hicAdjustMatrix --regions ${Arms} --action mask -m Control_100K_SmKr.h5 -o Control_100K_SmKrPeri.h5 
hicAdjustMatrix --regions ${Arms} --action mask -m Control_Rep1_100K_SmKr.h5 -o Control_Rep1_100K_SmKrPeri.h5 
hicAdjustMatrix --regions ${Arms} --action mask -m Control_Rep2_100K_SmKr.h5 -o Control_Rep2_100K_SmKrPeri.h5 
hicAdjustMatrix --regions ${Arms} --action mask -m Heat_100K_SmKr.h5 -o Heat_100K_SmKrPeri.h5 
hicAdjustMatrix --regions ${Arms} --action mask -m Heat_Rep1_100K_SmKr.h5 -o Heat_Rep1_100K_SmKrPeri.h5 
hicAdjustMatrix --regions ${Arms} --action mask -m Heat_Rep2_100K_SmKr.h5 -o Heat_Rep2_100K_SmKrPeri.h5 

## 4. hicPlotMatrix 

for i in Control_100K_SmKr.h5 Control_Rep1_100K_SmKr.h5 Control_Rep2_100K_SmKr.h5 Heat_100K_SmKr.h5 Heat_Rep1_100K_SmKr.h5 Heat_Rep2_100K_SmKr.h5 Control_100K_SmKrArms.h5 Control_Rep1_100K_SmKrArms.h5 Control_Rep2_100K_SmKrArms.h5 Heat_100K_SmKrArms.h5 Heat_Rep1_100K_SmKrArms.h5 Heat_Rep2_100K_SmKrArms.h5 Control_100K_SmKrPeri.h5 Control_Rep1_100K_SmKrPeri.h5 Control_Rep2_100K_SmKrPeri.h5 Heat_100K_SmKrPeri.h5 Heat_Rep1_100K_SmKrPeri.h5 Heat_Rep2_100K_SmKrPeri.h5
do 
	echo ${i}
	hicPlotMatrix --log1p --colorMap Reds -m ${i} -o $(basename ${i} .h5).pdf -t $(basename ${i} .h5) --dpi 600
	hicPlotMatrix --log1p --colorMap Reds -m ${i} -o $(basename ${i} .h5)PerChr.pdf -t $(basename ${i} .h5) --dpi 600 --perChromosome
done

## 5 PlotDistVsCounts

### 5.1 HiCPlotDistVsCounts by Chromosome

Samples=("Control_100K_SmKr.h5" "Control_Rep1_100K_SmKr.h5" "Control_Rep2_100K_SmKr.h5" "Heat_100K_SmKr.h5" "Heat_Rep1_100K_SmKr.h5" "Heat_Rep2_100K_SmKr.h5")
Labels=("Control" "Control_Rep1" "Control_Rep2" "Heat" "Heat_Rep1" "Heat_Rep2")

NoScaleHiCPlotDistVsCounts.py -m ${Samples[@]} --labels ${Labels[@]} --maxdepth 10000000 --plotsize 6 4 --skipDiagonal -o KR_10MToSmKr100K_AllChr.pdf --outFileData KR_10MToSmKr100K_AllChr.txt
NoScaleHiCPlotDistVsCounts.py -m ${Samples[@]} --labels ${Labels[@]} --maxdepth 10000000 --plotsize 6 4 --skipDiagonal -o KR_10MToSmKr100K_PerChr.pdf --outFileData KR_10MToSmKr100K_PerChr.txt --perchr

PlotDistVsCounts.R -i KR_10MToSmKr100K_AllChr.txt -o KR_10MToSmKr100K_AllChr
PlotDistVsCounts.R -i KR_10MToSmKr100K_PerChr.txt -o KR_10MToSmKr100K_PerChr -s NN

### 5.2 HiCPlotDistVsCounts by Arms

Samples=("Control_100K_SmKrArms.h5" "Control_Rep1_100K_SmKrArms.h5" "Control_Rep2_100K_SmKrArms.h5" "Heat_100K_SmKrArms.h5" "Heat_Rep1_100K_SmKrArms.h5" "Heat_Rep2_100K_SmKrArms.h5")
Labels=("Control" "Control_Rep1" "Control_Rep2" "Heat" "Heat_Rep1" "Heat_Rep2")

NoScaleHiCPlotDistVsCounts.py -m ${Samples[@]} --labels ${Labels[@]} --maxdepth 5000000 --plotsize 6 4 --skipDiagonal -o ArmsKR_5MToSmKr100K_AllChr.pdf --outFileData ArmsKR_5MToSmKr100K_AllChr.txt
NoScaleHiCPlotDistVsCounts.py -m ${Samples[@]} --labels ${Labels[@]} --maxdepth 5000000 --plotsize 6 4 --skipDiagonal -o ArmsKR_5MToSmKr100K_PerChr.pdf --outFileData ArmsKR_5MToSmKr100K_PerChr.txt --perchr

PlotDistVsCounts.R -i ArmsKR_5MToSmKr100K_AllChr.txt -o ArmsKR_5MToSmKr100K_AllChr
PlotDistVsCounts.R -i ArmsKR_5MToSmKr100K_PerChr.txt -o ArmsKR_5MToSmKr100K_PerChr -s NN

### 5.3 HiCPlotDistVsCounts by Peri

Samples=("Control_100K_SmKrPeri.h5" "Control_Rep1_100K_SmKrPeri.h5" "Control_Rep2_100K_SmKrPeri.h5" "Heat_100K_SmKrPeri.h5" "Heat_Rep1_100K_SmKrPeri.h5" "Heat_Rep2_100K_SmKrPeri.h5")
Labels=("Control" "Control_Rep1" "Control_Rep2" "Heat" "Heat_Rep1" "Heat_Rep2")

NoScaleHiCPlotDistVsCounts.py -m ${Samples[@]} --labels ${Labels[@]} --maxdepth 1000000 --plotsize 6 4 --skipDiagonal -o PeriKR_1MToSmKr100K_AllChr.pdf --outFileData PeriKR_1MToSmKr100K_AllChr.txt
NoScaleHiCPlotDistVsCounts.py -m ${Samples[@]} --labels ${Labels[@]} --maxdepth 1000000 --plotsize 6 4 --skipDiagonal -o PeriKR_1MToSmKr100K_PerChr.pdf --outFileData PeriKR_1MToSmKr100K_PerChr.txt --perchr

PlotDistVsCounts.R -i PeriKR_1MToSmKr100K_AllChr.txt -o PeriKR_1MToSmKr100K_AllChr
PlotDistVsCounts.R -i PeriKR_1MToSmKr100K_PerChr.txt -o PeriKR_1MToSmKr100K_PerChr -s NN

