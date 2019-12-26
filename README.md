# Bioinformatic analysis code of Heat Stress Project
Bioinformatic analysis of RNA-Seq, WGBS, and Hi-C analysis in *Arabidopsis thaliana*

Date: 2019-12-26

## RNA-Seq (TE expression analysis)

## WGBS analysis (TE methylation analaysis)

## Hi-C analysis and TE expression
### 1. Hi-C Seqencing data mapping and contact matrix construction
### 2. Hi-C data visuzalition and IDE analysis
### 3. KEE analysis
### 4. PC1 of Hi-C and intergration with TE expression
### 5. Local level chromatin organization analysis with TE expression 

## Overall structure

```
|-- Hi-C_Processing
|   |-- BatchCompareHi-C_100K.sh
|   |-- BatchPlotPGT.sh
|   |-- CompareHeatAndControl_CompartmentStrength_SaddlePlot.R
|   |-- CompareHi-C_WithPreviousPublishedHi-C_ByFPC.R
|   |-- ComparePC1BedGraph_WithABCompartmentChangeInfo.R
|   |-- CompareTwoHiCmatrixV2BasedOnHiCdat.R
|   |-- CrossCompare_Two_Hi-C_Rep12.R
|   |-- Heat-activated_TEs_PlusPC1Info.R
|   |-- HiC_Heatmap_Generator.R
|   |-- Hi-C_ScalingPlotGeneratorFromHiC-ProUsingHiCE_001kb.sh
|   |-- Hi-C_ScalingPlotGeneratorFromHiC-ProUsingHiCE_100kb.sh
|   |-- KEExKEEyInteractionsPersp3D.R
|   |-- NoScaleHiCPlotDistVsCounts.py
|   |-- PlotDistVsCounts.R
|   |-- PlotLocalRescaled2DMetaChromatinOrganizationOverHeatActivatedTEs.R
|   |-- RepHIC.ini
|   |-- TEKEE34_MinMax.ini
|   |-- theme_linhua.R
|   `-- UpdatedPlotDistVsCounts.R
|-- NotUsed.tgz
|-- Others
|   |-- HeatStressNuclearSizeDistributionCompare.R
|   `-- PlotKEEsDistancesDistributionMBPVersion20191028.R
|-- README.md
|-- RNA-Seq_Processing
|   |-- 01_Heat_polyA-RNA-Seq_DifferentialAnalysis_PlusExpressionInfo_2019.R
|   |-- 02_Heat_polyA-RNA-Seq_PCG_Scatterplot.R
|   |-- 03_Heat_polyA-RNA-Seq_TE_Chromosome-wide-Plot.R
|   |-- 04_Heat_polyA-RNA-Seq_TE_Expression_Dynamics_Analysis.R
|   |-- 05_Heat_polyA-RNA-Seq_TE_Combine_Expression_Distribution.R
|   |-- 4VS4_RNA-seq-polyA-FeatureCounts-Batch-DEG-Analysis-TE-Version.R
|   |-- EpiFeaturesEnrcihmentAnalaysisHeatmap.R
|   |-- Heat_Activated_TEs_Detailed_Analysis-Enrichment.R
|   |-- Heat_Activated_TEs_MainEnrichment.R
|   |-- SingleGroupTEBasicFeaturesEnrichmentAnalysis_20191119.R
|   `-- ss_SE_US_RNA-Seq_Ath.sh
`-- WGBS_Processing
    `-- PE_WGBS_MethylPy.sh
```
