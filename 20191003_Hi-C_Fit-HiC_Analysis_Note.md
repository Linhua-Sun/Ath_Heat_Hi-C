bash createFitHiCContacts-hic.sh [Juicer's dump command] [chr1] [chr2] [Output file name]
        
[Juicer's dump command]    Full path to the output of Juicer's dump command
[chr1]                     Chromosome 1 of the argument used in Juicer's dump command
[chr2]                     Chromosome 2 of the argument used in Juicer's dump command
[Output file name]          Name of output file



bash createFitHiCContacts-hic.sh 


Control_Chr1_NoPRs_FRAG.tsv.gz 
Control_Chr1_NoPRs_INT.tsv.gz 
Control_Chr1_NoPRs_bias.gz

## 测试FitHiC对单一染色体的效果如何？


> 注意我们自己的数据，和Chang Liu的数据依然有所差别？为什么？是测序的问题吗？(valid interactions?)

`fithic -i Control_Chr1_NoPRs_INT.tsv.gz -f Control_Chr1_NoPRs_FRAG.tsv.gz -t Control_Chr1_NoPRs_bias.gz -r 0 -U 25000 -L 2000 -r 0 -o ControlChr1NoPRULr0_Out`
fithic -i Control_Chr1_NoPRs_INT.tsv.gz -f Control_Chr1_NoPRs_FRAG.tsv.gz -t Control_Chr1_NoPRs_bias.gz -v -r 0 -U 25000 -L 2000 -r 0 -o ControlChr1NoPRULr0_v_Out > ControlChr1NoPRULr0_v_Out.log 2>&1

## ----

关键NNromd的脚本来自于`FitHiC`: HiCKRy="/data1/linhua/software/TEST/fithic/fithic/utils/HiCKRy.py"
脚本存储的位置(来自于GITHUB): `/data1/linhua/software/TEST/fithic/fithic/utils/`


#### 对全基因组互作数据进行Norm



```
python /data1/linhua/software/TEST/fithic/fithic/utils/HiCKRy.py -i Control_NoPRs_INT.tsv.gz -f Control_NoPRs_FRAG.tsv.gz -o Control_NoPRs_bias.tsv.gz
```
Control_NoPRs_INT.tsv.gz
Control_NoPRs_FRAG.tsv.gz
Control_NoPRs_bias.tsv.gz

nohup fithic -v -r 0 -U 25000 -L 2000 -r 0 -b 50 -i Control_NoPRs_INT.tsv.gz      -f Control_NoPRs_FRAG.tsv.gz      -t Control_NoPRs_bias.tsv.gz   -o AllControlNoPRULr0_v_OutBin50   > AllControlNoPRULr0_v_OutBin50.log 2>&1 &

nohup fithic -v -r 0 -U 25000 -L 2000 -r 0 -b 50 -i Control_Chr1_NoPRs_INT.tsv.gz -f Control_Chr1_NoPRs_FRAG.tsv.gz -t Control_Chr1_NoPRs_bias.gz  -o ControlChr1NoPRULr0_v_OutBin50  > ControlChr1NoPRULr0_v_OutBin50.log 2>&1 &

nohup fithic -v -r 0 -U 25000 -L 2000 -r 0 -b 50 \
-i Control_Chr1_NoPRs_INT.tsv.gz \
-f Control_Chr1_NoPRs_FRAG.tsv.gz \
-t Control_Chr1_NoPRs_bias.gz  \
-o ControlChr1NoPRULr0_v_OutBin50  \
> ControlChr1NoPRULr0_v_OutBin50.log 2>&1 &

nohup fithic -v -r 0 -U 25000 -L 2000 -r 0 -b 50 -tU 4 -tL 0 \
-i Control_Chr1_NoPRs_INT.tsv.gz \
-f Control_Chr1_NoPRs_FRAG.tsv.gz \
-t Control_Chr1_NoPRs_bias.gz  \
-o tControlChr1NoPRULr0_v_OutBin50  \
> tControlChr1NoPRULr0_v_OutBin50.log 2>&1 &

nohup fithic -tU 4 -tL 0 -v -r 0 -U 25000 -L 2000 -r 0 -b 50 -i Control_NoPRs_INT.tsv.gz      -f Control_NoPRs_FRAG.tsv.gz      -t Control_NoPRs_bias.tsv.gz   -o tAllControlNoPRULr0_v_OutBin50   > tAllControlNoPRULr0_v_OutBin50.log 2>&1 &

**为什么R version的FitHiC 和 Python3 版本的FitHiC的结果差异很大？对于Ath C忽然** 

#### 测试FitHiC on Whole Genome Datasets including cis and trans chromatin interactions datasets


## 为什么我们的数据和Chang Liu的数据Call出来的差距很大？数量是接近的，但是位置上差异很大？
## 我们用所有染色体的数据一起分析，单一染色体的是有点不一样的吗？
## 如何处理PR区域的情况？
## bias 似乎差异不是很大，如果控制下会有什么变化呢？似乎有效的Loops会多一些
## BIN到50,似乎是比较合适？Call出来的Loops会多一些？







KRy='/data1/linhua/QIANLAB/PROJECT/hic/fithic_compare/RF/test/fithic/fithic/utils/HiCKRy.py'

python ${KRy} -i Liu_Col_Chr1_NoPRs_INT.tsv.gz -f Liu_Col_Chr1_NoPRs_FRAG.tsv.gz -o Liu_Col_Chr1_NoPRs_bias_0.05.gz



GetFitHiC='/data1/linhua/QIANLAB/PROJECT/hic/Loop_fithic/TSV_Fit/GetFitHiC.R'

Rscript ${GetFitHiC} Liu_Col_Chr1_NoPRs_FRAG.tsv.gz Liu_Col_Chr1_NoPRs_INT.tsv.gz Liu_Col_Chr1_NoPRs_bias_0.05.gz Liu_Col_Chr1_NoPRs_fithic_0.05






nohup fithic -v -r 0 -U 25000 -L 2000 -r 0 -b 50 -tU 4 -tL 0 \
-i Liu_Col_Chr1_NoPRs_INT.tsv.gz \
-f Liu_Col_Chr1_NoPRs_FRAG.tsv.gz \
-t Liu_Col_Chr1_NoPRs_bias.gz \
-o Liu_Col_Chr1_NoPRs_bias0.001_vr0U25000L2000r0b50tU4tL0  \
> Liu_Col_Chr1_NoPRs_bias0.001_vr0U25000L2000r0b50tU4tL0.log 2>&1 &

nohup fithic -v -r 0 -U 25000 -L 2000 -r 0 -b 50 -tU 4 -tL 0 \
-i Liu_Col_Chr1_NoPRs_INT.tsv.gz \
-f Liu_Col_Chr1_NoPRs_FRAG.tsv.gz \
-t Liu_Col_Chr1_NoPRs_bias_0.05.gz \
-o Liu_Col_Chr1_NoPRs_bias0.05_vr0U25000L2000r0b50tU4tL0  \
> Liu_Col_Chr1_NoPRs_bias0.05_vr0U25000L2000r0b50tU4tL0.log 2>&1 &



zcat Liu_Col_Chr1_NoPRs_bias.gz	|awk '{OFS="\t"; print $1,$2,$2+1,$3}' > Liu_Col_Chr1_NoPRs_bias.bedGraph &
zcat Liu_Col_Chr1_NoPRs_bias_0.05.gz	|awk '{OFS="\t"; print $1,$2,$2+1,$3}' > Liu_Col_Chr1_NoPRs_bias_0.05.bedGraph &