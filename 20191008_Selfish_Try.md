## selfish test at 20 kb
```
selfish -m1 /data1/linhua/QIANLAB/PROJECT/hic/HeatProjectHi-C-Results-ToBeSubmitted/Control_20000.matrix \
		-m2 /data1/linhua/QIANLAB/PROJECT/hic/HeatProjectHi-C-Results-ToBeSubmitted/Heat_20000.matrix \
	-ch Chr1 \
	-r 20kb -o ./output.npy
```

```
for i in 1 2 3 4 5
do
	selfish \
	-m1 Control_20000.matrix \
	-m2 Heat_20000.matrix \
	-ch ${i} \
	-r 20kb \
	-o Chr${i}CtrlVsHeat20kSelfish.npy \
	-bed1 TAIR10_BIN20000.bed \
	-bed2 TAIR10_BIN20000.bed
done
```

