#!/bin/bash

mkdir ./output_fig15

dataSet="ctBones"
file=../data/"$dataSet".vti
echo "    prog $dataSet"
for i in {0..8}
do
	dl=$((8-$i))
	echo $dl
	../install/bin/ttkPersistenceDiagramCmd -i $file -P 1 -S 8 -E $dl -d 4 -PP 1 -t 1 -o ./output_fig15/output_temp
	sleep 1
	mv ./output_fig15/output_temp_port\#0.vtu ./output_fig15/pd_"$dataSet"."$i".vtu
	sleep 1
done

