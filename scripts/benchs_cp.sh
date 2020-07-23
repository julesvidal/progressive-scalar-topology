#!/bin/sh
mkdir ./outputs_cp
while read dataSet
do
        echo $dataSet
        file=../data/"$dataSet".vti
        echo "sequential"
        for try in {0..11}
        do
                echo "try $try"
                echo "   ftm $dataSet"
                sleep 1
                ../install/bin/ttkScalarFieldCriticalPointsCmd -i $file -t 1  > ./outputs_cp/ttk_"$dataSet"_try"$try"
                sleep 1
                echo "   non prog $dataSet"
                ../install/bin/ttkScalarFieldCriticalPointsCmd -i $file -P 1 -S 0 -E 0 -d 4 -t 1 -I 0 > ./outputs_cp/nonProg_"$dataSet"_try"$try"
                echo "    prog $dataSet"
                sleep 1
                ../install/bin/ttkScalarFieldCriticalPointsCmd -i $file -P 1 -S 10 -E 0 -d 4 -t 1 -I 0 > ./outputs_cp/prog_"$dataSet"_try"$try"
        done



        echo "parallel"
        for try in {0..11}
        do
                echo "try $try"
                echo "   ttk"
                sleep 1
                ../install/bin/ttkScalarFieldCriticalPointsCmd -i $file -t 16  > ./outputs_cp/ttk_"$dataSet"_para_try"$try"
                sleep 1
                echo "   non prog"
                ../install/bin/ttkScalarFieldCriticalPointsCmd -i $file -P 1 -S 0 -E 0 -d 4 -I 0 -t 16  > ./outputs_cp/nonProg_"$dataSet"_para_try"$try"
                echo "    prog"
                sleep 1
                ../install/bin/ttkScalarFieldCriticalPointsCmd -i $file -P 1 -S 10 -E 0 -d 4 -I 0 -t 16 > ./outputs_cp/prog_"$dataSet"_para_try"$try"
        done
done < ../data/data_list
