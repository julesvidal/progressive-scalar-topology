#!/bin/sh
mkdir ./outputs_pd
while read dataSet
do
        # dataSet=$line
        file=../data/"$dataSet".vti
        echo $dataSet
        echo "sequential"
        for try in {0..11}
        do
                echo "try $try"
                echo "   ftm $dataSet"
                sleep 1
                ../install/bin/ttkFTMTreeCmd -i $file -t 1 -T 0  > ./outputs_pd/ftm_"$dataSet"_try"$try"
                sleep 1
                ../install/bin/ttkFTMTreeCmd -i $file -t 1 -T 1  >> ./outputs_pd/ftm_"$dataSet"_try"$try"
                echo "    prog $dataSet"
                sleep 1
                ../install/bin/ttkPersistenceDiagramCmd -i $file -P 1 -S 8 -E 0 -d 4 -PP 1 -t 1 > ./outputs_pd/prog_"$dataSet"_try"$try"
                sleep 1
                echo "   non prog $dataSet"
                ../install/bin/ttkPersistenceDiagramCmd -i $file -P 1 -S 0 -E 0 -d 4 -PP 1 -t 1  > ./outputs_pd/nonProg_"$dataSet"_try"$try"
        done



        echo "parallel"
        for try in {0..11}
        do
                echo "try $try"
                echo "   ftm"
                sleep 1
                ../install/bin/ttkFTMTreeCmd -i $file -t 16 -T 0  > ./outputs_pd/ftm_"$dataSet"_para_try"$try"
                sleep 1
                ../install/bin/ttkFTMTreeCmd -i $file -t 16 -T 1  >> ./outputs_pd/ftm_"$dataSet"_para_try"$try"
                sleep 1
                echo "   non prog"
                ../install/bin/ttkPersistenceDiagramCmd -i $file -P 1 -S 0 -E 0 -d 4 -PP 1 -t 16  > ./outputs_pd/nonProg_"$dataSet"_para_try"$try"
                echo "    prog"
                sleep 1
                ../install/bin/ttkPersistenceDiagramCmd -i $file -P 1 -S 10 -E 0 -d 4 -PP 1 -t 16 > ./outputs_pd/prog_"$dataSet"_para_try"$try"
        done

done < ../data/data_list

