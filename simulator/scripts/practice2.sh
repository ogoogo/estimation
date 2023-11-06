#!/bin/sh

input_csv="./informations/information1004.csv"
# /Applications/MATLAB_R2022b.app/bin/matlab -nodesktop -nosplash -r 'build; exit'
#  /usr/local/MATLAB/R2023a/bin/matlab -nodesktop -nosplash -r 'build; exit'


id=65439

csv_count = $((id / 10000 + 1))

for j in `seq $csv_count $61`
do

    input_csv = "./informations/output/output_${j}.csv"


    tail -n +$((id % 10000 + 1)) "$input_csv" | while IFS= read -r line; do
   
        echo "$line" > "inforow.txt"
        povray povray.pov +W659 +H494 Output_File_Name=./output/$id
        echo "ID{$id}の画像生成"
        id=$((id + 1))

    done



done



