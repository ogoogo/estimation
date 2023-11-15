#!/bin/sh




id=1

csv_count=$((id / 10000 + 1))
echo "$csv_count"

# for j in `seq $csv_count 1`
# do

#     input_csv="../../output/1102/informations/output_${j}.csv"
#     echo "$input_csv"

#     tail -n +$((id % 10000)) "$input_csv" | while IFS= read -r line; do
   
#         echo "$line" > "inforow.txt"
#         povray povray.pov +W659 +H494 Output_File_Name=../../output/1102/images/raw/$id
#         echo "ID{$id}の画像生成"
#         id=$((id + 1))

#     done



# done

input_csv="../../output/5/information.csv"
echo "$input_csv"

tail -n +$((id % 10000)) "$input_csv" | while IFS= read -r line; do

    echo "$line" > "inforow.txt"
    povray povray.pov +W659 +H494 Output_File_Name=../../output/5/images/raw/$id
    echo "ID{$id}の画像生成"
    id=$((id + 1))

done



