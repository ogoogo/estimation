#!/bin/sh




id=21

csv_count=$((id / 10000 + 1))
echo "$csv_count"


input_csv="../../output/1115/information.csv"
echo "$input_csv"

awk 'NR == 21' $input_csv > inforow.txt
povray povray.pov +W659 +H494 Output_File_Name=../../output/1115/images/raw/21
for i in {1..41}; do
    cp ../../output/1115/images/raw/21.png "../../output/1115/images/raw/$i.png"
done
awk 'NR == 62' $input_csv > inforow.txt
povray povray.pov +W659 +H494 Output_File_Name=../../output/1115/images/raw/62
for i in {42..82}; do
    cp ../../output/1115/images/raw/62.png "../../output/1115/images/raw/$i.png"
done
awk 'NR == 103' $input_csv > inforow.txt
povray povray.pov +W659 +H494 Output_File_Name=../../output/1115/images/raw/103
for i in {83..123}; do
    cp ../../output/1115/images/raw/103.png "../../output/1115/images/raw/$i.png"
done





