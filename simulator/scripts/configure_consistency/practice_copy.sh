#!/bin/sh


a=0
echo "hi"
matlab -r example 
# example
# while [ $a -lt $1 ]
# do
#     /Applications/MATLAB_R2022b.app/bin/matlab -nodesktop -nosplash -r 'example2; exit'

#     while read line
#     do
#     id=`echo $line | cut -d , -f 1`
#     echo $id
#     povray povray.pov +W659 +H494 Output_File_Name=./output/$id
#     echo "ID{$id}の画像生成"
#     done < inforow.txt
    
#     let a++

#     echo $a
    
# done