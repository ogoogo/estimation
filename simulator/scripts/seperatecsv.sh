#!/bin/zsh

# 入力CSVファイル
input_csv="../../output/1102/information1102.csv"

# 出力ディレクトリ
output_dir="../../output/1102/informations/"

# 1ファイルあたりの行数
lines_per_file=10000

# 出力ファイルのカウンタ
file_counter=0

# 出力ディレクトリが存在しない場合は作成
mkdir -p "$output_dir"

# CSVファイルを1万行ずつのファイルに分割
awk -v lines_per_file="$lines_per_file" -v output_dir="$output_dir" '
    BEGIN { 
        file_counter = 0
        file_name = output_dir "output_" file_counter ".csv"
        print > file_name
    }
    {
        if (NR % lines_per_file == 1) {
            close(file_name)
            file_counter++
            file_name = output_dir "output_" file_counter ".csv"
            print > file_name
        } else {
            print > file_name
        }
    }
' "$input_csv"

echo "分割が完了しました。"
