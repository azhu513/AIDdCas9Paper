#!/bin/bash

# First base directory path
base_dir1="scratch/VEP_annotated"

# First array of file names
files=("Sample_01_annotated.txt.modified" "Sample_02_annotated.txt.modified" "Sample_03_annotated.txt.modified" "Sample_04_annotated.txt.modified" "Sample_05_annotated.txt.modified" "Sample_06_annotated.txt.modified" "Sample_07_annotated.txt.modified" "Sample_08_annotated.txt.modified" "Sample_09_annotated.txt.modified" "Sample_10_annotated.txt.modified" "Sample_11_annotated.txt.modified" "Sample_12_annotated.txt.modified" "Sample_13_annotated.txt.modified" "Sample_14_annotated.txt.modified" "Sample_15_annotated.txt.modified" "Sample_16_annotated.txt.modified")

# Concatenate all sample 05 to sample 10 to one big file
cat "${base_dir1}/${files[5]}" "${base_dir1}/${files[6]}" "${base_dir1}/${files[7]}" "${base_dir1}/${files[8]}" "${base_dir1}/${files[9]}" "${base_dir1}/${files[10]}" > "${base_dir1}/All_DOX_Samples.txt.modified"

files1="All_DOX_Samples.txt.modified"

# Second base directory path
base_dir2="scratch/data"

# Second array of file names
files2=("sample_01.txt" "sample_02.txt" "sample_03.txt" "sample_04.txt" "sample_05.txt" "sample_06.txt" "sample_07.txt" "sample_08.txt" "sample_09.txt" "sample_10.txt" "sample_11.txt" "sample_12.txt" "sample_13.txt" "sample_14.txt" "sample_15.txt" "sample_16.txt")

# Iterate over each file
for file in "${files2[@]}"
do
# Construct the full file path
file_path="${base_dir2}/${file}"
ref_file_path="${base_dir1}/${files1}"

echo "Processing $file_path..."

# Remove first 1 line of file 1
tail -n +2 "$file_path" > "${base_dir2}/temp.txt"

# sort both files by key columns
sort -t $'\t' -k2,2 -k3,3 "${base_dir2}/temp.txt" > "${base_dir2}/sorted_file1.txt"
sort -t $'\t' -k2,2 -k3,3 "$ref_file_path" > "${base_dir2}/sorted_file2.txt"

# Merge two sorted files and include specific columns in the outcome
awk 'BEGIN {FS=OFS="\t"} NR==FNR{key=$2"\t"$3; chrompos[key]=$4; gene[key]=$5; type[key]=$6; protein[key]=$7; cons[key]=$8; position[key]=$9; next} {print $0, chrompos[$2"\t"$3], gene[$2"\t"$3], type[$2"\t"$3], protein[$2"\t"$3], cons[$2"\t"$3], position[$2"\t"$3]}' "${base_dir2}/sorted_file2.txt" "${base_dir2}/sorted_file1.txt" > "$file_path.annotated"

rm "${base_dir2}/temp.txt" "${base_dir2}/sorted_file1.txt" "${base_dir2}/sorted_file2.txt"
echo "Done processing $file_path."
done

echo "All files processed."
