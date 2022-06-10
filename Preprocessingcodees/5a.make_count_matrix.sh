# !/bin/bash
###

#AUTHOR: Leanne Whitmore 
#IMPORTANT NOTE: this makes a count matrix from reads aligned by starv2.7.5a with genome 10.99 and reads counted 
#by htseq 0.13.5 
#First move files over to final folder 3.counts_results0.13.5/ two count folders because of issue with naming convention 
#and in the first pass did not get all S samples
echo 'starting run'

COUNT_SOURCEDIR='./counts0.13.5Rhesus/'

array=("$COUNT_SOURCEDIR"*.txt)

echo "Number of count files = ${#array[*]} should be 138"
#echo "Number of samples in target file ${#target_array[*]}"

# read in names from first file (they all should be the same)
cat ${array[0]} | awk '{print $1}' > count_matrix_Rhesus.txt

endoffile='_nohmrRNA_noglobin_counts.txt'
SAMPLEARRAY=()

for item in ${array[@]} #Using target list to not include duplicates of the ABL1 study and non infected samples
do
    # read in second column of each file and add it is a new column
    echo $item
    ##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
    file_name="${item}"
    ##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
    sample_name=${file_name#$COUNT_SOURCEDIR}
    sample_name=${sample_name%$endoffile}
    ##Fill array with sample names without path 
    
    SAMPLEARRAY+=($sample_name)
    eval "cat '$file_name' | awk '{print \$2}' | paste count_matrix_Rhesus.txt - > output.txt"

    mv output.txt count_matrix_Rhesus.txt
done

# delete last 5 lines
head -n -5 count_matrix_Rhesus.txt > tmp.txt && mv tmp.txt count_matrix_Rhesus.txt

# add header
output=$(printf "\t%s" "${SAMPLEARRAY[@]}")
echo -e "Name${output}" | cat - count_matrix_Rhesus.txt > tmp.txt && mv tmp.txt count_matrix_Rhesus.txt
