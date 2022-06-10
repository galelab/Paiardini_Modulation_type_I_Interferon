#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 68
#SBATCH -t 5800-12
#SBATCH --partition HoldingPen
#SBATCH -o slurm-%j.out


#AUTHOR: Leanne Whitmore 
# Filter globin and RNA reads
module load bowtie2/2.4.2

bowtie2 --version

BOWTIE2_LIBRARIES="./bowtie2_libraries_hmrRNAglobin/human_mouse_rhesus_05262021_rRNA_v2"
RAW_SOURCEDIR='./trimgalore_results/'
R1_files=("$RAW_SOURCEDIR"*R1*.fq.gz)
R2_files=("$RAW_SOURCEDIR"*R2*.fq.gz)

echo "Number samples = ${#R1_files[*]} R1 files"
echo "Number samples = ${#R2_files[*]} R2 files"


mkdir nohmrRNA_noglobin
mkdir logs
mkdir logs/rRNAfilter
mkdir logs/rRNA_globinfilter

counter=0

while [ "$counter" -lt  ${#R1_files[*]} ]
do

	echo "Counter variable $counter"
        
        if [ "$counter" -lt ${#R1_files[*]} ]
        then
		
		sample_name=${R1_files[$counter]}
		sample_name=${sample_name%_R*}
		sample_name=${sample_name#$RAW_SOURCEDIR}

		echo "Sample being processed $sample_name"
                echo "Read 1 file ${R1_files[$counter]}"
                echo "Read 2 file ${R2_files[$counter]}"

		srun --error ./logs/rRNA_globinfilter/"$sample_name"_rRNA_globinfilter.log -c 17 bowtie2  -p 16 -5 1 -x $BOWTIE2_LIBRARIES --un-conc-gz ./nohmrRNA_noglobin/"$sample_name"_nohmrRNA_noglobin.fastq.gz -1 ${R1_files[$counter]} -2 ${R2_files[$counter]} -S ./nohmrRNA_noglobin/"$sample_name".sam  &
		counter=$((counter+1))
	fi
        if [ "$counter" -lt ${#R1_files[*]} ]
        then
                sample_name=${R1_files[$counter]}
                sample_name=${sample_name%_R*}
                sample_name=${sample_name#$RAW_SOURCEDIR}

                echo "Sample being processed $sample_name"
                echo "Read 1 file ${R1_files[$counter]}"
                echo "Read 2 file ${R2_files[$counter]}"

                srun --error ./logs/rRNA_globinfilter/"$sample_name"_rRNA_globinfilter.log -c 17 bowtie2  -p 16 -5 1 -x $BOWTIE2_LIBRARIES --un-conc-gz ./nohmrRNA_noglobin/"$sample_name"_nohmrRNA_noglobin.fastq.gz -1 ${R1_files[$counter]} -2 ${R2_files[$counter]} -S ./nohmrRNA_noglobin/"$sample_name".sam  &
                counter=$((counter+1))
        fi
        if [ "$counter" -lt ${#R1_files[*]} ]
        then
                sample_name=${R1_files[$counter]}
                sample_name=${sample_name%_R*}
                sample_name=${sample_name#$RAW_SOURCEDIR}

                echo "Sample being processed $sample_name"
                echo "Read 1 file ${R1_files[$counter]}"
                echo "Read 2 file ${R2_files[$counter]}"

                srun --error ./logs/rRNA_globinfilter/"$sample_name"_rRNA_globinfilter.log -c 17 bowtie2  -p 16 -5 1 -x $BOWTIE2_LIBRARIES --un-conc-gz ./nohmrRNA_noglobin/"$sample_name"_nohmrRNA_noglobin.fastq.gz -1 ${R1_files[$counter]} -2 ${R2_files[$counter]} -S ./nohmrRNA_noglobin/"$sample_name".sam  &
        	counter=$((counter+1))
        fi

        if [ "$counter" -lt ${#R1_files[*]} ]
        then
                sample_name=${R1_files[$counter]}
                sample_name=${sample_name%_R*}
                sample_name=${sample_name#$RAW_SOURCEDIR}

                echo "Sample being processed $sample_name"
                echo "Read 1 file ${R1_files[$counter]}"
                echo "Read 2 file ${R2_files[$counter]}"

		srun --error ./logs/rRNA_globinfilter/"$sample_name"_rRNA_globinfilter.log -c 17 bowtie2  -p 16 -5 1 -x $BOWTIE2_LIBRARIES --un-conc-gz ./nohmrRNA_noglobin/"$sample_name"_nohmrRNA_noglobin.fastq.gz -1 ${R1_files[$counter]} -2 ${R2_files[$counter]} -S ./nohmrRNA_noglobin/"$sample_name".sam  &
        fi
	counter=$((counter+1))
	wait
done

mutt -s "rRNA and Globin removal completed CovTEN  of trim galore results" lwhitmo@uw.edu -c lwhitmo@uw.edu 

