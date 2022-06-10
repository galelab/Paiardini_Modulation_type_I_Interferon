#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -t 5800-12
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH --partition HoldingPen

#AUTHOR: Leanne Whitmore 

#NOTE using samtools 1.10 version already installed on server (at this point this is the most recent version)
#NOTE using htseq LW local installed version which is version 0.13.5

source /share/lwhitmo/pythonvenvs/py3venv/bin/activate #loads htseqv0.13.5
htseq-count -h

module load samtools/1.10
samtools --version

processors=25
mkdir "sorted_mappingRhesus"
mkdir "sorted_mappingRhesus/logs"
mkdir "counts0.13.5Rhesus/"

Al_files=()
mappingresultspath="./mapping2.7.5aRhesus/"

Al_files=("$mappingresultspath"*.out.bam)

GTF_FILE='/vol01/genome/Macaca_mulatta/10.100_STAR_2.7.5a/Macaca_mulatta.Mmul_10.100.gtf'

echo "Number of Alignment files = ${#Al_files[*]} should be 138"

counter=0
while [ "$counter" -lt  ${#Al_files[*]} ]
do
    echo "Counter $counter"
    #SAMTOOLS
    ##-o output file name (note: SBN means sorted by name)
    ##-O output file format 
    ##-n sort by name
    ##-T where to write temporary file 
    ##--threads number of threads to use 
    echo "Sorting $sample_name..."
    if [ "$counter" -lt ${#Al_files[*]} ]
    then
    	##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
    	sample_name=${Al_files[$counter]}
    	sample_name=${sample_name%Aligned.out.bam}

    	##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)}
        sample_name=${sample_name#$mappingresultspath}
        echo "$sample_name"
    	srun -c 10 samtools sort --threads 10 -o ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out.bam -O bam -n -T ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out "${Al_files[$counter]}" 1>> ./sorted_mappingRhesus/logs/"$sample_name"_mapping.log 2>&1 &
	    counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#Al_files[*]} ]
    then
        ##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
        sample_name=${Al_files[$counter]}
    	sample_name=${sample_name%Aligned.out.bam}
        ##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        sample_name=${sample_name#$mappingresultspath}
        echo "$sample_name"
	    srun -c 10 samtools sort --threads 10 -o ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out.bam -O bam -n -T ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out "${Al_files[$counter]}" 1>> ./sorted_mappingRhesus/logs/"$sample_name"_mapping.log 2>&1 &
         counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#Al_files[*]} ]
    then
        ##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
        sample_name=${Al_files[$counter]}
    	sample_name=${sample_name%Aligned.out.bam}
        ##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        sample_name=${sample_name#$mappingresultspath}
        echo "$sample_name"
        srun -c 10 samtools sort --threads 10 -o ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out.bam -O bam -n -T ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out "${Al_files[$counter]}" 1>> ./sorted_mappingRhesus/logs/"$sample_name"_mapping.log 2>&1 &
        counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#Al_files[*]} ]
    then
        ##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
        sample_name=${Al_files[$counter]}
    	sample_name=${sample_name%Aligned.out.bam}
        #2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        sample_name=${sample_name#$mappingresultspath}
        echo "$sample_name"
    	srun -c 10 samtools sort --threads 10 -o ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out.bam -O bam -n -T ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out "${Al_files[$counter]}" 1>> ./sorted_mappingRhesus/logs/"$sample_name"_mapping.log 2>&1 &
    fi
    counter=$((counter+1))
    wait 
done 

counter=0
while [ "$counter" -lt  ${#Al_files[*]} ]
do 
    ###Htseq-count PARAMETERS        
    ##--stranded - whether sequencing was done strand specific (generally this is reverse with Gale lab protocal)
    ##--format - format of input files 
    ##--mode (intersection-nonempty) - is the less strigent parameter for assigning a read to a gene 
    ##--idattr - gff/gtf attribute to be used as feature ID (defualt is gene_id for gtf and RNA seq)
    ##--type (Default for gtf files is exon so not specified here)
    ##--minqual (skip reads with alignment score of 10 (default) or less)
    
    if [ "$counter" -lt ${#Al_files[*]} ]
    then
   	 ##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
    	sample_name=${Al_files[$counter]}
    	sample_name=${sample_name%Aligned.out.bam}
    	##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        sample_name=${sample_name#$mappingresultspath}
    	echo "Running htseq-count on $sample_name..."
    	srun -c 1 htseq-count -n 1 --stranded=reverse --format=bam --mode=intersection-nonempty --idattr=gene_id ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out.bam "$GTF_FILE"> ./counts0.13.5Rhesus/"$sample_name"_counts.txt &
	counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#Al_files[*]} ]
    then
         ##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
        sample_name=${Al_files[$counter]}
    	sample_name=${sample_name%Aligned.out.bam}
        ##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        sample_name=${sample_name#$mappingresultspath}
        echo "Running htseq-count on $sample_name..."
        srun -c 1 htseq-count -n 1 --stranded=reverse --format=bam --mode=intersection-nonempty --idattr=gene_id ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out.bam "$GTF_FILE"> ./counts0.13.5Rhesus/"$sample_name"_counts.txt &
        counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#Al_files[*]} ]
    then
         ##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
        sample_name=${Al_files[$counter]}
    	sample_name=${sample_name%Aligned.out.bam}
        ##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        sample_name=${sample_name#$mappingresultspath}
        
        echo "Running htseq-count on $sample_name..."
        srun -c 1 htseq-count -n 1 --stranded=reverse --format=bam --mode=intersection-nonempty --idattr=gene_id ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out.bam "$GTF_FILE"> ./counts0.13.5Rhesus/"$sample_name"_counts.txt &
        counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#Al_files[*]} ]
    then
         ##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
        sample_name=${Al_files[$counter]}
    	sample_name=${sample_name%Aligned.out.bam}
        ##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        sample_name=${sample_name#$mappingresultspath}
        
        echo "Running htseq-count on $sample_name..."
        srun -c 1 htseq-count -n 1 --stranded=reverse --format=bam --mode=intersection-nonempty --idattr=gene_id ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out.bam "$GTF_FILE"> ./counts0.13.5Rhesus/"$sample_name"_counts.txt &
        counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#Al_files[*]} ]
    then
         ##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
        sample_name=${Al_files[$counter]}
    	sample_name=${sample_name%Aligned.out.bam}
        ##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        sample_name=${sample_name#$mappingresultspath}
        
        echo "Running htseq-count on $sample_name..."
        srun -c 1 htseq-count -n 1 --stranded=reverse --format=bam --mode=intersection-nonempty --idattr=gene_id ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out.bam "$GTF_FILE"> ./counts0.13.5Rhesus/"$sample_name"_counts.txt &
        counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#Al_files[*]} ]
    then
         ##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
        sample_name=${Al_files[$counter]}
    	sample_name=${sample_name%Aligned.out.bam}
        ##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        sample_name=${sample_name#$mappingresultspath}
        
        echo "Running htseq-count on $sample_name..."
        srun -c 1 htseq-count -n 1 --stranded=reverse --format=bam --mode=intersection-nonempty --idattr=gene_id ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out.bam "$GTF_FILE"> ./counts0.13.5Rhesus/"$sample_name"_counts.txt &
        counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#Al_files[*]} ]
    then
         ##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
        sample_name=${Al_files[$counter]}
    	sample_name=${sample_name%Aligned.out.bam}
        ##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        sample_name=${sample_name#$mappingresultspath}
        
        echo "Running htseq-count on $sample_name..."
        srun -c 1 htseq-count -n 1 --stranded=reverse --format=bam --mode=intersection-nonempty --idattr=gene_id ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out.bam "$GTF_FILE"> ./counts0.13.5Rhesus/"$sample_name"_counts.txt &
        counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#Al_files[*]} ]
    then
         ##1.Removes read and file information from file name (i.e will remove _R2_001.fastq.gz)
        sample_name=${Al_files[$counter]}
    	sample_name=${sample_name%Aligned.out.bam}
        ##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        sample_name=${sample_name#$mappingresultspath}

        echo "Running htseq-count on $sample_name..."
        srun -c 1 htseq-count -n 1 --stranded=reverse --format=bam --mode=intersection-nonempty --idattr=gene_id ./sorted_mappingRhesus/"$sample_name"Aligned.SBN.out.bam "$GTF_FILE"> ./counts0.13.5Rhesus/"$sample_name"_counts.txt &
    fi
    counter=$((counter+1))
    wait
done

mutt -s "HTSEQv0.13.5 genome v10.100 Rhesus counts done for NHP Zika " lwhitmo@uw.edu -c lwhitmo@uw.edu
