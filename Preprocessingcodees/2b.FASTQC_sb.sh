#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 5800-12
#SBATCH --partition HoldingPen
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
PATH_FASTQ="/vol08/ngs/P51/CoVTEN/CoVTEN-01/WhitmoreAnalysis/nohmrRNA_noglobin/"
samples=("$PATH_FASTQ"*.fastq.*.gz)

mkdir 'FASTQC_nohmrRNA_noglobin'

for sample in  ${samples[*]}
do
	srun -c 16 /vol01/ngs_tools/FastQC/fastqc "$sample" --noextract -t 16 -o /vol08/ngs/P51/CoVTEN/CoVTEN-01/WhitmoreAnalysis/FASTQC_nohmrRNA_noglobin/
	wait
done
