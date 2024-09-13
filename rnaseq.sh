#!/bin/bash
#SBATCH --job-name=RNASEQ
#SBATCH --mail-type=ALL
#SBATCH --mail-user=o.akinsuyi@ufl.edu
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=40
#SBATCH --output=serial_test_%j.log
#SBATCH --account=microbiology-dept
#SBATCH --qos=microbiology-dept-b


module load trimmomatic

#Paired end reads
trimmomatic PE GSM7011594_1.fastq GSM7011594_2.fastq GSM7011594_1.paired.fastq GSM7011594_1.unpaired.fastq GSM7011594_2.paired.fastq GSM7011594_2.unpaired.fastq SLIDINGWINDOW:4:20 

###download the reference genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
gunzip GRCm38.primary_assembly.genome.fa.gz

module load  hisat2

#use hisat2 to read alignments. First index reference genome
hisat2-build GRCm38.primary_assembly.genome.fa reference_index

#now map reference genome to create an alignment file as sam file
hisat2 --max-intronlen 50000 -x reference_index -1 GSM7011594_1.paired.fastq -2 GSM7011594_2.paired.fastq -S alignment.sam 

module load samtools

##### Convert sam file to bam file
samtools view -bS -h alignment.sam > alignment.bam
samtools sort alignment.bam -o alignment.sorted.bam
samtools index alignment.sorted.bam

#####download the gtf file
wget  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz
gunzip gencode.vM23.primary_assembly.annotation.gtf.gz

####Measure the counts
module load subread
featureCounts -T 24 -p -t exon -g gene_name -a gencode.vM23.primary_assembly.annotation.gtf -o counts.txt alignment.sorted.bam 

