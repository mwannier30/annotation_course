#!/usr/bin/env bash

#SBATCH --cpus-per-task=15
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=star
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/data/users/mwannier/star_%j.o
#SBATCH --error=/data/users/mwannier/star_%j.e
#SBATCH --partition=pall

#import modules
module load UHTS/Assembler/cufflinks/2.2.1;
module add UHTS/Aligner/STAR/2.7.9a

#set variables
GFF=/data/users/mwannier/annotation/MAKER/flye/flye.maker.output/flye_functional.gff
GENOMEDIR=/data/users/mwannier/annotation/genome_directory/flye
READS=/data/users/mwannier/FS22_Assembly/participant_2/Illumina
FLYE=/data/users/mwannier/FS22_Assembly/polishing/flye/pilon_flye/flye.fasta

#move to working directory
cd ../genome_directory/flye

#convert gff to gtf
gffread -E $GFF -T -o  flye_functional.gtf

#unzip illumina reads
zcat $READS/SRR1584462_1.fastq.gz > read1.fasta
zcat $READS/SRR1584462_2.fastq.gz > read2.fasta

#generate genome index
STAR \
    --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir $GENOMEDIR \
    --genomeFastaFiles $FLYE \
    --sjdbGTFfile flye_functional.gtf

#run star
STAR --runThreadN 4 --genomeDir $GENOMEDIR --readFilesIn read1.fasta read2.fasta --outFileNamePrefix star --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 10 --outFilterMismatchNoverLmax 0.01 --alignIntronMax 60000 