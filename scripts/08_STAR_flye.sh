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

module load UHTS/Assembler/cufflinks/2.2.1;

GFF=/data/users/mwannier/annotation/MAKER/flye/flye.maker.output/flye_functional.gff

mkdir ../genome_directory/flye; 

cd ../genome_directory/flye

gffread -E $GFF -T -o  flye_functional.gtf

GENOMEDIR=/data/users/mwannier/annotation/genome_directory/flye
READS=/data/users/mwannier/FS22_Assembly/participant_2/Illumina
FLYE=/data/users/mwannier/FS22_Assembly/polishing/flye/pilon_flye/flye.fasta

module add UHTS/Aligner/STAR/2.7.9a

zcat $READS/SRR1584462_1.fastq.gz > read1.fasta
zcat $READS/SRR1584462_2.fastq.gz > read2.fasta


STAR \
    --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir $GENOMEDIR \
    --genomeFastaFiles $FLYE \
    --sjdbGTFfile flye_functional.gtf


STAR --runThreadN 4 --genomeDir $GENOMEDIR --readFilesIn read1.fasta read2.fasta --outFileNamePrefix star --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 10 --outFilterMismatchNoverLmax 0.01 --alignIntronMax 60000 