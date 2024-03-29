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

#import module
module load UHTS/Assembler/cufflinks/2.2.1;
module add UHTS/Aligner/STAR/2.7.9a

#set variables
GFF=/data/users/mwannier/annotation/MAKER/canu/canu.maker.output/canu_functional.gff
GENOMEDIR=/data/users/mwannier/annotation/genome_directory/canu
READS=/data/users/mwannier/FS22_Assembly/participant_2/RNAseq
CANU=/data/users/mwannier/FS22_Assembly/polishing/canu/pilon_canon/canu.fasta

# move to working directory 
cd ../genome_directory/canu

#convert gff to gtf
gffread -E $GFF -T -o  canu_functional.gtf

#unzip illumina reads
zcat $READS/SRR1584462_1.fastq.gz > read1.fasta
zcat $READS/SRR1584462_2.fastq.gz > read2.fasta

#generate genome index
STAR \
    --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir $GENOMEDIR \
    --genomeFastaFiles $CANU \
    --sjdbGTFfile canu_functional.gtf \
    --genomeSAindexNbases 12

#run star
STAR --runThreadN 4 --genomeDir $GENOMEDIR --readFilesIn read1.fasta read2.fasta --outFileNamePrefix star --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 10 --outFilterMismatchNoverLmax 0.01 --alignIntronMax 60000 