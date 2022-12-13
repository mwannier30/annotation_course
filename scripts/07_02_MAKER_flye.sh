#!/usr/bin/env bash

#SBATCH --cpus-per-task=15
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00
#SBATCH --job-name=maker
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/data/users/mwannier/maker_%j.o
#SBATCH --error=/data/users/mwannier/maker_%j.e
#SBATCH --partition=pall

#import modules
module load SequenceAnalysis/GenePrediction/maker/2.31.9;
module load Blast/ncbi-blast/2.9.0+;

#set variables
WORKDIR=/data/users/mwannier/annotation/MAKER/flye
UNIPROT=/data/courses/assembly-annotation-course/CDS_annotation/uniprot-plant_reviewed.fasta

#move to working directory 
cd $WORKDIR

#run maker
maker

#move to output directory
cd flye.maker.output

#generate gff and fasta files
gff3_merge -d flye_master_datastore_index.log
fasta_merge -d flye_master_datastore_index.log

#Build shorter IDs/Names for MAKER genes and transcripts
maker_map_ids --prefix 'C24_' flye.all.gff > flye.all.id.map

#Map short IDs/Names to MAKER fasta and gff files
map_fasta_ids flye.all.id.map flye.all.maker.proteins.fasta
map_fasta_ids flye.all.id.map flye.all.maker.transcripts.fasta
map_gff_ids flye.all.id.map flye.all.gff

#create blast db
makeblastdb -in $UNIPROT -dbtype prot -out blast_db

##Blastp predicted proteins against the UniProt/SwissProt database
blastp -query flye.all.maker.proteins.fasta -db blast_db \
       -num_threads 10 -outfmt 6 -evalue 1e-10 > flye_blastp_output

#Map putative functions to the MAKER produced GFF3 and FASTA files:
maker_functional_fasta $UNIPROT flye_blastp_output \
        flye.all.maker.proteins.fasta > flye_functional.fasta
maker_functional_gff $UNIPROT flye_blastp_output flye.all.gff > flye_functional.gff