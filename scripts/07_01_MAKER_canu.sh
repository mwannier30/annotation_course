#!/usr/bin/env bash

#SBATCH --cpus-per-task=15
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00
#SBATCH --job-name=maker_canu
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/data/users/mwannier/maker_%j.o
#SBATCH --error=/data/users/mwannier/maker_%j.e
#SBATCH --partition=pall

#import modules 
module load SequenceAnalysis/GenePrediction/maker/2.31.9;
module load Blast/ncbi-blast/2.9.0+;

#set variables
WORKDIR=/data/users/mwannier/annotation/MAKER/canu
UNIPROT=/data/courses/assembly-annotation-course/CDS_annotation/uniprot-plant_reviewed.fasta

#move to working directory
cd $WORKDIR

#run maker
maker

#move to output directory
cd canu.maker.output

#Generate gff and fasta files
gff3_merge -d canu_master_datastore_index.log
fasta_merge -d canu_master_datastore_index.log

#Build shorter IDs/Names for MAKER genes and transcripts
maker_map_ids --prefix 'C24_' canu.all.gff > canu.all.id.map

#Map short IDs/Names to MAKER fasta and gff files
map_fasta_ids canu.all.id.map canu.all.maker.proteins.fasta
map_fasta_ids canu.all.id.map canu.all.maker.transcripts.fasta
map_gff_ids canu.all.id.map canu.all.gff

#create blast db
makeblastdb -in $UNIPROT -dbtype prot -out blast_db

#Blastp predicted proteins against the UniProt/SwissProt database
blastp -query canu.all.maker.proteins.fasta -db blast_db \
        -num_threads 10 -outfmt 6 -evalue 1e-10 > canu_blastp_output

#Map putative functions to the MAKER produced GFF3 and FASTA files:
maker_functional_fasta $UNIPROT canu_blastp_output \
        canu.all.maker.proteins.fasta > canu_functional.fasta
maker_functional_gff $UNIPROT canu_blastp_output canu.all.gff > canu_functional.gff