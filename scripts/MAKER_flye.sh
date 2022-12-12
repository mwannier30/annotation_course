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

module load SequenceAnalysis/GenePrediction/maker/2.31.9;
module load Blast/ncbi-blast/2.9.0+;

WORKDIR=/data/users/mwannier/annotation/MAKER/flye
UNIPROT=/data/courses/assembly-annotation-course/CDS_annotation/uniprot-plant_reviewed.fasta

cd $WORKDIR

#maker

cd flye.maker.output

#gff3_merge -d flye_master_datastore_index.log
#fasta_merge -d flye_master_datastore_index.log

#maker_map_ids --prefix 'C24_' flye.all.gff > flye.all.id.map

#map_fasta_ids flye.all.id.map flye.all.maker.proteins.fasta
map_fasta_ids flye.all.id.map flye.all.maker.transcripts.fasta

#map_gff_ids flye.all.id.map flye.all.gff

#makeblastdb -in $UNIPROT -dbtype prot -out blast_db

#blastp -query flye.all.maker.proteins.fasta -db blast_db \
 #       -num_threads 10 -outfmt 6 -evalue 1e-10 > flye_blastp_output

#maker_functional_fasta $UNIPROT flye_blastp_output \
        #flye.all.maker.proteins.fasta > flye_functional.fasta

#maker_functional_gff $UNIPROT flye_blastp_output flye.all.gff > flye_functional.gff