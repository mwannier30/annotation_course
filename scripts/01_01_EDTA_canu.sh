#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=EDTA
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/mwannier/output_EDTA_%j.o
#SBATCH --error=/data/users/mwannier/error_EDTA_%j.e
#SBATCH --partition=pcourseassembly

#set variables
CONTAINER=/data/courses/assembly-annotation-course/containers2
WORKDIR=/data/users/mwannier/annotation/EDTA
ASSEMBLYDIR=/data/users/mwannier/FS22_Assembly
GENOME=/data/users/mwannier/FS22_Assembly/polishing/canu/pilon_canon/canu.fasta
OUTDIR=$WORKDIR/canu

#move to output directory
cd $OUTDIR

#run EDTA on canu genome
singularity exec \
--bind $CONTAINER \
--bind $WORKDIR \
--bind $ASSEMBLYDIR \
$CONTAINER/EDTA_v1.9.6.sif \
EDTA.pl \
--genome $GENOME \
--species others \
--step all \
--cds $WORKDIR/TAIR10_cds_20110103_representative_gene_model \
--anno 1 \
--threads 16
