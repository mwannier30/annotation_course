#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --job-name=dating
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --output=/data/users/mwannier/dating_%j.o
#SBATCH --error=/data/users/mwannier/dating_%j.e
#SBATCH --partition=pall

#import modules
module add Emboss/EMBOSS/6.6.0;

#set variable
WORKDIR=/data/users/mwannier/annotation/dating

#move to flye working directory
cd $WORKDIR/flye/splitted

sed -i 's/:/_/g' ../flye.mod.EDTA.intact.LTR.fa

#use split_flat to split the fasta file into multiple ones.
../../scripts/split_flat ../flye.mod.EDTA.intact.LTR.fa 

#use the perl script LTR to align the two LTR sequences of each TE
../../scripts/LTR C24_ N

#use the perl script date_pair to estimate the divergence time between the LTR sequences of each TE.
../../scripts/date_pair 

#move to canu working directory
cd $WORKDIR/canu/splitted

sed -i 's/:/_/g' ../canu.mod.EDTA.intact.LTR.fa

#use split_flat to split the fasta file into multiple ones.
../../scripts/split_flat ../canu.mod.EDTA.intact.LTR.fa 

#use the perl script LTR to align the two LTR sequences of each TE
../../scripts/LTR C24_ N

#use the perl script date_pair to estimate the divergence time between the LTR sequences of each TE.
../../scripts/date_pair 