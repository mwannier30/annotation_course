#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --job-name=dating
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --output=/data/users/mwannier/dating_%j.o
#SBATCH --error=/data/users/mwannier/dating_%j.e
#SBATCH --partition=pall

module add Emboss/EMBOSS/6.6.0;

WORKDIR=/data/users/mwannier/annotation/dating

cd $WORKDIR/flye/splitted

sed -i 's/:/_/g' ../flye.mod.EDTA.intact.LTR.fa

../../scripts/split_flat ../flye.mod.EDTA.intact.LTR.fa 

../../scripts/LTR C24_ N

../../scripts/date_pair 


cd $WORKDIR/canu/splitted

sed -i 's/:/_/g' ../canu.mod.EDTA.intact.LTR.fa

../../scripts/split_flat ../canu.mod.EDTA.intact.LTR.fa 

../../scripts/LTR C24_ N

../../scripts/date_pair 