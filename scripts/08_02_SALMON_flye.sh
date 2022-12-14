#!/usr/bin/env bash

#SBATCH --cpus-per-task=15
#SBATCH --mem=4G
#SBATCH --time=05:00:00
#SBATCH --job-name=salmon
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/data/users/mwannier/salmon_flye_%j.o
#SBATCH --error=/data/users/mwannier/salmon_flye_%j.e
#SBATCH --partition=pall

#import module
module load UHTS/Analysis/busco/4.1.4;

#set variables
OUTDIR=/data/users/mwannier/annotation_course/Salmon/flye
BUSCO=/data/users/mwannier/annotation_course/BUSCO/flye
SALMONDIR=/data/users/mwannier/annotation_course/salmon-1.9.0_linux_x86_64/bin
MAKEOUT=/data/users/mwannier/annotation_course/MAKER/flye/flye.maker.output
SALMONINDEX=$OUTDIR/salmon_index
READS=/data/users/mwannier/FS22_Assembly/participant_2/RNAseq

#go to busco directory
cd $BUSCO

#run busco
busco -i $MAKEOUT/flye.all.maker.proteins.fasta -l brassicales_odb10 -m proteins -c 4 --out busco_out

#create salmon index
$SALMONDIR/salmon index -t $MAKEOUT/flye.all.maker.transcripts.fasta -i $SALMONINDEX -k 31

#run salmon quantification
$SALMONDIR/salmon quant -i $SALMONINDEX -l A -1 $READS/SRR1584462_1.fastq.gz -2 $READS/SRR1584462_2.fastq.gz \
                       --validateMappings -o $OUTDIR/flye_transcripts_quant
