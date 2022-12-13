#!/usr/bin/env bash

#SBATCH --cpus-per-task=15
#SBATCH --mem=4G
#SBATCH --time=05:00:00
#SBATCH --job-name=salmon
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/data/users/mwannier/salmon_%j.o
#SBATCH --error=/data/users/mwannier/salmon_%j.e
#SBATCH --partition=pall

#import modules
module load UHTS/Analysis/busco/4.1.4;

#set variables
OUTDIR=/data/users/mwannier/annotation/Salmon/canu
BUSCO=/data/users/mwannier/annotation/BUSCO/canu
SALMONDIR=/data/users/mwannier/annotation/salmon-1.9.0_linux_x86_64/bin
MAKEOUT=/data/users/mwannier/annotation/MAKER/canu/canu.maker.output
SALMONINDEX=$OUTDIR/salmon_index
READS=/data/users/mwannier/FS22_Assembly/participant_2/Illumina

#go to busco directory
cd $BUSCO

#run busco
busco -i $MAKEOUT/canu.all.maker.proteins.fasta -l brassicales_odb10 -m proteins -c 4 --out busco_out

#create salmon index
$SALMONDIR/salmon index -t $MAKEOUT/canu.all.maker.transcripts.fasta -i $SALMONINDEX -k 31

#run salmon quantification
$SALMONDIR/salmon quant -i $SALMONINDEX -l A -1 $READS/ERR3624577_1.fastq.gz -2 $READS/ERR3624577_2.fastq.gz \
                       --validateMappings -o $OUTDIR/canu_transcripts_quant