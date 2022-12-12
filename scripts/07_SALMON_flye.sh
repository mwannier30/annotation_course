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

module load UHTS/Analysis/busco/4.1.4;

OUTDIR=/data/users/mwannier/annotation/Salmon/flye
BUSCO=/data/users/mwannier/annotation/BUSCO/flye
SALMONDIR=/data/users/mwannier/annotation/salmon-1.9.0_linux_x86_64/bin
MAKEOUT=/data/users/mwannier/annotation/MAKER/flye/flye.maker.output
SALMONINDEX=$OUTDIR/salmon_index
READS=/data/users/mwannier/FS22_Assembly/participant_2/Illumina


#cd $BUSCO

#busco -i $MAKEOUT/flye.all.maker.proteins.fasta -l brassicales_odb10 -m proteins -c 4 --out busco_out


#$SALMONDIR/salmon index -t $MAKEOUT/flye.all.maker.transcripts.fasta -i $SALMONINDEX -k 31


$SALMONDIR/salmon quant -i $SALMONINDEX -l A -1 $READS/ERR3624577_1.fastq.gz -2 $READS/ERR3624577_2.fastq.gz \
                       --validateMappings -o $OUTDIR/flye_transcripts_quant
