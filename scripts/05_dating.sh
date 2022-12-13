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
module load UHTS/Analysis/BEDTools/2.29.2;

#set variables
CANU=/data/users/mwannier/annotation/EDTA/canu/canu.fasta.mod.EDTA.anno/
CANU_ref=/data/users/mwannier/FS22_Assembly/polishing/canu/pilon_canon/canu.fasta
FLYE=/data/users/mwannier/annotation/EDTA/flye/flye.fasta.mod.EDTA.anno/
FLYE_ref=/data/users/mwannier/FS22_Assembly/polishing/flye/pilon_flye/flye.fasta
BRASSICACEA=/data/users/mwannier/annotation/EDTA/canu/TEsorter/brassicacea/
WORKDIR=/data/users/mwannier/annotation/dating

#move to working directory
cd $WORKDIR

#create a GFF with only LTR retrotransposons (this time from $genome.mod.EDTA.intact.gff3):
awk '$3~/retrotransposon/' $CANU/canu.fasta.mod.EDTA.intact.gff3 > canu/canu.mod.EDTA.intact.gff3_LTR
awk '$3~/retrotransposon/' $FLYE/flye.fasta.mod.EDTA.homo.gff3 > flye/flye.mod.EDTA.intact.gff3_LTR

#reformat gff (shorten IDs and rearrange column order):
sed -i 's/ID.\+Name=//' canu/canu.mod.EDTA.intact.gff3_LTR
sed -i 's/;.\+//' canu/canu.mod.EDTA.intact.gff3_LTR
sed -i 's/_pi\s/_pilon\t/g' canu/canu.mod.EDTA.intact.gff3_LTR
sed -i 's/_pil\s/_pilon\t/g' canu/canu.mod.EDTA.intact.gff3_LTR
sed -i 's/_pilo\s/_pilon\t/g' canu/canu.mod.EDTA.intact.gff3_LTR
awk '{print($1,$2,$9,$4,$5,$6,$7,$8,$3)}' \
canu/canu.mod.EDTA.intact.gff3_LTR |sed 's/ /\t/g' > canu/canu.mod.EDTA.intact.gff3_LTR_ref

sed -i 's/ID.\+Name=//' flye/flye.mod.EDTA.intact.gff3_LTR
sed -i 's/;.\+//' flye/flye.mod.EDTA.intact.gff3_LTR
sed -i 's/_pi\s/_pilon\t/g' flye/flye.mod.EDTA.intact.gff3_LTR
sed -i 's/_pil\s/_pilon\t/g' flye/flye.mod.EDTA.intact.gff3_LTR
sed -i 's/_pilo\s/_pilon\t/g' flye/flye.mod.EDTA.intact.gff3_LTR
awk '{print($1,$2,$9,$4,$5,$6,$7,$8,$3)}' \
flye/flye.mod.EDTA.intact.gff3_LTR |sed 's/ /\t/g' > flye/flye.mod.EDTA.intact.gff3_LTR_ref


#extract sequences with bedtools getfasta:
bedtools getfasta -fi $CANU_ref -bed \
canu/canu.mod.EDTA.intact.gff3_LTR_ref -name \
> canu/canu.mod.EDTA.intact.LTR.fa

bedtools getfasta -fi $FLYE_ref -bed \
flye/flye.mod.EDTA.intact.gff3_LTR_ref -name \
> flye/flye.mod.EDTA.intact.LTR.fa

#Rename LTR retrotransposons

#replace ":" with "_"
sed -i 's/:/_/g' canu/canu.mod.EDTA.intact.LTR.fa

sed -i 's/:/_/g' flye/flye.mod.EDTA.intact.LTR.fa

#All LTR retrotransposons need a "common pattern" inside their fasta header. For
#example you can put the abbreviation of your Arabidopsis accession (in this case Ler) in
#front of the TE names:
sed -i 's/>/>C24_/' canu/canu.mod.EDTA.intact.LTR.fa
sed -i 's/>/>C24_/' flye/flye.mod.EDTA.intact.LTR.fa


