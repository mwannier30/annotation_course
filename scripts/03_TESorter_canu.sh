#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --job-name=TESorter
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --output=/data/users/mwannier/output_tesorter_%j.o
#SBATCH --error=/data/users/mwannier/error_tesorter_%j.e
#SBATCH --partition=pall

module load UHTS/Analysis/BEDTools/2.29.2;

PROJECTDIR=/data/users/mwannier/annotation

OUTDIR=$PROJECTDIR/EDTA/
WORKDIR=$OUTDIR/canu
CANU=/data/users/mwannier/FS22_Assembly/polishing/canu/pilon_canon/canu.fasta
BRASSICACEA=/data/courses/assembly-annotation-course/CDS_annotation/Brassicaceae_repbase_all_march2019.fasta

cd $WORKDIR

awk '$3~/retrotransposon/' canu.fasta.mod.EDTA.TEanno.gff3 > canu.fasta.mod.EDTA.TEanno.gff3_edited

grep -v LTR canu.fasta.mod.EDTA.TEanno.gff3 >> canu.fasta.mod.EDTA.TEanno.gff3_edited


sed 's/\(ID=.*;Name=\)\(\w*\)\(;.*\)/\2/g' canu.fasta.mod.EDTA.TEanno.gff3_edited \
| awk '{printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $9, $4, $5, $6, $7, $8, $3)}' \
> parsed.tab

sed -i 's/_pi\s/_pilon\t/g' $WORKDIR/parsed.tab
sed -i 's/_pil\s/_pilon\t/g' $WORKDIR/parsed.tab
sed -i 's/_pilo\s/_pilon\t/g' $WORKDIR/parsed.tab

bedtools getfasta -s \
-name \
-fi $CANU \
-bed parsed.tab \
-fo getfasta_output

OUT=/data/users/mwannier/annotation/EDTA/canu/getfasta_output


cd TEsorter

singularity exec \
--bind $PROJECTDIR \
/data/courses/assembly-annotation-course/containers2/TEsorter_1.3.0.sif \
TEsorter $OUT -db rexdb-plant

cd brassicacea 

singularity exec \
--bind $PROJECTDIR \
--bind $BRASSICACEA \
/data/courses/assembly-annotation-course/containers2/TEsorter_1.3.0.sif \
TEsorter $BRASSICACEA -db rexdb-plant