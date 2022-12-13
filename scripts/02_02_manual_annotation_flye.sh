#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --job-name=manual_annot_flye
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/mwannier/output_manuel_annot_flye_%j.o
#SBATCH --error=/data/users/mwannier/error_manuel_annot_flye_%j.e
#SBATCH --partition=pshort

#load modules
module load UHTS/Analysis/SeqKit/0.13.2
module load Blast/ncbi-blast/2.9.0+
module load SequenceAnalysis/MultipleSequenceAlignment/clustalw2/2.1
module load Emboss/EMBOSS/6.6.0

#set variables
FLYE=/data/users/mwannier/FS22_Assembly/polishing/flye/pilon_flye/flye.fasta
OUTDIR=/data/users/mwannier/annotation/manual/flye

#move to output directory
cd $OUTDIR

#Extract a long contig from the genome assembly and save it to a new fasta file
seqkit fx2tab -l -n $FLYE | head -n 1 | cut -f 1 > ID.txt #print sequences length

seqkit grep -n -f ID.txt $FLYE -o $OUTDIR/FLYE_contig.fa #extract sequence from genome

#split the contig into 500bp windows
seqkit sliding -s 500 -W 500 -g $OUTDIR/FLYE_contig.fa > FLYE_contig_windows.fa

#create directory for blast databases
mkdir blast_db

#Blast 500bp-windows against themselves
makeblastdb -in FLYE_contig_windows.fa -dbtype nucl -out blast_db/FLYE_contig_windows
blastn \
    -query FLYE_contig_windows.fa \
    -db blast_db/FLYE_contig_windows \
    -num_threads 10 \
    -outfmt 6 \
    -perc_identity 80 \
    -max_hsps 1 > FLYE_contig_windows.blastn

#find the 50 most abundant 500bp-windows
awk '{print $1}' $OUTDIR/FLYE_contig_windows.blastn > $OUTDIR/blast_id.txt
sort $OUTDIR/blast_id.txt > $OUTDIR/blast_id_sorted.txt
uniq -c $OUTDIR/blast_id_sorted.txt > $OUTDIR/blast_id_cnt.txt
sort -n -r -k 1 $OUTDIR/blast_id_cnt.txt > $OUTDIR/blast_id_cnt_sorted.txt
awk '1<=NR && NR<=50{print $2}' $OUTDIR/blast_id_cnt_sorted.txt > $OUTDIR/blast_50max.txt

#save the 50 most abundant 500bp-windows into a new fasta file
seqkit \
    grep -n -f $OUTDIR/blast_50max.txt \
    $OUTDIR/FLYE_contig_windows.fa \
    -o $OUTDIR/flye_contig_windows_TOP50.fa


#Blast a window containing a candidate TE sequence against the genome
makeblastdb \
    -in $FLYE \
    -dbtype nucl \
    -out $OUTDIR/blast_db/flye


WINDOW='contig_10_pilon_sliding:1034001-1034500'
awk -v seq=$WINDOW \
    -v RS='>' '$1 == seq {print RS $0}' $OUTDIR/flye_contig_windows_TOP50.fa > $OUTDIR/flye_window.fa


blastn \
    -query $OUTDIR/flye_window.fa \
    -db $OUTDIR/blast_db/flye \
    -outfmt 6 \
    -num_threads $SLURM_CPUS_PER_TASK \
    -perc_identity 80 \
    -qcov_hsp_perc 80 \
    -max_hsps 1 \
    > $OUTDIR/flye_window.blastn

#convert blast output to bed format
awk '$10 < $9 {print($2"\t"$10-1"\t"$9)}' $OUTDIR/flye_window.blastn > $OUTDIR/flye_window.bed
awk '$10 > $9 {print($2"\t"$9-1"\t"$10)}' $OUTDIR/flye_window.blastn >> $OUTDIR/flye_window.bed

#extract homologous TE sequences from the genome
seqkit subseq \
    --bed $OUTDIR/flye_window.bed $FLYE \
    -u 2000 \
    -d 2000 \
    > $OUTDIR/flye_window.bed.fa

#Align TE sequences
clustalw2 -INFILE=$OUTDIR/flye_window.bed.fa -OUTFILE=$OUTDIR/flye_window.bed.aln

#write consensus sequence
cons \
    -sequence $OUTDIR/flye_window.bed.aln \
    -outseq $OUTDIR/TE_family.cons \
    -name TE_cons


makeblastdb \
    -in /data/courses/assembly-annotation-course/Brassicaceae_repbase_all_march2019.fasta \
    -dbtype nucl \
    -out $OUTDIR/blast_db/Brassicaceae

blastn \
    -query $OUTDIR/TE_family.cons \
    -db $OUTDIR/blast_db/Brassicaceae \
    -outfmt 6 \
    -num_threads $SLURM_CPUS_PER_TASK \
    -max_hsps 1 \
    > $OUTDIR/TE_consensus_blasted.blastn