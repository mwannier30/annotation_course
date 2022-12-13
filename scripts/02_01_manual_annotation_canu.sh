#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --job-name=manual_annot_canu
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/mwannier/output_manuel_annot_canu_%j.o
#SBATCH --error=/data/users/mwannier/error_manuel_annot_canu_%j.e
#SBATCH --partition=pall

#load modules
module load UHTS/Analysis/SeqKit/0.13.2
module load Blast/ncbi-blast/2.9.0+
module load SequenceAnalysis/MultipleSequenceAlignment/clustalw2/2.1
module load Emboss/EMBOSS/6.6.0

#set variables
CANU=/data/users/mwannier/FS22_Assembly/polishing/canu/pilon_canon/canu.fasta
OUTDIR=/data/users/mwannier/annotation/manual/canu

#move to output directory
cd $OUTDIR

#Extract a long contig from the genome assembly and save it to a new fasta file
seqkit fx2tab -l -n $CANU | head -n 1 | cut -f 1 > ID.txt #print sequences length

seqkit grep -n -f ID.txt $CANU -o $OUTDIR/CANU_contig.fa #extract sequence from genome

#split the contig into 500bp windows
seqkit sliding -s 500 -W 500 -g $OUTDIR/CANU_contig.fa > CANU_contig_windows.fa

#create directory for blast databases
mkdir blast_db

#Blast 500bp-windows against themselves
makeblastdb -in CANU_contig_windows.fa -dbtype nucl -out blast_db/CANU_contig_windows
blastn \
    -query CANU_contig_windows.fa \
    -db blast_db/CANU_contig_windows \
    -num_threads 10 \
    -outfmt 6 \
    -perc_identity 80 \
    -max_hsps 1 > CANU_contig_windows.blastn

#find the 50 most abundant 500bp-windows
awk '{print $1}' $OUTDIR/CANU_contig_windows.blastn > $OUTDIR/blast_id.txt
sort $OUTDIR/blast_id.txt > $OUTDIR/blast_id_sorted.txt
uniq -c $OUTDIR/blast_id_sorted.txt > $OUTDIR/blast_id_cnt.txt
sort -n -r -k 1 $OUTDIR/blast_id_cnt.txt > $OUTDIR/blast_id_cnt_sorted.txt
awk '1<=NR && NR<=50{print $2}' $OUTDIR/blast_id_cnt_sorted.txt > $OUTDIR/blast_50max.txt

#save the 50 most abundant 500bp-windows into a new fasta file
seqkit \
    grep -n -f $OUTDIR/blast_50max.txt \
    $OUTDIR/CANU_contig_windows.fa \
    -o $OUTDIR/CANU_contig_windows_TOP50.fa


#Blast a window containing a candidate TE sequence against the genome
makeblastdb \
    -in $CANU \
    -dbtype nucl \
    -out $OUTDIR/blast_db/CANU


WINDOW='tig00000003_pilon_sliding:304001-304500'
awk -v seq=$WINDOW \
    -v RS='>' '$1 == seq {print RS $0}' $OUTDIR/CANU_contig_windows_TOP50.fa > $OUTDIR/CANU_window.fa


blastn \
    -query $OUTDIR/CANU_window.fa \
    -db $OUTDIR/blast_db/CANU \
    -outfmt 6 \
    -num_threads $SLURM_CPUS_PER_TASK \
    -perc_identity 80 \
    -max_hsps 1 \
    > $OUTDIR/CANU_window.blastn

#convert blast output to bed format
awk '$10 < $9 {print($2"\t"$10-1"\t"$9)}' $OUTDIR/CANU_window.blastn > $OUTDIR/CANU_window.bed
awk '$10 > $9 {print($2"\t"$9-1"\t"$10)}' $OUTDIR/CANU_window.blastn >> $OUTDIR/CANU_window.bed

#extract homologous TE sequences from the genome
seqkit subseq \
    --bed $OUTDIR/CANU_window.bed $CANU \
    -u 2000 \
    -d 2000 \
    > $OUTDIR/CANU_window.bed.fa

#Align TE sequences
clustalw2 -INFILE=$OUTDIR/CANU_window.bed.fa -OUTFILE=$OUTDIR/CANU_window.bed.aln

#write consensus sequence
cons \
    -sequence $OUTDIR/CANU_window.bed.aln \
    -outseq $OUTDIR/TE_family.cons \
    -name TE_cons

makeblastdb \
    -in /data/courses/assembly-annotation-course/Brassicaceae_repbase_all_march2019.fasta \
    -dbtype nucl \
    -out $WORKDIR/blast_db/Brassicaceae

blastn \
    -query $WORKDIR/TE5_family.cons \
    -db $WORKDIR/blast_db/Brassicaceae \
    -outfmt 6 \
    -num_threads $SLURM_CPUS_PER_TASK \
    -max_hsps 1 \
    > $WORKDIR/TE5_consensus_blasted.blastn