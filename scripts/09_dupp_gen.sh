#!/usr/bin/env bash

#SBATCH --cpus-per-task=15
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=dupp_gen
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/data/users/mwannier/dupp_gen_%j.o
#SBATCH --error=/data/users/mwannier/dupp_gen_%j.e
#SBATCH --partition=pall

module load UHTS/Analysis/SeqKit/0.13.2
module add Emboss/EMBOSS/6.6.0
module load Phylogeny/paml/4.9j
module load Blast/ncbi-blast/2.10.1+;

NNU=/data/courses/assembly-annotation-course/CDS_annotation
OUTDIR=/data/users/mwannier/annotation/dupp_gen
ATH_GFF=/data/users/mwannier/annotation/MAKER/canu/canu.maker.output/canu.all.gff
CANU=/data/users/mwannier/annotation/MAKER/canu/canu.maker.output
WORKDIR=/data/users/mwannier/annotation/dupp_gen/canu
DUPGEN=/data/users/mwannier/annotation/DupGen_finder
SCRIPTS=/data/courses/assembly-annotation-course/CDS_annotation/scripts

cd $WORKDIR

# 1) Prepare input files

#gff file

#awk '{sub(/ID=/, "", $9); sub(/;.*/, "", $9); sub(/:.*/, "", $9); printf ("%s\t%s\t%s\t%s\n", $1, $9, $4, $5)}' $ATH_GFF | \
#grep '\-RA' \
#> Ath.gff

#cat Ath.gff $NNU/NNU_mRNA_single_model.gff > Ath_Nnu.gff

#blast file
# Create a reference database
#makeblastdb -in $CANU/canu.all.maker.proteins.fasta -dbtype prot -title Ath -out $WORKDIR/blast_db/Ath
# Align protein query sequences against the reference database
#blastp -query $CANU/canu.all.maker.proteins.fasta -db $WORKDIR/blast_db/Ath -evalue 1e-10 -max_target_seqs 5 -outfmt 6 -num_threads $SLURM_CPUS_PER_TASK -out $WORKDIR/Ath.blast

#blastp -query $NNU/NNU.pep.fa.ref.single_model -db $WORKDIR/blast_db/Ath -evalue 1e-10 -max_target_seqs 5 -outfmt 6 -num_threads $SLURM_CPUS_PER_TASK -out $WORKDIR/Ath_Nnu.blast

#2) Run DupGen_finder
#$DUPGEN/DupGen_finder.pl -i $WORKDIR -t Ath -c Nnu -o $OUTDIR/canu



#Estimate genome duplication events
#1. Create two fasta files with CDS and protein sequences of WGD genes.

#cut -f 1 $WORKDIR/Ath.wgd.genes > ID.txt # create a list of sequences to subset
#seqkit grep -f ID.txt $CANU/canu.all.maker.transcripts.fasta -o C24.wgd.genes.fa
#seqkit translate C24.wgd.genes.fa -o C24.wgd.genes.fa.proteins
#sed -i 's/_frame=1/_p/' C24.wgd.genes.fa.proteins # modify ID of protein sequences

#2. create a list of WGD gene pairs

#cut -f 1,3 $WORKDIR/Ath.wgd.pairs > wgd_pairs.csv
#sed -i 's/RA/RA_p/g' wgd_pairs.csv

#3. Copy wgd_pairs.csv, Ler.wgd.genes.fa, and Ler.wgd.genes.fa.proteins to a new directory

#cp wgd_pairs.csv $WORKDIR/estimation/
#cp C24.wgd.genes.fa $WORKDIR/estimation/
#cp C24.wgd.genes.fa.proteins $WORKDIR/estimation/

cd $WORKDIR/estimation

#4. Split fasta files

#$SCRIPTS/split_flat C24.wgd.genes.fa
#$SCRIPTS/split_flat C24.wgd.genes.fa.proteins

#5. Align proteins with bestflash_from_list. Output: protein alignments (.pair)

#$SCRIPTS/bestflash_from_list wgd_pairs.csv

#6. Generate codon-by-codon CDS alignments based on the protein alignments with pair_to_CDS_paml_pair. Output: CDS alignment in PAML format (.CDS_paml_pair)

#$SCRIPTS/pair_to_CDS_paml_pair C24 pair

#7. Use yn00 from PAML (PAML_yn00_from_CDS_pair script) on all WGD pair alignments to calculate Ka/Ks (omega) values

#$SCRIPTS/PAML_yn00_from_CDS_pair C24 > PAML_yn00_results 


#Calculate average Ks value for each syntenic block
#1. Parse PAML_yn00_results

#awk '{print($1,$1,$6,$7,$5)}' PAML_yn00_results|sed 's/ /\t/g' | \
 #sed 's/__x__/\t/'|sed 's/_p//g'|cut -f 1,2,4,5,6|sed 's/dN=//'| \
 #sed 's/dS=//'|sed 's/omega=//'|awk '$4<5' > C24.wgd.kaks

#2. Combine <Ler>.wgd.kaks and <Ler>.collinearity file using add_ka_ks_to_collinearity_file.pl (both files must be in the same directory)

perl $SCRIPTS/add_ka_ks_to_collinearity_file.pl C24

#3. Calculate average Ks value for each syntenic block using compute_ks_for_synteny_blocks.pl

perl $SCRIPTS/compute_ks_for_synteny_blocks.pl C24.collinearity.kaks

#4. Estimate Ks peaks from Ks distribution using plot_syntenic_blocks_ks_distri.py

python $SCRIPTS/plot_syntenic_blocks_ks_distri.py C24.synteny.blocks.ks.info 2 C24
