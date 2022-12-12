#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --job-name=phylo_analysis
#SBATCH --mail-user=maelle.wannier@students.unibe.ch
#SBATCH --output=/data/users/mwannier/phylo_analysis_%j.o
#SBATCH --error=/data/users/mwannier/phylo_analysis_%j.e
#SBATCH --partition=pall

module load UHTS/Analysis/SeqKit/0.13.2;
module load SequenceAnalysis/MultipleSequenceAlignment/clustal-omega/1.2.4;
module load Phylogeny/FastTree/2.1.10

OUTDIR=/data/users/mwannier/annotation/phylogeny

CANU_TESORTER=/data/users/mwannier/annotation/EDTA/canu/TEsorter/getfasta_output.rexdb-plant.cls.pep
FLYE_TESORTER=/data/users/mwannier/annotation/EDTA/flye/TEsorter/getfasta_output.rexdb-plant.cls.pep
BRASSICACEA_TESORTER=/data/users/mwannier/annotation/EDTA/canu/TEsorter/brassicacea/Brassicaceae_repbase_all_march2019.fasta.rexdb-plant.cls.pep

CANU_TSV=/data/users/mwannier/annotation/EDTA/canu/TEsorter/getfasta_output.rexdb-plant.cls.tsv
FLYE_TSV=/data/users/mwannier/annotation/EDTA/canu/TEsorter/brassicacea/Brassicaceae_repbase_all_march2019.fasta.rexdb-plant.cls.tsv
BRA_TSV=/data/users/mwannier/annotation/EDTA/flye/TEsorter/getfasta_output.rexdb-plant.cls.tsv

cd $OUTDIR

#Create IDs files

grep Ty3-RT $CANU_TESORTER | sed 's/>//' | sed 's/ .\+//' > ID_gypsy_canu.txt
grep Ty1-RT $CANU_TESORTER | sed 's/>//' | sed 's/ .\+//' > ID_copia_canu.txt

grep Ty3-RT $FLYE_TESORTER | sed 's/>//' | sed 's/ .\+//' > ID_gypsy_flye.txt
grep Ty1-RT $FLYE_TESORTER | sed 's/>//' | sed 's/ .\+//' > ID_copia_flye.txt

grep Ty3-RT $BRASSICACEA_TESORTER | sed 's/>//' | sed 's/ .\+//' > ID_gypsy_bra.txt
grep Ty1-RT $BRASSICACEA_TESORTER | sed 's/>//' | sed 's/ .\+//' > ID_copia_bra.txt

#Extract RT protein sequences from $.rexdb-plant.cls.pep using seqkit grep
seqkit grep -f ID_copia_canu.txt $CANU_TESORTER -o canu/result_canu_copia
seqkit grep -f ID_gypsy_canu.txt $CANU_TESORTER -o canu/result_canu_gypsy

seqkit grep -f ID_copia_flye.txt $FLYE_TESORTER -o flye/result_flye_copia
seqkit grep -f ID_gypsy_flye.txt $FLYE_TESORTER -o flye/result_flye_gypsy

seqkit grep -f ID_copia_bra.txt $BRASSICACEA_TESORTER -o brassicacea/result_brassicacea_copia
seqkit grep -f ID_gypsy_bra.txt $BRASSICACEA_TESORTER -o brassicacea/result_brassicacea_gypsy

#Concatenate RTs from both Brassicaceae and Arabidopsis TEs into one fasta file
cat canu/result_canu_copia brassicacea/result_brassicacea_copia > fasta_files/canu_copia.fa
cat canu/result_canu_gypsy brassicacea/result_brassicacea_gypsy > fasta_files/canu_gypsy.fa

cat flye/result_flye_copia brassicacea/result_brassicacea_copia > fasta_files/flye_copia.fa
cat flye/result_flye_gypsy brassicacea/result_brassicacea_gypsy > fasta_files/flye_gypsy.fa

#shorten identifiers of RT sequences and replace ":" with "_"
sed -i 's/#.\+//' fasta_files/canu_copia.fa 
sed -i 's/:/_/' fasta_files/canu_copia.fa

sed -i 's/#.\+//' fasta_files/canu_gypsy.fa
sed -i 's/:/_/' fasta_files/canu_gypsy.fa

sed -i 's/#.\+//' fasta_files/flye_copia.fa
sed -i 's/:/_/' fasta_files/flye_copia.fa

sed -i 's/#.\+//' fasta_files/flye_gypsy.fa
sed -i 's/:/_/' fasta_files/flye_gypsy.fa

#Align the sequences with clustal omega
clustalo -i fasta_files/canu_copia.fa  -o alignment/canu_copia_aligned.fa
clustalo -i fasta_files/canu_gypsy.fa  -o alignment/canu_gypsy_aligned.fa
clustalo -i fasta_files/flye_copia.fa  -o alignment/flye_copia_aligned.fa
clustalo -i fasta_files/flye_gypsy.fa -o alignment/flye_gypsy_aligned.fa

#Infer approximately-maximum-likelihood phylogenetic tree with FastTree
FastTree -out tree/canu_copia_tree alignment/canu_copia_aligned.fa
FastTree -out tree/canu_gypsy_tree alignment/canu_gypsy_aligned.fa
FastTree -out tree/flye_copia_tree alignment/flye_copia_aligned.fa
FastTree -out tree/flye_gypsy_tree alignment/flye_gypsy_aligned.fa

#create a list with "TE_ID" and "color_ID" separated by a space:
grep -e "Retand RT" $CANU_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' > annot/ID_canu.txt
grep -e "Retand RT" $FLYE_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' > annot/ID_flye.txt
grep -e "Retand RT" $BRA_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' > annot/ID_bra.txt

grep -e "Tekay RT" $CANU_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_canu.txt
grep -e "Tekay RT" $FLYE_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_flye.txt
grep -e "Tekay RT" $BRA_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_bra.txt

grep -e "Athila RT" $CANU_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_canu.txt
grep -e "Athila RT" $FLYE_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_flye.txt
grep -e "Athila RT" $BRA_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_bra.txt

grep -e "Reina RT" $CANU_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_canu.txt
grep -e "Reina RT" $FLYE_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_flye.txt
grep -e "Reina RT" $BRA_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_bra.txt

grep -e "Galadriel RT" $CANU_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_canu.txt
grep -e "Galadriel RT" $FLYE_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_flye.txt
grep -e "Galadriel RT" $BRA_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_bra.txt

grep -e "TatII RT" $CANU_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_canu.txt
grep -e "TatII RT" $FLYE_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_flye.txt
grep -e "TatII RT" $BRA_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_bra.txt

grep -e "CRM RT" $CANU_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_canu.txt
grep -e "CRM RT" $FLYE_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_flye.txt
grep -e "CRM RT" $BRA_TSV | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' >> annot/ID_bra.txt
