#!/usr/bin/bash

# Set base URLs
FASTA_URL="https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.chr.gtf.gz"

mkdir -p references

# Download the FASTA chr1 file
curl -L "${FASTA_URL}" -o references/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz

# Download the GTF file
curl -L "${GTF_URL}" -o references/Homo_sapiens.GRCh38.116.chr.gtf.gz

# Unzip the files

gunzip references/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
gunzip references/Homo_sapiens.GRCh38.116.chr.gtf.gz

# separate only chr1 annotations from the GTF file
awk '$1 == "1"' references/Homo_sapiens.GRCh38.116.chr.gtf > references/Homo_sapiens.GRCh38.116.chr1.gtf

# remove complete gtf file
rm references/Homo_sapiens.GRCh38.116.chr.gtf

echo -e "\nHuman reference chromosome 1 - genome and annotation files downloaded.\n"

# Download the fastq files

mkdir -p samples

FASTQ_URL="https://www.encodeproject.org/files/"

curl -L "${FASTQ_URL}"ENCFF827DUW/@@download/ENCFF827DUW.fastq.gz -o samples/ENCFF827DUW.fastq.gz
curl -L "${FASTQ_URL}"ENCFF708BOP/@@download/ENCFF708BOP.fastq.gz -o samples/ENCFF708BOP.fastq.gz
curl -L "${FASTQ_URL}"ENCFF785KVJ/@@download/ENCFF785KVJ.fastq.gz -o samples/ENCFF785KVJ.fastq.gz
curl -L "${FASTQ_URL}"ENCFF260AWP/@@download/ENCFF260AWP.fastq.gz -o samples/ENCFF260AWP.fastq.gz

echo -e "\nENCODE Project RUSH AD PacBio - Downloaded ENCFF827DUW, ENCFF708BOP, ENCFF785KVJ and ENCFF260AWP fastq files for dorsolateral pre-frontal cortex samples.\n"
