#!/bin/bash -e

# Download hg19 chr10 fasta
rm -rf chr15.fa chr15.fa.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr15.fa.gz
gunzip chr15.fa.gz

# SNP analysis using VCF
../snappymeth.py --region_bams PrEC.SNRPN.vcf.gz PrEC.SNRPN.bam chr15.fa PrEC.SNRPN.VCF

# CpG site analysis
../snappymeth.py --input_type CpGs --min_per_allele 10 --region_bams PrEC.SNRPN.CpGs.csv.gz PrEC.SNRPN.bam chr15.fa PrEC.SNRPN.VCF
