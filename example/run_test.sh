#!/bin/bash -e

check_no_lines_in_file() {
	if [ `wc -l < $1` -ne $2 ]
	then
		exit 1
	fi
}

# Clean up
rm -rf chr15.fa chr15.fa.gz output

# Download hg19 chr10 fasta
mkdir output
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr15.fa.gz
gunzip chr15.fa.gz

# SNP analysis using VCF
../snappymeth.py --region_bams PrEC.SNRPN.vcf.gz PrEC.SNRPN.bam chr15.fa output/PrEC.SNRPN.VCF
## check 13 sites 2 regions (plus 1 line for header)
check_no_lines_in_file output/PrEC.SNRPN.VCF.sites.csv 14
check_no_lines_in_file output/PrEC.SNRPN.VCF.regions.csv 3

# CpG site analysis
../snappymeth.py --input_type CpGs --min_per_allele 10 --region_bams PrEC.SNRPN.CpGs.csv.gz PrEC.SNRPN.bam chr15.fa output/PrEC.SNRPN.CpGs
## check 22 sites 3 regions (plus 1 line for header)
check_no_lines_in_file output/PrEC.SNRPN.CpGs.sites.csv 23
check_no_lines_in_file output/PrEC.SNRPN.CpGs.regions.csv 4

echo "Tests passed!"
