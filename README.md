snappymeth.py
---

[![Build Status](https://travis-ci.org/astatham/snappymeth.svg?branch=master)](https://travis-ci.org/astatham/snappymeth)

*snappymeth.py* was written to discover allele-specific methylation (ASM) of CpG sites and small regions from whole genome bisulfite sequencing (WGBS) data.

Two approaches have been implemented:

1. Using heterozygous SNPs (supplied in a VCF) to separate read pairs correspoing to each allele.
2. Using intermediately methylated CpG sites as pseudo-heterozygous SNPs to separate read pairs.

Counts of methylated and unmethylated bases sequenced in each alleles reads are summed at each CpG site surrounding the SNP and a fisher exact test p-value (two-tailed) calculated. If enough CpG sites are present then a regional analysis is performed by summing the counts per-allele at all covered CpG sites. If desired, reads for each allele at regions that meet a p-value cutoff are exported as separate BAM files, and IGV screenshots of these reads can be automatically exported.

Input data is:

1. BAM aligned using [bwa-meth](https://github.com/brentp/bwa-meth)
2. Sites of potential allele-specific methylation; either
	* a VCF created from that BAM using [BisSNP](http://epigenome.usc.edu/publicationdata/bissnp2011/) (via [bwa-meth tabulate](https://github.com/brentp/bwa-meth))
	* a CSV file containing 4 columns summarising the methylation state of all CpG sites generated from that BAM
		1. chr - chromosome 
		2. position - 0-based position of the CpG site
		3. M - count of methylated reads at this CpG site
		4. U - count of unmethylated reads at this CpG site
3. FASTA of the reference sequence

Requirements
---

    pip install pysam pyvcf pyfasta fisher pandas

Usage
---
	usage: snappymeth.py [-h] [--input_type {VCF,CpGs}] [--VCF_sample VCF_SAMPLE]
	                     [--pair_distance PAIR_DISTANCE] [--max_depth MAX_DEPTH]
	                     [--min_per_allele MIN_PER_ALLELE]
	                     [--min_sites_in_region MIN_SITES_IN_REGION]
	                     [--min_mapping_quality MIN_MAPPING_QUALITY]
	                     [--min_base_quality MIN_BASE_QUALITY] [--region_bams]
	                     [--fisher_cutoff FISHER_CUTOFF] [--IGV_screenshot]
	                     input_file input_bam reference prefix

	snappymeth.py - Discover sites and regions of allele specific methylation from
	whole genome bisulfite sequencing data by counting CpG methylation on alleles
	separately. Reads can be separated by either a heterozygous SNP (when a VCF is
	supplied), or by the methylation status of a single CpG site. Both analyses
	modes require sufficient sequencing coverage of both alleles (default is 10x).

	positional arguments:
	  input_file            Input VCF/CpG sites file, gzip compressed.
	  input_bam             Input BAM file
	  reference             Reference FASTA file
	  prefix                Prefix for the sites and regions output csvs

	optional arguments:
	  -h, --help            show this help message and exit
	  --input_type {VCF,CpGs}
	                        Whether the input_file is a VCF (default) or a csv of
	                        methylation counts at CpG sites with the format
	                        'chr,position,M,U' where the fields are chromosome
	                        name, 0-based position of the CpG site, count of
	                        methylated bases sequenced at this site and count of
	                        unmethylated bases sequenced.
	  --VCF_sample VCF_SAMPLE
	                        The sample in the VCF to be processed - either as the
	                        sample name or numeric index (0-based). Default is 0,
	                        the first sample.
	  --pair_distance PAIR_DISTANCE
	                        The distance in basepairs to search up and downstream
	                        from each position (default is 500).
	  --max_depth MAX_DEPTH
	                        Maximum number of reads to process at each position
	                        (default is 8000).
	  --min_per_allele MIN_PER_ALLELE
	                        Minimum number of reads containing each allele to
	                        process a position.
	  --min_sites_in_region MIN_SITES_IN_REGION
	                        Minimum number of CpG sites linked to a SNP to perform
	                        a regional analysis.
	  --min_mapping_quality MIN_MAPPING_QUALITY
	                        Minimum mapping quality score for a read to be
	                        considered.
	  --min_base_quality MIN_BASE_QUALITY
	                        Minimum basecall quality score at the SNP for a read
	                        to be considered.
	  --region_bams         Specity to output BAM files per allele when the
	                        regional fisher exact p-value is less than the cutoff
	                        specified by --fisher_cutoff.
	  --fisher_cutoff FISHER_CUTOFF
	                        Regional fisher exact p-value cutoff for a regional
	                        BAM to be created/IGV screenshot be taken (default is
	                        0.0001).
	  --IGV_screenshot      Specity to take IGV screenshots of each region that
	                        passes --fisher_cutoff. Requires that IGV be running
	                        on the local machine and listening on port 60151

Output
---

Two comma separated value files:

1. *prefix*.sites.csv - One line per CpG site linked to a heterozygous SNP (or intermediately methylated CpG site), containing the following fields:
	* SNP.chr - SNP chromosome
	* SNP.pos - SNP position
	* SNP.ref - reference allele
	* SNP.alt - alternate allele
	* CpG.pos - position of the CpG site linked to the SNP (via a read-pair)
	* ref.A/C/G/T/N - count of A/C/G/T/N bases sequenced at the CpG methylation site in reference allele reads
	* alt.A/C/G/T/N - count of A/C/G/T/N bases sequenced at the CpG methylation site in alternate allele reads
	* ref.cov/alt.cov - sequencing coverage (only counting Cs and Ts) for reference and alternate alleles
	* ref.meth/alt.meth - methylation ratio for reference and alternate alleles
	* meth.diff - difference between ref.meth and alt.meth
	* p.val - two-tailed fisher exact test p-value of count of Cs and Ts being different between ref and alt alleles
2. *prefix*.regions.csv - One line per region around each heterozygous SNP containing *min\_sites\_in\_region* CpG sites. This file contains similar fields to *prefix*.sites.cov - with counts being the sum of all CpG sites contained in the region. Additional columns are:
	* first.CpG - position of the first CpG site found in this region
	* last.CpG - position of the last CpG site found in this region
	* nCG - number of CpG sites found in this region

If *--region_bams* is specified, then for each region which meets the *--fisher_cutoff* two (indexed) BAM files are created containing the reads for each allele. Additionally if *--IGV_screenshot* is specified, then these BAM files are loaded into a local running instance of IGV and a screenshot is saved.

Example usage of snappymeth.py
---

A subset of a WGBS analysis of prostate epithelial cells (PrEC) for the [SNRPN imprinted locus](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2892089/) (hg19 chr15:25197573-25203941) is supplied in the wxample/ directory. The following command will use the genotyping data from BisSNP find 2 heterozygous SNPs in this region that separate the methylated and unmethylated alleles.

    ../snappymeth.py --region_bams PrEC.SNRPN.vcf.gz PrEC.SNRPN.bam hg19.fa PrEC.SNRPN.VCF

![Example of ASM at the SNRPN imprinted region](http://i.imgur.com/Pg6CP5H.png)
Visualised using [IGV](https://www.broadinstitute.org/igv/)

TODO
---

* Perform IGV screenshots in a separate process
* Link adjacent ASM regions from CpG sites analysis into methyl-haplotypes

Acknowledgements
---
* IGV.py included in the repository was copied from [Brent Pedersen's](https://github.com/brentp) [bio-playground github repository](https://github.com/brentp/bio-playground/blob/master/igv/igv.py). Cheers!

