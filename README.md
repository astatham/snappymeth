snappymeth.py
---

Simple attempt to call allele-specific methylation (ASM) at CpG sitesfrom whole genome bisulfite data. Input data is

1. BAM aligned using [bwa-meth](https://github.com/brentp/bwa-meth)
2. VCF created from that BAM using [BisSNP](http://epigenome.usc.edu/publicationdata/bissnp2011/) (via [bwa-meth tabulate](https://github.com/brentp/bwa-meth))
3. FASTA of the reference sequence

Requirements
---

    pip install pysam pyvcf pyfasta fisher

Usage
---
	usage: snappymeth.py [-h] [--pair_distance PAIR_DISTANCE]
	                     [--max_depth MAX_DEPTH] [--min_per_allele MIN_PER_ALLELE]
	                     [--min_sites_in_region MIN_SITES_IN_REGION]
	                     [--min_mapping_quality MIN_MAPPING_QUALITY]
	                     [--min_base_quality MIN_BASE_QUALITY]
	                     input_vcf input_bam reference prefix

	snappymeth.py - Parses a paired bam and vcf and looks for regions of allele-
	specific methylation.

	positional arguments:
	  input_vcf             Input VCF file
	  input_bam             Input BAM file
	  reference             Reference FASTA file
	  prefix                Prefix for the sites and regions output csvs

	optional arguments:
	  -h, --help            show this help message and exit
	  --pair_distance PAIR_DISTANCE
	                        The distance in basepairs to search up and downstream
	                        from each SNP (default is 500)
	  --max_depth MAX_DEPTH
	                        Maximum number of reads to process at each SNP
	                        position (default is 8000)
	  --min_per_allele MIN_PER_ALLELE
	                        Minimum number of reads containing each allele to
	                        process a SNP position
	  --min_sites_in_region MIN_SITES_IN_REGION
	                        Minimum number of CpG sites linked to a SNP to perform
	                        a regional analysis
	  --min_mapping_quality MIN_MAPPING_QUALITY
	                        Minimum mapping quality score for a read to be
	                        considered
	  --min_base_quality MIN_BASE_QUALITY
	                        Minimum basecall quality score at the SNP for a read
	                        to be considered

Output
---

Two comma separated value files:

1. *prefix*.sites.csv - One line per CpG site linked to a heterozygous SNP, containing the following fields:
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

Example of ASM found using snappymeth.py
---

The output of snappymeth.py was sorted by region p-value, and then reads around the most significant SNP separated into separate bams by [splitSNP.py](https://github.com/astatham/splitSNP) then visualised with [IGV](https://www.broadinstitute.org/igv/).

![Example of ASM](http://i.imgur.com/8loeIYk.png)
