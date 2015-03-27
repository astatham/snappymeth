snappymeth.py
---

Simple attempt to call allele-specific methylation (ASM) at CpG sitesfrom whole genome bisulfite data. Input data is

1. BAM aligned using [bwa-meth](https://github.com/brentp/bwa-meth)
2. VCF created from that BAM using [BisSNP](http://epigenome.usc.edu/publicationdata/bissnp2011/) (via [bwa-meth tabulate](https://github.com/brentp/bwa-meth))
3. FASTA of the reference sequence

Usage
---
	usage: snappymeth.py [-h] [--pair_distance PAIR_DISTANCE]
	                     [--max_depth MAX_DEPTH]
	                     input_vcf input_bam reference output

	snappymeth.py - Parses a paired bam and vcf and looks for regions of allele-
	specific methylation.

	positional arguments:
	  input_vcf             Input VCF file
	  input_bam             Input BAM file
	  reference             Reference FASTA file
	  output                Filename for the output table (csv format)

	optional arguments:
	  -h, --help            show this help message and exit
	  --pair_distance PAIR_DISTANCE
	                        The distance in basepairs to search up and downstream
	                        from each SNP (default is 500)
	  --max_depth MAX_DEPTH
	                        Maximum number of reads to process at the specified
	                        SNP position (default is 1000000)


Output
---

A comma separated value file with the following columns (one per CpG site):

* SNP.chr - SNP chromosome
* SNP.pos - SNP position
* SNP.ref - reference allele
* SNP.alt - alternate allele
* CpG.pos - position of the CpG site linked to the SNP (via a read-pair)
* ref.A - count of A bases sequenced at the CpG methylation site in reference allele reads
* ref.C
* ref.G
* ref.T
* ref.N
* alt.A - count of A bases sequenced at the CpG methylation site in alternate allele reads
* alt.C
* alt.G
* alt.T
* alt.N

Example of ASM found using snappymeth.py
---

The output of snappymeth.py was post-processed in R (code coming soon), and then reads around this example SNP separated into separate bams by [splitSNP.py](https://github.com/astatham/splitSNP) then visualised with [IGV](https://www.broadinstitute.org/igv/).

![Example of ASM](http://i.imgur.com/8loeIYk.png)
