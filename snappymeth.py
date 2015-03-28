#! /usr/bin/env python
from __future__ import division

def main():
    import argparse
    import pysam
    import vcf
    from pyfasta import Fasta
    import os
    import tempfile
    import re
    import pandas
    from collections import OrderedDict
    from fisher import pvalue

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    def can_create_file(folder_path):
        try:
            tempfile.TemporaryFile(dir=folder_path)
            return True
        except OSError:
            return False

    def findCpGs(fafile, chrom, pos, distance):
        sequence = fafile[chrom][pos-distance:pos+distance]
        CpGs = [m.start() for m in re.finditer('CG', sequence)]
        return [x+pos-distance for x in CpGs]

    def processSNP(chrom, pos, ref, alt, CpGs, min_coverage, min_region_CpGs, cutoff_mapq, cutoff_baseq):
        # PASS 1 - find readnames of all reads with ref/alt alleles
        alleles = list()
        ref_readnames = set()
        alt_readnames = set()
        n_mapq = 0
        n_baseq = 0
        for pileup in samfile.pileup(chrom, pos-1, pos, max_depth=args.max_depth):
            if pileup.reference_pos == pos-1:  # filter for position of interest
                print "Processing %s reads covering SNP position %s:%s" % (
                    len(pileup.pileups), chrom, pos)
                for read in pileup.pileups:
                    # read mapping quality filter
                    if read.alignment.mapping_quality >= cutoff_mapq:
                        n_mapq += 1
                        SNP_base = read.alignment.query_sequence[read.query_position]
                        SNP_qual = read.alignment.query_qualities[read.query_position]
                        # base quality score filter @ SNP position
                        if SNP_qual >= cutoff_baseq:
                            n_baseq += 1
                            alleles.append(SNP_base)
                            if SNP_base == ref:
                                ref_readnames.add(read.alignment.query_name)
                            elif SNP_base == alt:
                                alt_readnames.add(read.alignment.query_name)
        print " - Found %s reads passing mapping quality filter of %s" % (n_mapq, cutoff_mapq)
        print " - Found %s reads passing base quality filter of %s" % (n_baseq, cutoff_baseq)
        print " - Found %s reads matching '%s' REF allele" % (len(ref_readnames), ref)
        print " - Found %s reads matching '%s' ALT allele" % (len(alt_readnames), alt)

        # Remove reads in both
        ref_and_alt = ref_readnames.intersection(alt_readnames)
        if len(ref_and_alt) > 0:
            print " - %s reads discarded for being ambiguous" % len(ref_and_alt)
            ref_readnames = ref_readnames.difference(ref_and_alt)
            alt_readnames = alt_readnames.difference(ref_and_alt)

        # PASS 2 - iterate through the CpG sites around the SNP, and count C/Ts in REF/ALT in reads
        CpGs_bases = pandas.DataFrame(OrderedDict([
            ('SNP.chr', chrom),
            ('SNP.pos', pos),
            ('SNP.ref', ref),
            ('SNP.alt', alt),
            ('CpG.pos', CpGs),
            ('ref.A', 0),
            ('ref.C', 0),
            ('ref.G', 0),
            ('ref.T', 0),
            ('ref.N', 0),
            ('alt.A', 0),
            ('alt.C', 0),
            ('alt.G', 0),
            ('alt.T', 0),
            ('alt.N', 0)]))
        for read in samfile.fetch(chrom, CpGs[0]-1, CpGs[len(CpGs)-1]+1):  # extra 1bp buffer
            # Is this a REF, ALT or neither read?
            if read.query_name in ref_readnames:
                read_type = 'ref.'
            elif read.query_name in alt_readnames:
                read_type = 'alt.'
            else:
                read_type = None
            if read_type is not None:
                # Work out where the methylation information is in the CpG site, and whether to complement it
                # Depends on read1/read2 and forward/reverse status
                if read.is_read1 and not read.is_reverse:  # First, forward
                    offset = 0
                    to_complement = False
                elif not read.is_read1 and read.is_reverse:  # Second, reverse
                    offset = 0
                    to_complement = False
                elif read.is_read1 and read.is_reverse:  # First, reverse
                    offset = 1
                    to_complement = True
                elif not read.is_read1 and not read.is_reverse:  # Second, forward
                    offset = 1
                    to_complement = True
                # Iterate through all aligned read positions, and store methylation calls
                for pair in read.get_aligned_pairs():
                    if pair[0] is not None and pair[1] is not None:
                        try:
                            i = CpGs.index(pair[1]-offset)
                            this_base = read.query_sequence[pair[0]]
                            if to_complement:
                                this_base = complement[this_base]
                            CpGs_bases.ix[i, read_type+this_base] += 1
                        except ValueError:
                            pass
        # Subset to rows with minimum coverage
        # Calculate coverage and methylation per CpG site
        CpGs_bases["ref.cov"] = [CpGs_bases.loc[i, ["ref.C", "ref.T"]].sum() for i in CpGs_bases.index]
        CpGs_bases["alt.cov"] = [CpGs_bases.loc[i, ["alt.C", "alt.T"]].sum() for i in CpGs_bases.index]
        CpGs_bases = CpGs_bases[CpGs_bases["ref.cov"] >= min_coverage]
        CpGs_bases = CpGs_bases[CpGs_bases["alt.cov"] >= min_coverage]
        if len(CpGs_bases.index)>0:  # If rows are left
            CpGs_bases["ref.meth"] = [CpGs_bases["ref.C"][i] / CpGs_bases["ref.cov"][i] for i in CpGs_bases.index]
            CpGs_bases["alt.meth"] = [CpGs_bases["alt.C"][i] / CpGs_bases["alt.cov"][i] for i in CpGs_bases.index]
            CpGs_bases["meth.diff"] = [CpGs_bases["ref.meth"][i] - CpGs_bases["alt.meth"][i] for i in CpGs_bases.index]

            # Calculate fisher pvalue per CpG site
            CpGs_bases["pvalue"] =  [pvalue(*CpGs_bases.loc[i, ["ref.C", "ref.T", "alt.C", "alt.T"]].tolist()).two_tail for i in CpGs_bases.index]

            # Export sites table
            CpGs_bases.to_csv(out_sites, header=False, index=False)
            
        if len(CpGs_bases.index) >= min_region_CpGs:  # If >=3 CpG sites, calculate fisher pvalue for pool for region and export
            output = "%s,%s,%s,%s,%s,%s,%s," % (
                chrom, pos, ref, alt,  # SNP position
                CpGs_bases["CpG.pos"].tolist()[0], CpGs_bases["CpG.pos"].tolist()[-1],  # First and last CpG sites of region
                len(CpGs_bases.index))  # Number of CpG sites in region

            # Sums of counts across the region
            CpGs_sums = CpGs_bases[["ref.C", "ref.T", "alt.C", "alt.T", "ref.cov", "alt.cov"]].sum(0).tolist()
            output += "%s,%s,%s,%s,%s,%s," % tuple(CpGs_sums)

            # Methylation ratios and pvalue
            ref_meth = CpGs_sums[0] / CpGs_sums[4]
            alt_meth = CpGs_sums[1] / CpGs_sums[5]
            meth_diff = ref_meth-alt_meth
            p_value = pvalue(*CpGs_sums[0:4]).two_tail
            output += "%s,%s,%s,%s\n" % (ref_meth, alt_meth, meth_diff, p_value)

            # Export row for this region
            out_regions.write(output)

    parser = argparse.ArgumentParser(description="snappymeth.py - "
        "Parses a paired bam and vcf and looks for regions of allele-specific methylation.")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("reference", help="Reference FASTA file")
    parser.add_argument("prefix", help="Prefix for the sites and regions output csvs")
    parser.add_argument("--pair_distance", type=int, default=500, help="The distance in "
        "basepairs to search up and downstream from each SNP (default is 500)")
    parser.add_argument("--max_depth", type=int, default=8000, help="Maximum number "
        "of reads to process at each SNP position (default is 8000)")
    parser.add_argument("--min_per_allele", type=int, default=5, help="Minimum number "
        "of reads containing each allele to process a SNP position")
    parser.add_argument("--min_sites_in_region", type=int, default=3, help="Minimum number "
        "of CpG sites linked to a SNP to perform a regional analysis")
    parser.add_argument("--min_mapping_quality", type=int, default=40, help="Minimum mapping "
        "quality score for a read to be considered")
    parser.add_argument("--min_base_quality", type=int, default=30, help="Minimum basecall "
        "quality score at the SNP for a read to be considered")
    args = parser.parse_args()

    if args.max_depth < 8000:
        print("Specified max_depth is too low - changing to 8000")
        args.max_depth = 8000

    # Check input files exists, and thet output folder is writeable
    if not os.path.isfile(args.input_vcf):
        print("Input VCF file %s does not exist!" % args.input_vcf)
        return

    if not os.path.isfile(args.input_bam):
        print("Input BAM file %s does not exist!" % args.input_bam)
        return

    if not os.path.isfile(args.reference):
        print("Reference FASTA file %s does not exist!" % args.reference)
        return

    if not can_create_file(os.path.dirname(args.prefix)):
        print("Output %s is not writable!" % os.path.dirname(args.prefix))
        return

    # Open the reference fasta file
    print "Loading %s" % args.reference
    fafile = Fasta(args.reference)

    # Index samfile if one does not already exist
    samfile = pysam.AlignmentFile(args.input_bam, "rb")
    if not samfile._hasIndex():
        print("BAM file '%s' does not have an index, creating one..." % args.input_bam)
        samfile.close()
        pysam.index(args.input_bam)
        samfile = pysam.AlignmentFile(args.input_bam, "rb")

    # Open the VCF file
    vcffile = vcf.Reader(filename=args.input_vcf, compressed=True)

    # Open the output files and write headers
    out_sites = open(args.prefix + ".sites.csv", "w")
    out_sites.write("SNP.chr,SNP.pos,SNP.ref,SNP.alt,CpG.pos,ref.A,ref.C,ref.G,ref.T,ref.N,"
        "alt.A,alt.C,alt.G,alt.T,alt.N,ref.cov,alt.cov,ref.meth,alt.meth,meth.diff,p.value\n")
    out_regions = open(args.prefix + ".regions.csv", "w")
    out_regions.write("SNP.chr,SNP.pos,SNP.ref,SNP.alt,first.CpG,last.CpG,nCG,ref.C,ref.T,alt.C,alt.T,"
        "ref.cov,alt.cov,ref.meth,alt.meth,meth.diff,p.val\n")

    # Iterate through the VCF
    for record in vcffile:
        call = record.samples[0]
        if call.is_het:
            n_ref = call['DP4'][0] + call['DP4'][1]
            n_alt = call['DP4'][2] + call['DP4'][3]
            if n_ref >= args.min_per_allele and n_alt >= args.min_per_allele:
                CpGs = findCpGs(fafile, record.CHROM, record.POS, args.pair_distance)
                if len(CpGs) > 0:
                    processSNP(record.CHROM, record.POS, record.REF, record.ALT[0].sequence, CpGs,
                        args.min_per_allele, args.min_sites_in_region, args.min_mapping_quality,
                        args.min_base_quality)

    samfile.close()
    out_sites.close()
    out_regions.close()


if __name__ == '__main__':
    main()
