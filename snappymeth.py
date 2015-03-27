#! /usr/bin/env python
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

    def processSNP(chrom, pos, ref, alt, CpGs, cutoff_mapq=40, cutoff_baseq=30, min_coverage=10):
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
        # Only output rows with minimum coverage
        CpGs_bases[CpGs_bases.iloc[:, 5:].sum(1) >= min_coverage].to_csv(outfile, header=False, index=False)

    parser = argparse.ArgumentParser(description="snappymeth.py - "
        "Parses a paired bam and vcf and looks for regions of allele-specific methylation.")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("reference", help="Reference FASTA file")
    parser.add_argument("output", help="Filename for the output table (csv format)")
    parser.add_argument("--pair_distance", type=int, default=500, help="The distance in "
        "basepairs to search up and downstream from each SNP (default is 500)")
    parser.add_argument("--max_depth", type=int, default=1000000, help="Maximum number "
        "of reads to process at the specified SNP position (default is 1000000)")
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

    if not can_create_file(os.path.dirname(args.output)):
        print("Output %s is not writable!" % os.path.dirname(args.output))
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

    # Open the output file and write header
    outfile = open(args.output, "w")
    outfile.write("SNP.chr,SNP.pos,SNP.ref,SNP.alt,CpG.pos,ref.A,ref.C,ref.G,ref.T,ref.N,alt.A,alt.C,alt.G,alt.T,alt.N\n")
    for record in vcffile:
        call = record.samples[0]
        if call.is_het:
            n_ref = call['DP4'][0] + call['DP4'][1]
            n_alt = call['DP4'][2] + call['DP4'][3]
            if n_ref >= 10 and n_alt >= 10:
                CpGs = findCpGs(fafile, record.CHROM, record.POS, args.pair_distance)
                if len(CpGs) > 0:
                    processSNP(record.CHROM, record.POS, record.REF, record.ALT[0].sequence, CpGs)

    samfile.close()
    outfile.close()


if __name__ == '__main__':
    main()
