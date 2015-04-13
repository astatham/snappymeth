#! /usr/bin/env python
from __future__ import division
from __future__ import print_function

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
    import sys
    import gzip
    import csv
    from IGV import IGV
    from multiprocessing import Process, Queue

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    def can_create_file(folder_path):
        try:
            tempfile.TemporaryFile(dir=folder_path)
            return True
        except OSError:
            return False

    def findCpGs(fafile, chrom, pos, distance):
        minpos = 0 if pos<distance else pos-distance
        sequence = fafile[chrom][minpos:pos+distance]
        CpGs = [m.start() for m in re.finditer('CG', sequence, flags=re.IGNORECASE)]
        return [x+minpos for x in CpGs]

    def type_of_read(read):
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
        return offset, to_complement
    
    def IGV_reader(queue):
        ## Read from the queue
        while True:
            msg = queue.get()         # Read from the queue and do nothing
            if (msg == 'DONE'):
                break
            chrom, pos, ref, alt, ref_filename, alt_filename = msg.split(",")
            pos = int(pos)
            png_filename = os.path.basename("%s.%s.%s.%s.%s.png" % (args.prefix, chrom, pos, ref, alt))
            igv.load("file://"+os.path.abspath(ref_filename))
            igv.load("file://"+os.path.abspath(alt_filename))
            igv.go("%s:%s-%s" % (chrom, pos-250, pos+250))
            igv.send("collapse")
            igv.send("region %s %s %s" % (chrom, pos+1, pos+2))  # 1 based co-ordinates for IGV
            igv.save(png_filename)
            igv.clear()

    def exportBAMs(chrom, pos, ref, alt, minpos, maxpos, ref_readnames, alt_readnames):
        ref_filename = "%s.%s.%s.ref.%s.bam" % (args.prefix, chrom, pos, ref)
        ref_bam = pysam.AlignmentFile(ref_filename, "wb", template=samfile)
        alt_filename = "%s.%s.%s.alt.%s.bam" % (args.prefix, chrom, pos, alt)
        alt_bam = pysam.AlignmentFile(alt_filename, "wb", template=samfile)

        for read in samfile.fetch(chrom, minpos, maxpos):
            if read.query_name in ref_readnames:
                ref_bam.write(read)
            elif read.query_name in alt_readnames:
                alt_bam.write(read)

        ref_bam.close()
        alt_bam.close()
        pysam.index(ref_filename)
        pysam.index(alt_filename)

        if args.IGV_screenshot:
            IGV_queue.put("%s,%s,%s,%s,%s,%s" % (chrom, pos, ref, alt, ref_filename, alt_filename))

    def processReadsAtPosition(chrom, pos, ref, alt, CpGs, ref_readnames, alt_readnames, min_coverage, 
        min_region_CpGs):
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
                offset, to_complement = type_of_read(read)
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

            # Export BAM per allele if feature is turned on and region meets fisher_cutoff
            if args.region_bams and p_value <= args.fisher_cutoff:
                print(" - Regional fisher exact p_value: %s - exporting BAMs" % p_value)
                exportBAMs(chrom, pos, ref, alt, CpGs[0]-1, CpGs[len(CpGs)-1]+1,
                    ref_readnames, alt_readnames)

    def processCpG(chrom, pos, cutoff_mapq, cutoff_baseq):
        """
        Find readnames of all reads that are meth or unmeth at the specified CpG position
        """
        M_readnames = set()
        U_readnames = set()
        n_mapq = 0
        n_baseq = 0
        for pileup in samfile.pileup(chrom, pos, pos+1):
            if pileup.reference_pos == pos:  # filter for position of interest
                print("Processing %s reads covering CpG position %s:%s" % (
                    len(pileup.pileups), chrom, pos))
                for read in pileup.pileups:
                    # read mapping quality filter
                    if read.alignment.mapping_quality >= cutoff_mapq:
                        n_mapq += 1
                        offset, to_complement = type_of_read(read.alignment)
                        if read.query_position + offset < len(read.alignment.query_sequence):
                            CpG_base = read.alignment.query_sequence[read.query_position + offset]
                            if to_complement:
                                CpG_base = complement[CpG_base]
                            CpG_qual = read.alignment.query_qualities[read.query_position + offset]
                            # base quality score filter @ SNP position
                            if CpG_qual >= cutoff_baseq:
                                n_baseq += 1
                                if CpG_base == "C":
                                    M_readnames.add(read.alignment.query_name)
                                elif CpG_base == "T":
                                    U_readnames.add(read.alignment.query_name)
        print(" - Found %s reads passing mapping quality filter of %s" % (n_mapq, cutoff_mapq))
        print(" - Found %s reads passing base quality filter of %s" % (n_baseq, cutoff_baseq))
        print(" - Found %s reads with M allele" % len(M_readnames))
        print(" - Found %s reads with U allele" % len(U_readnames))

        # Remove reads in both
        M_and_U = M_readnames.intersection(U_readnames)
        if len(M_and_U) > 0:
            print(" - %s reads discarded for being ambiguous" % len(M_and_U))
            M_readnames = M_readnames.difference(M_and_U)
            U_readnames = U_readnames.difference(M_and_U)
        return M_readnames, U_readnames

    def processSNP(chrom, pos, ref, alt, cutoff_mapq, cutoff_baseq):
        """
        Find readnames of all reads with REF and ALT alleles
        """
        ref_readnames = set()
        alt_readnames = set()
        n_mapq = 0
        n_baseq = 0
        for pileup in samfile.pileup(chrom, pos, pos+1):
            if pileup.reference_pos == pos:  # filter for position of interest
                print("Processing %s reads covering SNP position %s:%s" % (
                    len(pileup.pileups), chrom, pos))
                for read in pileup.pileups:
                    # read mapping quality filter
                    if read.alignment.mapping_quality >= cutoff_mapq:
                        n_mapq += 1
                        SNP_base = read.alignment.query_sequence[read.query_position]
                        SNP_qual = read.alignment.query_qualities[read.query_position]
                        # base quality score filter @ SNP position
                        if SNP_qual >= cutoff_baseq:
                            n_baseq += 1
                            if SNP_base == ref:
                                ref_readnames.add(read.alignment.query_name)
                            elif SNP_base == alt:
                                alt_readnames.add(read.alignment.query_name)
        print(" - Found %s reads passing mapping quality filter of %s" % (n_mapq, cutoff_mapq))
        print(" - Found %s reads passing base quality filter of %s" % (n_baseq, cutoff_baseq))
        print(" - Found %s reads matching '%s' REF allele" % (len(ref_readnames), ref))
        print(" - Found %s reads matching '%s' ALT allele" % (len(alt_readnames), alt))

        # Remove reads in both
        ref_and_alt = ref_readnames.intersection(alt_readnames)
        if len(ref_and_alt) > 0:
            print(" - %s reads discarded for being ambiguous" % len(ref_and_alt))
            ref_readnames = ref_readnames.difference(ref_and_alt)
            alt_readnames = alt_readnames.difference(ref_and_alt)
        return ref_readnames, alt_readnames

    ## Entry point
    parser = argparse.ArgumentParser(description="snappymeth.py - "
        "Discover sites and regions of allele specific methylation from whole genome bisulfite "
        "sequencing data by counting CpG methylation on alleles separately. Reads can be "
        "separated by either a heterozygous SNP (when a VCF is supplied), or by the methylation "
        "status of a single CpG site. Both analyses modes require sufficient sequencing coverage "
        "of both alleles (default is 10x).")

    parser.add_argument("input_file", help="Input VCF/CpG sites file, gzip compressed." )
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("reference", help="Reference FASTA file")
    parser.add_argument("prefix", help="Prefix for all output files - the sites and regions output csvs, "
        "regional BAMs and IGV screenshots")

    parser.add_argument("--input_type", choices=("VCF", "CpGs"), default="VCF", help="Whether the "
        "input_file is a VCF (default) or a csv of methylation counts at CpG sites with the format "
        "'chr,position,M,U' where the fields are chromosome name, 0-based position of the CpG site, "
        "count of methylated bases sequenced at this site and count of unmethylated bases sequenced.")
    parser.add_argument("--VCF_sample", default="0", help="The sample in the VCF to be processed - "
        "either as the sample name or numeric index (0-based). Default is 0, the first sample.")
    parser.add_argument("--pair_distance", type=int, default=500, help="The distance in "
        "basepairs to search up and downstream from each position (default is 500).")
    parser.add_argument("--max_depth", type=int, default=100, help="Maximum number "
        "of reads allowed at a position to try and filter out repeat reads (default is 100)..")
    parser.add_argument("--min_per_allele", type=int, default=5, help="Minimum number "
        "of reads containing each allele to process a position.")
    parser.add_argument("--min_sites_in_region", type=int, default=3, help="Minimum number "
        "of CpG sites linked to a SNP to perform a regional analysis.")
    parser.add_argument("--min_mapping_quality", type=int, default=40, help="Minimum mapping "
        "quality score for a read to be considered.")
    parser.add_argument("--min_base_quality", type=int, default=30, help="Minimum basecall "
        "quality score at the SNP for a read to be considered.")
    parser.add_argument("--region_bams", default=False, action='store_true', help="Specity to output "
        "BAM files per allele when the regional fisher exact p-value is less than the cutoff "
        "specified by --fisher_cutoff.")
    parser.add_argument("--fisher_cutoff", type=float, default=0.0001, help="Regional fisher exact "
        "p-value cutoff for a regional BAM to be created/IGV screenshot be taken (default is 0.0001).")
    parser.add_argument("--IGV_screenshot", default=False, action='store_true', help="Specity to take "
        "IGV screenshots of each region that passes --fisher_cutoff. Requires that IGV be running on "
        "the local machine and listening on port 60151")
    args = parser.parse_args()

    # Check input files exists, and thet output folder is writeable
    if not os.path.isfile(args.input_file):
        print("Input file %s does not exist!" % args.input_file)
        return

    if not os.path.isfile(args.input_bam):
        print("Input BAM file %s does not exist!" % args.input_bam)
        return

    if not os.path.isfile(args.reference):
        print("Reference FASTA file %s does not exist!" % args.reference)
        return

    if not can_create_file(os.path.dirname(args.prefix)):
        print("Output directory %s/ is not writable!" % os.path.dirname(args.prefix))
        return

    # Setup for IGV
    if args.IGV_screenshot:
        args.region_bams = True
        igv = IGV()
        igv.clear()
        print("BAMs and IGV screenshots will be saved in %s" % os.path.dirname(os.path.abspath(args.prefix)))
        igv.set_path(os.path.dirname(os.path.abspath(args.prefix)))
        # Setup queue for IGV screenshots in separate process
        print("Starting separate process for IGV screenshots")
        IGV_queue = Queue()
        reader_process = Process(target=IGV_reader, args=((IGV_queue),))
        reader_process.daemon = True
        reader_process.start()        # Launch IGV_reader() as a separate python process

    # Open the reference fasta file
    print("Loading %s" % args.reference)
    fafile = Fasta(args.reference)

    # Index samfile if one does not already exist
    samfile = pysam.AlignmentFile(args.input_bam, "rb")
    if not samfile._hasIndex():
        print("BAM file '%s' does not have an index, creating one..." % args.input_bam)
        samfile.close()
        pysam.index(args.input_bam)
        samfile = pysam.AlignmentFile(args.input_bam, "rb")

    # Open the output files and write headers
    out_sites = open(args.prefix + ".sites.csv", "w")
    out_sites.write("SNP.chr,SNP.pos,SNP.ref,SNP.alt,CpG.pos,ref.A,ref.C,ref.G,ref.T,ref.N,"
        "alt.A,alt.C,alt.G,alt.T,alt.N,ref.cov,alt.cov,ref.meth,alt.meth,meth.diff,p.value\n")
    out_regions = open(args.prefix + ".regions.csv", "w")
    out_regions.write("SNP.chr,SNP.pos,SNP.ref,SNP.alt,first.CpG,last.CpG,nCG,ref.C,ref.T,alt.C,alt.T,"
        "ref.cov,alt.cov,ref.meth,alt.meth,meth.diff,p.val\n")

    if args.input_type=="VCF":  # VCF analysis
        # Open the VCF file
        vcffile = vcf.Reader(filename=args.input_file, compressed=True)

        # Check VCF_sample validity
        if args.VCF_sample.isdigit():  # If a number convert to int
            args.VCF_sample = int(args.VCF_sample)
        if isinstance(args.VCF_sample, basestring):
            try:
                sample_no = vcffile.samples.index(args.VCF_sample)
            except ValueError:
                sys.exit("Sample %s not found in VCF!" % args.VCF_sample)
        elif not args.VCF_sample < len(vcffile.samples):
                sys.exit("Sample number %s not found in VCF!" % args.VCF_sample)
        else:
            sample_no = args.VCF_sample

        print("Processing sample no %s (%s) from VCF" % (sample_no, vcffile.samples[sample_no]))

        # Iterate through the VCF
        for record in vcffile:
            call = record.samples[sample_no]
            if call.is_het:
                n_ref = call['DP4'][0] + call['DP4'][1]
                n_alt = call['DP4'][2] + call['DP4'][3]
                if n_ref >= args.min_per_allele and n_alt >= args.min_per_allele and (n_ref + n_alt) <= args.max_depth:
                    # record.POS-1 as VCFs are 1 based and everything is 0 based
                    CpGs = findCpGs(fafile, record.CHROM, record.POS-1, args.pair_distance)
                    # If SNP overlaps a CpG site, remove
                    for site in range(record.POS-2, record.POS+1):
                        if site in CpGs:
                            CpGs.remove(site)
                    if len(CpGs) > 0:  # If there are any CpG sites in the vicinity
                        ref_reads, alt_reads = processSNP(record.CHROM, record.POS-1, record.REF,
                            record.ALT[0].sequence, args.min_mapping_quality, args.min_base_quality)
                        if len(ref_reads) + len(alt_reads) <= args.max_depth:
                            processReadsAtPosition(record.CHROM, record.POS-1, record.REF,
                                record.ALT[0].sequence, CpGs, ref_reads, alt_reads, args.min_per_allele,
                                args.min_sites_in_region)

    else:  ## CpG sites analysis
        with gzip.open(args.input_file, "r") as f:
            CpGreader = csv.DictReader(f)
            if CpGreader.fieldnames != ['chr', 'position', 'M', 'U']:
                sys.exit("Field names in %s must be 'chr,position,M,U'" % args.input_file)
            for CpG in CpGreader:
                if int(CpG["M"]) >= args.min_per_allele and int(CpG["U"]) >= args.min_per_allele and (int(CpG["M"]) + int(CpG["U"])) <= args.max_depth:
                    CpGs = findCpGs(fafile, CpG["chr"], int(CpG["position"]), args.pair_distance)
                    try:
                        CpGs.remove(int(CpG["position"]))  # Remove the CpG site we are processing
                    except ValueError:
                        sys.exit("Input file CpG site at '%s:%s' is a '%s' in reference. Are you sure your input file coordinates are 0-based?" % (CpG["chr"], CpG["position"], fafile[CpG["chr"]][int(CpG["position"]):int(CpG["position"])+2]))
                    if len(CpGs) > 0:  # If there are any other CpG sites in the vicinity
                        M_reads, U_reads = processCpG(CpG["chr"], int(CpG["position"]),
                            args.min_mapping_quality, args.min_base_quality)
                        if len(M_reads) + len(U_reads) <= args.max_depth:
                            processReadsAtPosition(CpG["chr"], int(CpG["position"]), "M", "U", CpGs,
                                M_reads, U_reads, args.min_per_allele, args.min_sites_in_region)

    # Close down IGV process
    if args.IGV_screenshot:
        IGV_queue.put("DONE")
        print("Waiting for IGV screenshots process to finish")
        reader_process.join()

    # close files
    samfile.close()
    out_sites.close()
    out_regions.close()

if __name__ == '__main__':
    main()
