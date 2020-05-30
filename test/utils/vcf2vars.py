#!/usr/bin/env python3
import pysam
import sys

if __name__ == '__main__':
    if len(sys.argv) > 2:
        sys.exit("Usage: {} [in.vcf]".format(sys.argv[0]))
    if len(sys.argv) > 1:
        vcf = pysam.VariantFile(sys.argv[1])
    else:
        vcf = pysam.VariantFile('-')
    pos_offset = 1000000
    for record in vcf:
        print("{}:{}-{}/{}".format(record.chrom,
                                   record.pos,
                                   record.ref,
                                   ','.join(record.alts)))
    vcf.close()
