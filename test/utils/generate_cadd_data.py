#!/usr/bin/env python3
import pysam
import sys
import random


def main(f):
    with pysam.VariantFile(f) as vcf:
        print('''##CADD style randomly generated numbers for testing
#Chrom  Pos     Ref     Alt     RawScore        PHRED''')
        for record in vcf:
            for alt in record.alts:
                if random.randint(0, 4):  # create score for 3/4 variants
                    raw = "{:f}".format(random.uniform(-2.5, 3.0))
                    phred = "{:g}".format(random.uniform(0, 40))
                    print("\t".join((record.chrom,
                                     str(record.pos),
                                     record.ref,
                                     alt,
                                     raw,
                                     phred)))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: {} ex1.vcf".format(sys.argv[0]))
    main(sys.argv[1])
