#!/usr/bin/env python3
import pysam
import sys
import random


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit("Usage: {} in.vcf sample1 [sample2...sampleN]".format(
            sys.argv[0]))
    vcf = pysam.VariantFile(sys.argv[1])
    samples = sys.argv[2:]
    vcf.subset_samples(samples)
    out = pysam.VariantFile("-", mode='w', header=vcf.header)
    pos_offset = 1000000
    try:
        first_chrom = vcf.index.keys()[0]
    except AttributeError:
        first_chrom = None
    for record in vcf:
        if first_chrom is None:
            if record.chrom.startswith('chr'):
                first_chrom = 'chr1'
            else:
                first_chrom =  '1'
        if any(y is not None and y > 0 for x in samples for y in
               record.samples[x].allele_indices):
            record.chrom = first_chrom
            record.pos = pos_offset + random.randint(10, 1e3)
            if random.randint(0, 10) % 2:
                record.id = "snp{}".format(random.randint(1, 1e5))
            else:
                record.id = '.'
            out.write(record)
            pos_offset = record.pos
    vcf.close()
    out.close()
