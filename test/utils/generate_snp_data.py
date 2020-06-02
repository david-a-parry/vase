#!/usr/bin/env python3
import pysam
import sys
import random

dbsnp = "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz"
nts = ['A', 'C', 'G', 'T']


def get_var_list(f):
    vcf = pysam.VariantFile(f)
    return [(record.pos, record.ref, record.alts) for record in vcf]


def random_nt(alleles, prefix=''):
    random.shuffle(nts)
    for n in nts:
        if n not in alleles:
            return n
    return random_nt(alleles, prefix=n)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit("Usage: {} ex1.vcf snp.vcf <chrom:start-end>".format(
            sys.argv[0]))
    var_list = get_var_list(sys.argv[1])
    vcf = pysam.VariantFile(sys.argv[2])
    chrom, region = sys.argv[3].split(':')
    start, stop = region.split('-')
    out = pysam.VariantFile("-", mode='w', header=vcf.header)
    pos_offset = 1000000
    for record in vcf.fetch(chrom, int(start) - 1, int(stop)):
        record.chrom = '1'
        if var_list and random.randint(0, 1):
            variant = var_list.pop(random.randint(0, len(var_list) - 1))
            record.pos, record.ref = variant[:2]
            new_alts = []
            for i in range(len(record.alts)):
                if i >= len(variant[2]):
                    new_alts.append(random_nt((record.ref,) + variant[2]))
                else:
                    new_alts.append(variant[2][i])
            record.alts = new_alts
        else:
            record.pos = pos_offset + random.randint(10, 1e3)
            pos_offset = record.pos
        out.write(record)
    vcf.close()
    out.close()
