#!/usr/bin/env python3
import pysam
import sys
import random
from vase.vcf_reader import VcfReader

header = '''##fileformat=VCFv4.2
##fileDate=20191004
##reference=GRCh38/hg38
##contig=<ID=1,length=248956422>
##contig=<ID=2,length=242193529>
##contig=<ID=3,length=198295559>
##contig=<ID=4,length=190214555>
##contig=<ID=5,length=181538259>
##contig=<ID=6,length=170805979>
##contig=<ID=7,length=159345973>
##contig=<ID=8,length=145138636>
##contig=<ID=9,length=138394717>
##contig=<ID=10,length=133797422>
##contig=<ID=11,length=135086622>
##contig=<ID=12,length=133275309>
##contig=<ID=13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=15,length=101991189>
##contig=<ID=16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=18,length=80373285>
##contig=<ID=19,length=58617616>
##contig=<ID=20,length=64444167>
##contig=<ID=21,length=46709983>
##contig=<ID=22,length=50818468>
##contig=<ID=X,length=156040895>
##contig=<ID=Y,length=57227415>
##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3 variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'''


def create_spliceai_score(symbol, alt):
    scores = []
    dists = []
    for i in range(4):
        scores.append("{:.2f}".format(abs(random.gauss(0, 0.2))))
        dists.append("{:d}".format(int(random.gauss(1, 50))))
    return "{}|{}|".format(alt, symbol,) + '|'.join(scores + dists)


def main(f):
    with VcfReader(f) as vcf:
        vcf.header.header.info.add(
            "SpliceAI",
            ".",
            "String",
            '''SpliceAIv1.3 variant annotation. These include delta scores \
(DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), \
donor gain (DG), and donor loss (DL). Format: \
ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL''')
        out = pysam.VariantFile("-", header=vcf.header.header, mode='w')
        for record in vcf:
            i = 0
            scores = []
            for alt in record.alts:
                i += 1
                if not random.randint(0, 4):  # create score for 3/4 variants
                    continue
                for symbol in set(x['SYMBOL'] for x in record.CSQ if
                                  x['alt_index'] == i and
                                  x['BIOTYPE'] == 'protein_coding'
                                  and 'stream' not in x['Consequence']):
                    scores.append(create_spliceai_score(symbol, alt))
            if scores:
                for f in (x for x in record.info if x != 'SpliceAI'):
                    del record.info[f]
                record.id = None
                record.record.qual = None
                record.filter.clear()
                record.info['SpliceAI'] = scores
                out.write(record.record)
    out.close()


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: {} ex1.vcf".format(sys.argv[0]))
    main(sys.argv[1])
