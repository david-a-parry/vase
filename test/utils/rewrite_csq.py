#!/usr/bin/env python3
import sys
import re
from parse_vcf import VcfReader

csq_fields = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene',
              'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON'
             ]
csq_format_re = re.compile(r'''.*Format:\s*((\S+\|)*\S+)"''')

transcripts = dict()
id_conversions = {'SYMBOL': dict(),
                  'Gene': dict(),
                  'Feature': dict()}
id_counts = {'SYMBOL': 0,
             'Gene': 0,
             'Feature': 0}

def rewrite_csq(rec):
    annots = []
    for c in rec.CSQ:
        for k, d in id_conversions.items():
            if c[k] in d:
                c[k] = d[c[k]]
            elif c[k]:
                id_counts[k] += 1
                n = "{}_{}".format(k.upper(), id_counts[k])
                id_conversions[k][c[k]] = n
                c[k] = n
        annots.append('|'.join(c[x] for x in c.keys() if x in csq_fields))
    return ','.join(annots)


if __name__ == '__main__':
    vcf = VcfReader(sys.argv[1])
    meta_data = []
    for x in vcf.header.meta_header:
        if x.startswith('##INFO=<ID=CSQ,'):
            d = vcf.header.metadata['INFO']['CSQ'][0]
            desc = d['Description']
            match = csq_format_re.match(desc)
            fields = match.group(1).split('|')
            new_f = [x for x in fields if x in csq_fields]
            new_desc = '"Consequence annotations from Ensembl VEP. Format: ' \
                       + '|'.join(new_f)
            s = str.join(',', ['##INFO' + "=<ID=" + 'CSQ'] +
                         ["Number=" + d['Number'],
                          "Type=" + d['Type'],
                          "Description=" + new_desc]) + "\">"
            meta_data.append(s)
        else:
            meta_data.append(x)
    vcf.header.meta_header = meta_data
    print(vcf.header, end='')
    for rec in vcf:
        csq = rewrite_csq(rec)
        rec.add_info_fields(info={'CSQ': csq})
        print(rec)
